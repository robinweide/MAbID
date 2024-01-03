#!/usr/bin/env python3

"""
Based on Franka's script:
FR181122.process_reads_and_count_gatcs_and_flatten_umis.py
- adds RE_agnosticity 
"""

import os
import sys
import logging
from collections import defaultdict, Counter, OrderedDict
from itertools import groupby, chain, product
from math import floor

import numpy as np
import h5py
import pysam

log = logging.getLogger(__name__)
logsh = logging.StreamHandler()
logsh.setLevel(logging.DEBUG)
logfmt = logging.Formatter('[%(asctime)s] [%(levelname)s] %(message)s')
logsh.setFormatter(logfmt)
log.addHandler(logsh)

DEFAULT_MIN_MAPQ = 0
DEFAULT_ELEMENT_LENGTH = 4

BASES = 'NACGT'
CHUNKSIZE_MAX = 256 * 1024

class InvalidOutfile():
    def __init__(self, fn):
        self.fn = fn

    def write(self, read):
        raise NotImplementedError

class InvalidReadBAMWriter(InvalidOutfile):
    def __init__(self, fn, inbf, *args, **kwargs):
        super().__init__(fn)
        self._bamfh = pysam.Samfile(fn, 'wb', header=inbf.header)

    def write(self, read):
        self._bamfh.write(read)

    def __del__(self):
        self._bamfh.close()

class InvalidReadTextWriter(InvalidOutfile):
    def __init__(self, fn, *args, **kwargs):
        super().__init__(fn)
        self._fh = open(fn, 'wt')

    def write(self, read):
        if read.is_secondary:
            return

        self._fh.write(read.qname)
        self._fh.write("\n")

    def __del__(self):
        self._fh.close()

class InvalidReadFastqWriter(InvalidOutfile):
    def __init__(self, fn, *args, **kwargs):
        super().__init__(fn)
        self._fh = open(fn, 'wt')

    def write(self, read):
        if read.is_secondary:
            return

        self._fh.write("@")
        self._fh.write(read.query_name)
        self._fh.write("\n")
        self._fh.write(read.query_sequence)
        self._fh.write("\n+\n")
        self._fh.write("".join(chr(x + 33) for x in read.query_qualities))
        self._fh.write("\n")

    def __del__(self):
        self._fh.close()

_EXT2WRITER = {".bam": InvalidReadBAMWriter, ".fastq": InvalidReadFastqWriter}

# function: generate dictionary for conversion of pos to index
def mappable_pos_file_to_posidx(fn):
    h5file = h5py.File(fn, 'r')
    mappable_pos = dict()
    for key in h5file:
        mappable_pos[key] = h5file[key][:].cumsum()
    return mappable_pos, h5file.attrs["element_length"]

# function: retrieve UMI seq from readname
def get_umi_seq(s):
    readname_fields = s.split(":")
    assert readname_fields[-2] == "UMI"
    return readname_fields[-1]

# function: get writer for invalid alignments
def get_invalid_outfh(fn):
    ext = os.path.splitext(fn)[1].lower()
    writercls = _EXT2WRITER.get(ext, InvalidReadTextWriter)
    return writercls

# function: calculate the hamming distance between two sequences
def hamming_dist(s1, s2):
    dist = abs(len(s1) - len(s2))
    for char1, char2 in zip(s1, s2):
        if char1 != char2:
            dist += 1
    return dist

# function: pick the top KEEPN umis with differing at least MINHAMDISTANCE and flatten counts
def reduce_umi_sequence_derivatives(counts, keepn, min_editdistance):
    reduced_counts = {}
    for k in counts:
        k_counts = counts[k].copy()
        picked = set()
        unpicked = set(k_counts.keys())
        while (unpicked) and (len(picked) < keepn):
            seq = sorted(sorted(k_counts.most_common()), key=lambda x: -x[1])[0][0]
            unpicked.remove(seq)
            del k_counts[seq]
            picked.add(seq)
            derivs = set()
            for seq2 in unpicked:
                dist = hamming_dist(seq, seq2)
                assert dist > 0
                if dist <= min_editdistance:
                    derivs.add(seq2)

            for seq2 in derivs:
                unpicked.remove(seq2)
                del k_counts[seq2]
        reduced_counts[k] = len(picked)

    return reduced_counts

# function: write count array to hdf5 file
def write_counts(outfn, chrom, counts):
    assert counts.shape[1] == 2, "assuming 2nd dimension == strand"

    itemsperrow = np.multiply.reduce(counts.shape[1:])
    itemsize = counts.itemsize
    chunkrows = max(1, floor(CHUNKSIZE_MAX / itemsperrow / itemsize))
    chunks = (min(counts.shape[0], chunkrows), ) + tuple(counts.shape[1:])

    with h5py.File(outfn, 'a') as f:
        f.create_dataset(chrom, data=counts, chunks=chunks, compression=3)
        f.flush()
    return

def count_reads_with_umis(bamfn, outfn, mappable_pos, min_mapq, element_length, umi_length, keepn, min_editdistance, offset, invalid_outfn=None):
    # open the BAM files
    bf = pysam.Samfile(bamfn, 'rb')

    # assert that all chromosomes in alignment file are also contained in mappable_pos
    bfchroms = set(bf.references)
    assert all(chrom in bfchroms for chrom in mappable_pos)

    # keep track of missing chromosomes by discounting seen chroms
    chroms = sorted(mappable_pos.keys())
    missing_chroms = set(chroms)

    # iniatie read & count stats
    if umi_length is not None:
        stats = OrderedDict([("total_reads", 0), ("unmapped_reads", 0), ("lowmapq_reads", 0), ("non_ligation_reads", 0), ("invalidumi_reads", 0), ("valid_reads", 0), ("unique_counts", 0)])
    else:
        stats = OrderedDict([("total_reads", 0), ("unmapped_reads", 0), ("lowmapq_reads", 0), ("non_ligation_reads", 0), ("valid_reads", 0), ("counts", 0)])

    i = 0 # counter
    for chrom, readiter_chrom in groupby(bf, key=lambda read: (read.reference_name if not read.is_unmapped else None)):
        # chroms not in mappable pos are "None", corresponding to unaligned reads
        # write these reads to invalid read file
        if chrom not in mappable_pos:
            log.warning("Skipping chrom: %s" % chrom)
            for read in readiter_chrom:
#                invalid_outfh.write(read)
                if not read.is_secondary:
                    stats["total_reads"] += 1
                    stats["unmapped_reads"] += 1
            continue

        missing_chroms.remove(chrom)

        # initiate a counter dictionary that keeps track of observed ligation positions, their strand, the observed UMIs and their count
        # if umi_length is not specified, UMIs are not taken into account
        if umi_length is not None:
            umi_counts = {s: defaultdict(lambda: Counter()) for s in ["+", "-"]}
        else:
            umi_counts = {s: Counter() for s in ["+", "-"]}

        # loop over all reads aligning to chrom
        for read in readiter_chrom:
            i += 1
            if i % 1000000 == 0:
                log.info("At alignment: %d" % i)

            if not read.is_secondary:
                stats["total_reads"] += 1

            # discard reads that are unaligned, secondary or that have a too low mapping quality
            if (read.mapq < min_mapq) or read.is_secondary or read.is_unmapped:
                #invalid_outfh.write(read)
                if not read.is_secondary:
                    if read.mapq < min_mapq:
                        stats["lowmapq_reads"] += 1
                    elif read.is_unmapped:
                        stats["unmapped_reads"] += 1
                continue

            if not read.is_reverse:
                ligationpos = read.pos + offset
                strand = "+"
            else:
                ligationpos = read.reference_end - (element_length + offset)
                strand = "-"

            # discard reads that do not align to a ligation
            if ligationpos not in mappable_pos[read.reference_name]:
                stats["non_ligation_reads"] += 1
                continue

            # add read to umi_count dict
            if umi_length is not None:
                umiseq = get_umi_seq(read.qname)
                # discard reads that have UMIs containing invalid bases
                if not all([base in BASES for base in umiseq]):
                    stats["invalidumi_reads"] += 1
                    continue
                else:
                    stats["valid_reads"] += 1
                    umi_counts[strand][ligationpos][umiseq] += 1
            else:
                stats["valid_reads"] += 1
                umi_counts[strand][ligationpos] += 1

        # reducing chromosome counts
        if umi_length is not None:
            counts = {s: reduce_umi_sequence_derivatives(umi_counts[s], keepn, min_editdistance) for s in ["+", "-"]}
            d = np.array([[counts[s].get(i, 0) for i in mappable_pos[chrom]] for s in ["+", "-"]])
            stats["unique_counts"] += d.sum()
        else:
            d = np.array([[umi_counts[s].get(i, 0) for i in mappable_pos[chrom]] for s in ["+", "-"]])
            stats["counts"] += d.sum()

        write_counts(outfn, chrom, np.transpose(d))

    for chrom in missing_chroms:
        zeros = np.zeros((len(mappable_pos[chrom]), 2), dtype="uint32")
        write_counts(outfn, chrom, zeros)

    return stats


def main():
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--verbose', '-v', action='count', default=0, required=False)
    ap.add_argument('--quiet', '-q', action="store_true", default=False, required=False)
    ap.add_argument('--offset', help="Set to 0 for GATC; -1 (default) for TTAA", default=-1, type=int, required=False)
    ap.add_argument('--outfile', help="HDF5 file to write counts to", required=True)
    ap.add_argument('--invalid-outfile', help="BAM file to write non-ligation reads to", required=False)
    ap.add_argument('--min-mapq', type=int, default=DEFAULT_MIN_MAPQ, required=False)
    ap.add_argument('--umi-length', type=int, default=None, required=False)
    ap.add_argument('--keep-n', help='How many distinct UMIs to keep. Default=all', type=int, required=False, default=None)
    ap.add_argument('--min-editdistance', type=int, required=False, default=1)
    ap.add_argument('--pos-file', required=True, metavar="HDF5_FILE", help="File with motif positions")
    ap.add_argument('infile', help="BAM file to be filtered. Use '-' for STDIN")

    args = ap.parse_args()

    if args.quiet:
        log.setLevel(logging.CRITICAL)
    else:
        print(args.verbose)
        log.setLevel([logging.WARNING, logging.INFO, logging.DEBUG][min(2, args.verbose)])
        #log.setLevel([logging.WARNING, logging.INFO, logging.DEBUG])

    for fn in (args.outfile, args.invalid_outfile):
        if fn is not None and os.path.exists(fn):
            ap.error("Output file exists (%s)" % fn)

    log.debug("Reading mappable pos file %s" % args.pos_file)
    mappable_pos, element_length = mappable_pos_file_to_posidx(args.pos_file)
    log.debug("Reading mappable pos file %s [DONE]" % args.pos_file)

    log.debug("Counting reads and filtering from file %s" % args.infile)

    stats = count_reads_with_umis(args.infile, args.outfile, mappable_pos, args.min_mapq, element_length, args.umi_length, args.keep_n, args.min_editdistance, args.offset, args.invalid_outfile)
    log.info("; ".join(map(lambda t: "%s: %d" % t, stats.items())))
    log.debug("Counting reads and filtering from file %s [DONE]" % args.infile)

    return

if __name__ == "__main__":
    main()


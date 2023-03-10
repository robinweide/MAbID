#!/hpc/hub_kind/rvanderweide/bin/miniconda3/bin/python

"""
"""

import os
from math import floor
import logging

import numpy as np
import h5py

from posarray import read_posfile


CHUNKSIZE_MAX = 256 * 1024


log = logging.getLogger(__name__)
logsh = logging.StreamHandler()
logsh.setLevel(logging.DEBUG)
logfmt = logging.Formatter('[%(asctime)s] [%(levelname)s] %(message)s')
logsh.setFormatter(logfmt)
log.addHandler(logsh)


def bin_countfiles(infiles, pos, binsize, mapab=None, keepn=None):
    assert any(infiles)
    fs = [h5py.File(fn, 'r') for fn in infiles]

    chroms = sorted(pos.keys())

    binned_pos = {chrom: pos[chrom] // binsize for chrom in chroms}
    binned_chromsizes = {chrom: int(binned_pos[chrom][-1]) + 1 for chrom in chroms}
    binned_counts = {chrom: np.zeros(binned_chromsizes[chrom], dtype=int) for chrom in chroms}
    for chrom in pos:
        log.debug("In chrom '%s'" % chrom)

        for f in fs:
            assert chrom in f

            d = f[chrom][:]
            if d.ndim == 3:
                # (positions, strands, UMIs) --> each value is a PCR dup count
                # this is a "counts" file, not a "unique_events" counts file...
                # this makes sense for bulk/population data
                # sum over all PCR duplicate counts
                d = d.sum(axis=2)
            else:
                assert d.ndim == 2  # (positions, strands)

            if keepn is not None:
                # keep at most `keepn` counts per position
                d = np.minimum(d, keepn)

            if mapab is not None:
                # set counts at unmappable positions to 0
                if chrom in mapab:
                    d = np.choose(
                        mapab[chrom][:] > 0,
                        (0, d),
                    )
                else:
                    # happens if the split-allele mappability file has no
                    # (unique) alignments on this chromosome (or contigs/ERCCs)
                    d[:] = 0

            d = d.sum(axis=1)  # sum across strands

            a = binned_counts[chrom]
            locs = binned_pos[chrom]

            np.add.at(a, locs, d)

    return binned_counts


def write_counts(outfn, counts):
    h5file = h5py.File(outfn, 'w')

    itemsize = next(iter(counts.values())).itemsize
    shape = next(iter(counts.values())).shape
    assert len(shape) == 1
    bytesperpos = itemsize
    chunkrows = max(1, floor(CHUNKSIZE_MAX / bytesperpos))

    for chrom in sorted(counts):
        rows = counts[chrom].shape[0]
        chunks = (min(rows, chunkrows), )
        h5file.create_dataset(chrom, data=counts[chrom], chunks=chunks, compression=3)
        h5file.flush()
    h5file.close()

    return


def main():
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--verbose', '-v', action='count', default=0, required=False)
    ap.add_argument('--quiet', '-q', action="store_true", default=False, required=False)
    ap.add_argument('--keepn', '-n', type=int, default=None, required=False, help="Keep at most `KEEPN` unique events per position. All if not given.")
    ap.add_argument('--outfile', help="HDF5 file to write counts to", required=True)
    ap.add_argument('--binsize', required=False, metavar="INT", type=int, default=2000)
    ap.add_argument('--posfile', help="File with GATC positions", required=True, metavar="HDF5_FILE")
    ap.add_argument('--mapfile', help='Mappability file', required=False, metavar="HDF5_FILE")
    ap.add_argument('infiles', metavar="FILE", nargs="+", help="Count file with number of unique events")

    args = ap.parse_args()

    if args.quiet:
        log.setLevel(logging.CRITICAL)
    else:
        log.setLevel([logging.WARNING, logging.INFO, logging.DEBUG][min(2, args.verbose)])

    if args.binsize <= 0:
        ap.error("Invalid binsize: %d" % args.binsize)

    if os.path.exists(args.outfile):
        ap.error("Output file exists (%s)" % args.outfile)

    for fn in args.infiles + [args.posfile]:
        if not os.access(fn, os.R_OK):
            ap.error("Can't read file: %s" % fn)

    if args.mapfile is not None:
        assert os.access(args.mapfile, os.R_OK)
        mapab = h5py.File(args.mapfile, 'r')
    else:
        mapab = None

    pos = read_posfile(args.posfile)
    binned_counts = bin_countfiles(args.infiles, pos, args.binsize, mapab, args.keepn)

    write_counts(args.outfile, binned_counts)

    return


if __name__ == "__main__":
    main()

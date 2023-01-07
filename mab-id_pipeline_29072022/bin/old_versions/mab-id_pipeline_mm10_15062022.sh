#!/bin/bash 

# MAb-ID simple pipeline June 2022
#
# This is a simple -but working— pipeline to map and parse MAb-ID data. It is in
# no way ready for production, so please contact Robin H. van der Weide before
# usage via e-mail: r.weide[at]hubrecht.eu.
#
# Expectations:
# > mouse data for mm10
# > single-end 100bp
# > one illumina-barcode

params_runid=index_07 # NEVER USE A DOT IN THIS !
# the sheet must contain the following columns:
# | illumina_index	ABBC_barcode	SBC_barcode	RS
# | index_07	      ABBC_001	    SBC_001	    TTAA
# | index_07	      ABBC_001	    SBC_002	    TTAA
params_sheet=index7.csv

params_abbc_file=/home/hub_kind/rvanderweide/projects/MAb-ID/current_pipeline/current_adapters/MAbID_ABBC_1-7_11-15.tsv
params_sbc_file=/home/hub_kind/rvanderweide/projects/MAb-ID/current_pipeline/current_adapters/scAbID_SCB_1-384.tsv

#!!!!!!!!!!!!!!!!!!!!!!!!! ONLY ADJUST ABOVE THIS LINE !!!!!!!!!!!!!!!!!!!!!!!!#

################################################################################ 
######################################################################### checks

# Absolute path to this script, e.g. /home/user/bin/foo.sh
SCRIPT=$(readlink -f "$0")
# Absolute path this script is in, thus /home/user/bin
SCRIPTPATH=$(dirname "$SCRIPT")


if [ ! -d "fastq" ] 
then
    echo "Input-directory fastq does not exists." 
    exit 9999 
fi
echo "[$params_runid]  Starting MM10-pipeline..."
mkdir $params_runid.filelists
mkdir joblog

################################################################################ 
################################################################## make adapters 
echo "[$params_runid]  - Generating adapters"
Rscript $SCRIPTPATH/mab-id_make_adapters_15062022.Rscript \
$params_sheet $params_sbc_file $params_abbc_file $params_runid.adapters.fa

################################################################################ 
################################################################### merge fastqs
echo "[$params_runid]  - Merging fastqs"
mkdir merged_lanes
cat fastq/*gz > merged_lanes/$params_runid.input.fasta.gz

################################################################################ 
#################################################################### demultiplex
echo "[$params_runid]  - Demultiplexing"
mkdir demux
sbatch -e /tmp/$params_runid.cutadapt.err -o /tmp/$params_runid.cutadapt.out -Q --wait --cpus-per-task 50 --time=08:15:00 --mem=100G --wrap "cutadapt -O 29 -g file:$params_runid.adapters.fa -e 2 -o demux/$params_runid.{name}.fq -j 50 --action=lowercase merged_lanes/$params_runid.input.fasta.gz"

################################################################################ 
#################################################################### parse demux
echo "[$params_runid]  - Cleaning up demux-output"
mkdir demux_parsed
for i in demux/$params_runid.*.fq; do echo "python $SCRIPTPATH/mab-id_parse_fastq_15062022.py $i $params_runid.adapters.fa demux_parsed/$(basename $i | sed 's/fq/fq.gz/')"; done > $params_runid.filelists/$params_runid.parse.filelist

N=$(wc -l $params_runid.filelists/$params_runid.parse.filelist | cut -f 1 -d " ")
sbatch -Q --wait -J $params_runid.parse --time=5:00:00 --mem=20G --array=[1-$N] /hpc/hub_kind/rvanderweide/bin/run_preformatted_cmd_slurm.sh $params_runid.filelists/$params_runid.parse.filelist $N

for i in demux_parsed/*gz; do echo $(basename $i .fq.gz); zcat $i | paste - - - - | wc -l; done | paste - - > parsed_counts.tsv

################################################################################ 
############################################################################ map 
echo "[$params_runid]  - Mapping to mm10"
mkdir aln
for i in demux_parsed/*.gz; do echo "/hpc/hub_kind/rvanderweide/bin/miniconda3/bin/bowtie2 --quiet --very-sensitive -N 1 -p 1 -x /home/hub_kind/rvanderweide/data/MAb-ID/mm10/mm10 -U $i | samtools view -ubh -F 4  - | samtools sort -l 9 -o aln/$(basename $i | sed 's/fq.gz/bam/') - 2> aln/$(basename $i | sed 's/fq.gz/err/')"; done > $params_runid.filelists/$params_runid.aln.filelist

N=$(wc -l $params_runid.filelists/$params_runid.aln.filelist | cut -f 1 -d " ")
sbatch -Q --wait -J $params_runid.aln --time=5:00:00 --mem=20G --array=[1-$N] /hpc/hub_kind/rvanderweide/bin/run_preformatted_cmd_slurm.sh $params_runid.filelists/$params_runid.aln.filelist $N

################################################################################ 
########################################################################## count 
echo "[$params_runid]  - Counting"
mkdir counts
## TTAA
for i in aln/*TTAA*.bam; do echo "python /home/hub_kind/rvanderweide/projects/AbID/AbID_pipeline/pipeline_parts/count_ligations_flatten_UMIs.py -vvv --offset -1 --outfile counts/$(basename $i | sed 's/bam/hdf5/') --min-mapq 10 --umi-length 6 --keep-n 1000 --min-editdistance 1 --pos-file /home/hub_kind/rvanderweide/data/MAb-ID/mm10/mm10_mabid.TTAA.posarray.hdf5 $i 2> counts/$(basename $i | sed 's/bam/log/')"; done > $params_runid.filelists/$params_runid.count.filelist
## GATC
for i in aln/*GATC*.bam; do echo "python /home/hub_kind/rvanderweide/projects/AbID/AbID_pipeline/pipeline_parts/count_ligations_flatten_UMIs.py -vvv --offset 0 --outfile counts/$(basename $i | sed 's/bam/hdf5/') --min-mapq 10 --umi-length 6 --keep-n 1000 --min-editdistance 1 --pos-file /home/hub_kind/rvanderweide/data/MAb-ID/mm10/mm10_mabid.GATC.posarray.hdf5 $i 2> counts/$(basename $i | sed 's/bam/log/')"; done >> $params_runid.filelists/$params_runid.count.filelist

grep -v "\\*" $params_runid.filelists/$params_runid.count.filelist > $params_runid.filelists/$params_runid.count.filelisttmp
mv $params_runid.filelists/$params_runid.count.filelisttmp $params_runid.filelists/$params_runid.count.filelist

N=$(wc -l $params_runid.filelists/$params_runid.count.filelist | cut -f 1 -d " ")
sbatch -Q --wait -J $params_runid.count --time=5:00:00 --mem=20G --array=[1-$N] /hpc/hub_kind/rvanderweide/bin/run_preformatted_cmd_slurm.sh $params_runid.filelists/$params_runid.count.filelist $N

for i in counts/*.log; do echo $i; grep unmapped_reads $i; done | paste - - | sed 's/\[INFO\] //g' | sed 's/\t/ /' | cut -f 1,4-17 -d " " | sed 's/.log//' | sed 's/;//g'> counting_stats.txt

################################################################################ 
############################################################################ bin 
echo "[$params_runid]  - Binning"
## TTAA
for i in counts/*TTAA*.hdf5; do echo "bin_countfile.py --binsize 10000 --posfile /home/hub_kind/rvanderweide/data/MAb-ID/mm10/mm10_mabid.TTAA.posarray.hdf5 --mapfile /home/hub_kind/rvanderweide/data/MAb-ID/mm10/mm10_mabid.TTAA.readlength_70.min_mapq_10.bowtie2.maparray.hdf5 --outfile counts/$(basename $i | sed 's/hdf5/.10k.hdf5/') $i"; done > $params_runid.filelists/$params_runid.bin.filelist

## GATC
for i in counts/*GATC*.hdf5; do echo "bin_countfile.py --binsize 10000 --posfile /home/hub_kind/rvanderweide/data/MAb-ID/mm10/mm10_mabid.GATC.posarray.hdf5 --mapfile /home/hub_kind/rvanderweide/data/MAb-ID/mm10/mm10_mabid.GATC.readlength_70.min_mapq_10.bowtie2.maparray.hdf5 --outfile counts/$(basename $i | sed 's/hdf5/.10k.hdf5/') $i"; done >> $params_runid.filelists/$params_runid.bin.filelist

grep -v "\\*" $params_runid.filelists/$params_runid.bin.filelist > $params_runid.filelists/$params_runid.bin.filelisttmp
mv $params_runid.filelists/$params_runid.bin.filelisttmp $params_runid.filelists/$params_runid.bin.filelist

N=$(wc -l $params_runid.filelists/$params_runid.bin.filelist | cut -f 1 -d " ")
sbatch -Q --wait -J $params_runid.bin --time=1:00:00 --mem=20G --array=[1-$N] /hpc/hub_kind/rvanderweide/bin/run_preformatted_cmd_slurm.sh $params_runid.filelists/$params_runid.bin.filelist $N

################################################################################ 
########################################################################## final
############################################################################ bin 
echo "[$params_runid]  + Done"



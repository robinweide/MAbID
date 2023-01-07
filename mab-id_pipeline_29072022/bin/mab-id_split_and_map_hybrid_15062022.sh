#!/bin/bash
# arg 1 must be a hybrid-mapped file !
# arg 2 must be the number of threads !

# get all reads mapping to hg19
samtools view -F 4 $1 | grep -v m. | grep -v ^@ | \
awk '{print "@"$1"\n"$10"\n+\n"$11}' > \
$(echo $1 | sed 's/.bam/.humanfq/')

# et all reads mapping to mm10
samtools view -F 4 $1 | grep m. | grep -v ^@ | \
awk '{print "@"$1"\n"$10"\n+\n"$11}' > \
$(echo $1 | sed 's/.bam/.mousefq/')

# map hg19_reads to hg19
/hpc/hub_kind/rvanderweide/bin/miniconda3/bin/bowtie2 --quiet --very-sensitive \
  -N 1 -p $2 -x /home/hub_kind/rvanderweide/data/MAb-ID/hg19/hg19 \
  -U $(echo $1 | sed 's/.bam/.humanfq/') | samtools view -ubh -F 4  - | \
  samtools sort -l 9 -o $(echo $1 | sed 's/hybrid.bam/hg19.bam/')

# map mm10_reads to mm10
/hpc/hub_kind/rvanderweide/bin/miniconda3/bin/bowtie2 --quiet --very-sensitive \
  -N 1 -p $2 -x /home/hub_kind/rvanderweide/data/MAb-ID/mm10/mm10 \
  -U $(echo $1 | sed 's/.bam/.mousefq/') | samtools view -ubh -F 4  - | \
  samtools sort -l 9 -o $(echo $1 | sed 's/hybrid.bam/mm10.bam/')
  
samtools view -c $(echo $1 | sed 's/hybrid.bam/hg19.bam/') > $(echo $1 | sed 's/hybrid.bam/hg19.stat/')
samtools view -c $(echo $1 | sed 's/hybrid.bam/mm10.bam/') > $(echo $1 | sed 's/hybrid.bam/mm10.stat/')

[[ -f $(echo $1 | sed 's/hybrid.bam/hg19.stat/') ]] || echo 0 > $(echo $1 | sed 's/hybrid.bam/hg19.stat/')
[[ -f $(echo $1 | sed 's/hybrid.bam/mm10.stat/') ]] || echo 0 > $(echo $1 | sed 's/hybrid.bam/mm10.stat/')

echo $1 | paste - $(echo $1 | sed 's/hybrid.bam/hg19.stat/') $(echo $1 | sed 's/hybrid.bam/mm10.stat/') > $(echo $1 | sed 's/hybrid.bam/hybrid.splitstat/')

rm $(echo $1 | sed 's/hybrid.bam/hg19.stat/')
rm $(echo $1 | sed 's/hybrid.bam/mm10.stat/')
rm $(echo $1 | sed 's/.bam/.humanfq/')
rm $(echo $1 | sed 's/.bam/.mousefq/')
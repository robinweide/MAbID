# Combinatorial single-cell profiling of all major chromatin types with MAbID

Thank you for your interest in MAbID. 

## Table of contents

1. mabid_pipeline
  - manual.html, containing the complete instructions to map raw MAbID data
  - bin, containing all scripts used to map raw MAbID data
2. MAbIDR_1.0.0
  - an R package to analyse MAbID data in R
  - contains help-pages for each tool
  - installation time: 3 minutes on a 2020 Macbook
3. demo
  - test datasets
  - Rmd and html example-code and results 

## Requirements

R (4.1)
Cutadapt (3.0)
Bowtie2 (2.2.9)
samtools (1.10)
scDamAndTools (1.0)

## Installation

To install MAbIDR, run install.packages("MAbIDR_1.0.0.tar.gz", type="source"). 
This should take less then 5 minutes.

To use the mapping pipeline, run one of the file starting with mab-id_pipeline_

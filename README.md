# nanoMD
<!-- TOC -->

- [nanoMD](#nanomd)
- [Overview](#overview)
- [Requirements](#requirements)
- [nanoMD modules](#nanomd-modules)
    - [nanomd gene](#nanomd-gene)
    - [Usage](#usage)
    - [nanomd count](#nanomd-count)
    - [Usage](#usage)
    - [nanomd ployA](#nanomd-ploya)
    - [Usage](#usage)
    - [nanomd matrix](#nanomd-matrix)
    - [Usage](#usage)
    - [nanomd detectMod](#nanomd-detectmod)
    - [Usage](#usage)
    - [nanomd nascentRNA](#nanomd-nascentrna)
    - [Usage](#usage)
- [Scripts](#scripts)
    - [detect5EU.py](#detect5eupy)
    - [Usage](#usage)
- [Docker](#docker)
- [Conda Environment](#conda-environment)
- [Cite nanoMD](#cite-nanomd)

<!-- /TOC -->


# Overview

nanoMD  is designed to work with the output of the ONT direct RNA sequencing platform and can be used to identify new m6A sites, genes, and isoforms, as well as detect and quantify new mRNA. 

It is recommended to use the docker or conda environment to run the pipeline.

# Requirements

1. Python 3.10+
2. Python modules:
    - pandas
    - numpy
    - scipy
    - sklearn
    - matplotlib
    - seaborn
    - pysam
3. minimap2
4. samtools
5. bedtools

# nanoMD modules

nanomd.py is the main module of nanoMD, which includes the following sub-modules:
- gene
- count
- ployA
- matrix
- detectMod
- nascentRNA

## nanomd gene

## Usage

`nanomd gene -i ../input/{}/pass.fq.gz -r ../reference/fasta/genome.fa -o ./{}_gene.sam --parms '--secondary=no --cs -a --sam-hit-only'`

## nanomd count

## Usage

`nanomd count -i ../input/{}/pass.fq.gz -r ../reference/fasta/transcripts.fa -p {}`

## nanomd ployA

## Usage

`nanomd ployA -i ./pass.fq.gz --transcriptome=$ref -o . -p Ctrl-1`

## nanomd matrix

## Usage

`nanomd matrix -i "*_ployA.tsv" -c "WT1_ployA.tsv,WT2_ployA.tsv" -p ployA -s human -t ployA --docker`

`nanomd matrix -i "*_quant" -c "NC1_quant,NC2_quant,NC3_quant,NC4_quant,NC5_quant" -p salmon -s human --docker`

## nanomd detectMod

## Usage

`nanomd detectMod -i ../input/{}/pass.fq.gz -s ../01_map_gene/{}_gene.sam -b ../reference/genes/genes.bed -r ../reference/genes/region_sizes.txt  -p {}`

## nanomd nascentRNA

## Usage

`nanomd nascentRNA -i ../input/{}/pass.fq.gz -s ../01_map_gene/{}_gene.sam -r ../reference/fasta/transcripts.fa -m ~/soft/newRNA.pkl -p {}`

# Scripts

We provide a set of standalone scripts for 5EU detection and quantification.

## detect5EU.py

This script detects 5' untranslated regions (5EU) from the ONT direct RNA sequencing data.

## Usage

`python detect5EU.py -i sample.fastq -o 5EU.bed`

# Docker

If the user has docker installed, the following command can be used to run the pipeline in a docker container:

```
docker run -v /path/to/data:/data -it nanomd/nanomd:latest /bin/bash
```

# Conda Environment

If the user has conda installed, the following command can be used to create a conda environment for nanoMD:

1. Install conda
2. Create a new conda environment: `conda create -n nanomd python=3.12`
3. Activate the environment: `conda activate nanomd`
4. Install the required packages: `conda install -c bioconda minimap2 samtools bedtools`
5. Install the required python packages: `pip install pandas numpy scipy sklearn matplotlib seaborn pysam`
6. Clone the nanoMD repository: `git clone https://github.com/LegendZDY/nanoMD.git`
7. Run the pipeline: `nanomd gene -i ../input/{}/pass.fq.gz -r ../reference/fasta/genome.fa -o ./{}_gene.sam --parms '--secondary=no --cs -a --sam-hit-only'`

# Cite nanoMD

If you use nanoMD in your research, please cite the following paper:
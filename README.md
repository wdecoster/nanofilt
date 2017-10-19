# Nanofilt
Filtering and trimming of Oxford Nanopore sequencing data.

[![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/wouter_decoster.svg?style=social&label=Follow%20%40wouter_decoster)](https://twitter.com/wouter_decoster)
[![conda badge](https://anaconda.org/bioconda/nanofilt/badges/installer/conda.svg)](https://anaconda.org/bioconda/nanofilt)
[![Build Status](https://travis-ci.org/wdecoster/nanofilt.svg?branch=master)](https://travis-ci.org/wdecoster/nanofilt) [![Code Health](https://landscape.io/github/wdecoster/nanofilt/master/landscape.svg?style=flat)](https://landscape.io/github/wdecoster/nanofilt/master)



Filtering on quality and/or read length, and optional trimming after passing filters.  
Reads from stdin, writes to stdout.  

Intended to be used:  
- directly after fastq extraction  
- prior to mapping  
- in a stream between extraction and mapping  

See also [my post about NanoFilt on my blog Gigabase or gigabyte](https://gigabaseorgigabyte.wordpress.com/2017/06/05/trimming-and-filtering-oxford-nanopore-sequencing-reads/).  
Due to [a discrepancy](https://gigabaseorgigabyte.wordpress.com/2017/07/14/calculated-average-quality-vs-albacore-summary/) between calculated read quality and the quality as summarized by albacore this script takes since v1.1.0 optionally also a `--summary` argument. Using this argument with the sequencing_summary.txt file from albacore will do the filtering using the quality scores from the summary. It's also faster.

### INSTALLATION AND UPGRADING:

```bash
pip install nanofilt
pip install nanofilt --upgrade
```
or

[![conda badge](https://anaconda.org/bioconda/nanofilt/badges/installer/conda.svg)](https://anaconda.org/bioconda/nanofilt)
```bash
conda install -c bioconda nanofilt
```
## STATUS
[![Build Status](https://travis-ci.org/wdecoster/nanofilt.svg?branch=master)](https://travis-ci.org/wdecoster/nanofilt) [![Code Health](https://landscape.io/github/wdecoster/nanofilt/master/landscape.svg?style=flat)](https://landscape.io/github/wdecoster/nanofilt/master)


NanoFilt is written for Python 3, and should also work for Python 2.7.

### USAGE:
```
NanoFilt [-h] [-q QUALITY] [-l LENGTH] [--headcrop HEADCROP] [--tailcrop TAILCROP]

optional arguments:  
  -h, --help            show this help message and exit  
  -s --summary SUMMARYFILE optional, the sequencing_summary file from albacore for extracting quality scores
  -q, --quality QUALITY  Filter on a minimum average read quality score  
  -l, --length LENGTH Filter on a minimum read length  
  --headcrop HEADCROP   Trim n nucleotides from start of read  
  --tailcrop TAILCROP   Trim n nucleotides from end of read
  --minGC MINGC         Sequences must have GC content >= to this. Float
                        between 0.0 and 1.0. Ignored if using summary file.
  --maxGC MAXGC         Sequences must have GC content <= to this. Float
                        between 0.0 and 1.0. Ignored if using summary file.
```

Example:
```bash
gunzip -c reads.fastq.gz | NanoFilt -q 10 -l 500 --headcrop 50 | minimap2 genome.fa - | samtools sort -O BAM -@24 -o alignment.bam -
gunzip -c reads.fastq.gz | NanoFilt -q 12 --headcrop 75 | gzip > trimmed-reads.fastq.gz
gunzip -c reads.fastq.gz | NanoFilt -q 10 | gzip > highQuality-reads.fastq.gz
```

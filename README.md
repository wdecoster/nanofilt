# Nanofilt
Filtering and trimming of long read sequencing data.

[![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/wouter_decoster.svg?style=social&label=Follow%20%40wouter_decoster)](https://twitter.com/wouter_decoster)
[![conda badge](https://anaconda.org/bioconda/nanofilt/badges/installer/conda.svg)](https://anaconda.org/bioconda/nanofilt)
[![Build Status](https://travis-ci.org/wdecoster/nanofilt.svg?branch=master)](https://travis-ci.org/wdecoster/nanofilt)



Filtering on quality and/or read length, and optional trimming after passing filters.  
Reads from stdin, writes to stdout.  Optionally reads directly from an uncompressed file specified on the command line.

Intended to be used:  
- directly after fastq extraction  
- prior to mapping  
- in a stream between extraction and mapping  

See also [my post about NanoFilt on my blog Gigabase or gigabyte](https://gigabaseorgigabyte.wordpress.com/2017/06/05/trimming-and-filtering-oxford-nanopore-sequencing-reads/).  
Due to [a discrepancy](https://gigabaseorgigabyte.wordpress.com/2017/07/14/calculated-average-quality-vs-albacore-summary/) between calculated read quality and the quality as summarized by albacore this script takes since v1.1.0 optionally also a `--summary` argument. Using this argument with the sequencing_summary.txt file from albacore will do the filtering using the quality scores from the summary. It's also faster.

### INSTALLATION AND UPGRADING:

`pip install nanofilt`  
`pip install nanofilt --upgrade`

or

`conda install -c bioconda nanofilt`

NanoFilt is written for Python 3.

### USAGE:
```
NanoFilt [-h] [-v] [--logfile LOGFILE] [-l LENGTH]
                [--maxlength MAXLENGTH] [-q QUALITY] [--minGC MINGC]
                [--maxGC MAXGC] [--headcrop HEADCROP] [--tailcrop TAILCROP]
                [-s SUMMARY] [--readtype {1D,2D,1D2}]
                [input]

Perform quality and/or length and/or GC filtering of (long read) fastq data. Reads on stdin.

General options:
  -h, --help            show the help and exit
  -v, --version         Print version and exit.
  --logfile LOGFILE     Specify the path and filename for the log file.
  input                 input, uncompressed fastq file (optional)

Options for filtering reads on.:
  -l, --length LENGTH   Filter on a minimum read length
  --maxlength MAXLENGTH Filter on a maximum read length
  -q, --quality QUALITY Filter on a minimum average read quality score
  --minGC MINGC         Sequences must have GC content >= to this. Float between 0.0 and 1.0. Ignored if
                        using summary file.
  --maxGC MAXGC         Sequences must have GC content <= to this. Float between 0.0 and 1.0. Ignored if
                        using summary file.

Options for trimming reads.:
  --headcrop HEADCROP   Trim n nucleotides from start of read
  --tailcrop TAILCROP   Trim n nucleotides from end of read

Input options.:
  -s, --summary SUMMARY Use albacore or guppy summary file for quality scores
  --readtype            Which read type to extract information about from summary. Options are 1D, 2D or 1D2
 ```

### EXAMPLES
```bash
gunzip -c reads.fastq.gz | NanoFilt -q 10 -l 500 --headcrop 50 | minimap2 genome.fa - | samtools sort -O BAM -@24 -o alignment.bam -
gunzip -c reads.fastq.gz | NanoFilt -q 12 --headcrop 75 | gzip > trimmed-reads.fastq.gz
gunzip -c reads.fastq.gz | NanoFilt -q 10 | gzip > highQuality-reads.fastq.gz
```

I welcome all suggestions, bug reports, feature requests and contributions. Please leave an [issue](https://github.com/wdecoster/nanofilt/issues) or open a pull request. I will usually respond within a day, or rarely within a few days.

## CITATION
If you use this tool, please consider citing our [publication](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty149/4934939).

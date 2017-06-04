# Nanofilt
Filtering and trimming of Oxford Nanopore sequencing data

Filtering on quality and/or read length, and optional trimming after passing filters
Reads from stdin, writes to stdout.  
Intended to be used:
- directly after fastq extraction
- prior to mapping
- in a stream between extraction and mapping

### USAGE:
```
usage: NanoFilt.py [-h] [-q QUALITY] [-l LENGTH] [--headcrop HEADCROP] [--tailcrop TAILCROP]

optional arguments:  
  -h, --help            show this help message and exit  
  -q QUALITY, --quality QUALITY  Filter on a minimum average read quality score  
  -l LENGTH, --length LENGTH Filter on a minimum read length  
  --headcrop HEADCROP   Trim n nucleotides from start of read  
  --tailcrop TAILCROP   Trim n nucleotides from end of read
```

Example:
```bash
zcat reads.fastq.gz | NanoFilt.py -q 10 -l 500 --headcrop 50 | bwa mem -t 48 -x ont2d genome.fa - | samtools sort -O BAM -@24 -o alignment.bam -
```

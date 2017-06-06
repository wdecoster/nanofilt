# Nanofilt
Filtering and trimming of Oxford Nanopore sequencing data.

Filtering on quality and/or read length, and optional trimming after passing filters
Reads from stdin, writes to stdout.  
Intended to be used:
- directly after fastq extraction
- prior to mapping
- in a stream between extraction and mapping

See also [my post about NanoFilt on my blog Gigabase or gigabyte](https://gigabaseorgigabyte.wordpress.com/2017/06/05/trimming-and-filtering-oxford-nanopore-sequencing-reads/).  

### INSTALLATION:

```bash
pip install nanofilt
```

Written for Python 3, might also work for python2.7 (untested).

### USAGE:
```
usage: NanoFilt [-h] [-q QUALITY] [-l LENGTH] [--headcrop HEADCROP] [--tailcrop TAILCROP]

optional arguments:  
  -h, --help            show this help message and exit  
  -q QUALITY, --quality QUALITY  Filter on a minimum average read quality score  
  -l LENGTH, --length LENGTH Filter on a minimum read length  
  --headcrop HEADCROP   Trim n nucleotides from start of read  
  --tailcrop TAILCROP   Trim n nucleotides from end of read
```

Example:
```bash
zcat reads.fastq.gz | NanoFilt -q 10 -l 500 --headcrop 50 | bwa mem -t 48 -x ont2d genome.fa - | samtools sort -O BAM -@24 -o alignment.bam -
zcat reads.fastq.gz | NanoFilt -q 12 --headcrop 75 | gzip > trimmed-reads.fastq.gz
```

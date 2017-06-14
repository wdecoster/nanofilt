#wdecoster
'''
Example usage:
zcat reads.fastq.gz | NanoFilt.py -q 10 -l 500 --headcrop 50 | bwa mem -t 48 -x ont2d genome.fa - | samtools sort -O BAM -@24 -o alignment.bam -
'''

from __future__ import print_function
from Bio import SeqIO
import argparse
import sys
from nanomath import aveQual
from nanofilt.version import __version__


def main():
    args = getArgs()
    filterstream(sys.stdin, args)


def getArgs():
    parser = argparse.ArgumentParser(description="Perform quality and or length filtering of Nanopore fastq data on stdin.")
    parser.add_argument("-v", "--version", help="Print version and exit.", action="version", version='NanoFilt {}'.format(__version__))
    parser.add_argument("-q", "--quality", help="Filter on a minimum average read quality score", default=0, type=int)
    parser.add_argument("-l", "--length", help="Filter on a minimum read length", default=0, type=int)
    parser.add_argument("--headcrop", help="Trim n nucleotides from start of read", default=None, type=int)
    parser.add_argument("--tailcrop", help="Trim n nucleotides from end of read", default=None, type=int)
    return parser.parse_args()


def filterstream(fq, args):
    '''
    If a fastq record passes quality filter (optional) and length filter (optional), print to stdout
    Optionally trim a number of nucleotides from beginning and end.
    '''
    for record in SeqIO.parse(fq, "fastq"):
        if aveQual(record.letter_annotations["phred_quality"]) > args.quality and len(record) > args.length:
            print(record[args.headcrop:args.tailcrop].format("fastq"), end="")

if __name__ == "__main__":
    main()

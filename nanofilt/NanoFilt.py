#wdecoster
'''
Example usage:
gunzip -c reads.fastq.gz | NanoFilt.py -q 10 -l 500 --headcrop 50 | bwa mem -t 48 -x ont2d genome.fa - | samtools sort -O BAM -@24 -o alignment.bam -
'''

from __future__ import print_function
from Bio import SeqIO
import argparse
import sys
from nanomath import aveQual
from nanoget import processSummary
from nanofilt.version import __version__


def main():
    args = getArgs()
    if args.tailcrop:
        args.tailcrop = -args.tailcrop
    if args.summary:
        filterusingSummary(sys.stdin, args)
    else:
        filterstream(sys.stdin, args)


def getArgs():
    parser = argparse.ArgumentParser(description="Perform quality and or length filtering of Nanopore fastq data on stdin.")
    parser.add_argument("-v", "--version",
                        help="Print version and exit.",
                        action="version",
                        version='NanoFilt {}'.format(__version__))
    parser.add_argument("-l", "--length",
                        help="Filter on a minimum read length",
                        default=0,
                        type=int)
    parser.add_argument("--headcrop",
                        help="Trim n nucleotides from start of read",
                        default=None,
                        type=int)
    parser.add_argument("--tailcrop",
                        help="Trim n nucleotides from end of read",
                        default=None,
                        type=int)
    parser.add_argument("-q", "--quality",
                        help="Filter on a minimum average read quality score",
                        default=0,
                        type=int)
    parser.add_argument("-s", "--summary",
                        help="Use summary file for quality scores")
    parser.add_argument("--readtype",
                        help="Which read type to extract information about from summary. Options are 1D or 2D",
                        default="1D",
                        choices=['1D', '2D'])
    return parser.parse_args()


def filterstream(fq, args):
    '''
    If a fastq record passes quality filter (optional) and length filter (optional), print to stdout
    Optionally trim a number of nucleotides from beginning and end.
    '''
    for record in SeqIO.parse(fq, "fastq"):
        if aveQual(record.letter_annotations["phred_quality"]) > args.quality and len(record) > args.length:
            print(record[args.headcrop:args.tailcrop].format("fastq"), end="")


def filterusingSummary(fq, args):
    '''
    Use the summary file from albacore for more accurate quality estimate
    Get the dataframe from nanoget, convert to dictionary
    '''
    data = {entry[0]: entry[1] for entry in processSummary(args.summary, args.readtype)[["readIDs", "quals"]].itertuples(index=False)}
    try:
        for record in SeqIO.parse(fq, "fastq"):
            if data[record.id] > args.quality and len(record) > args.length:
                print(record[args.headcrop:args.tailcrop].format("fastq"), end="")
    except KeyError:
        sys.exit('\nERROR: mismatch between sequencing_summary and fastq file:\n{} was not found in the summary file.\nQuitting.'.format(record.id))


if __name__ == "__main__":
    main()

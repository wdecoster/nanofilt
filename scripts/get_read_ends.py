from argparse import ArgumentParser
from Bio import SeqIO
import gzip


def main():
    args = get_args()
    for record in SeqIO.parse(gzip.open(args.fastq, 'rt'), "fastq"):
        print(record[-args.bases_from_end:].format("fastq"), end="")


def get_args():
    parser = ArgumentParser(description="Filter nanopore data based on time")
    parser.add_argument("fastq", help="input gzip compressed fastq file")
    parser.add_argument("--bases_from_end",
                        help="get a fragment of each read N bp from end",
                        default=100,
                        type=int)
    return parser.parse_args()


if __name__ == '__main__':
    main()

import textwrap as _textwrap
from argparse import ArgumentParser, ArgumentTypeError, HelpFormatter
from nanofilt.version import __version__
import sys
import logging


def get_args():
    epilog = "EXAMPLES:\n" \
        "  gunzip -c reads.fastq.gz | NanoFilt -q 10 -l 500 --headcrop 50 | " \
        "minimap2 genome.fa - | samtools sort -O BAM -@24 -o alignment.bam -\n" \
        "  gunzip -c reads.fastq.gz | NanoFilt -q 12 --headcrop 75 | " \
        "gzip > trimmed-reads.fastq.gz\n" \
        "  gunzip -c reads.fastq.gz | NanoFilt -q 10 | gzip > highQuality-reads.fastq.gz\n"
    parser = ArgumentParser(
        description="Perform quality and/or length and/or GC filtering of (long read) fastq data. \
          Reads on stdin.",
        epilog=epilog,
        formatter_class=custom_formatter,
        add_help=False)
    general = parser.add_argument_group(
        title='General options')
    general.add_argument("-h", "--help",
                         action="help",
                         help="show the help and exit")
    general.add_argument("-v", "--version",
                         help="Print version and exit.",
                         action="version",
                         version='NanoFilt {}'.format(__version__))
    general.add_argument("--logfile",
                         help="Specify the path and filename for the log file.",
                         default="NanoFilt.log")
    general.add_argument("input",
                         help="input, uncompressed fastq file",
                         default=sys.stdin,
                         nargs='?')
    filtering = parser.add_argument_group(
        title='Options for filtering reads on.')
    filtering.add_argument("-l", "--length",
                           help="Filter on a minimum read length",
                           default=1,
                           type=int)
    filtering.add_argument("--maxlength",
                           help="Filter on a maximum read length",
                           default=1e12,
                           type=int)
    filtering.add_argument("-q", "--quality",
                           help="Filter on a minimum average read quality score",
                           default=0,
                           type=int)
    filtering.add_argument("--minGC",
                           help="Sequences must have GC content >= to this.  Float between 0.0 and 1.0. \
                              Ignored if using summary file.",
                           default=0.0,
                           type=valid_GC)
    filtering.add_argument("--maxGC",
                           help="Sequences must have GC content <= to this.  Float between 0.0 and 1.0. \
                              Ignored if using summary file.",
                           default=1.0,
                           type=valid_GC)
    trimming = parser.add_argument_group(
        title='Options for trimming reads.')
    trimming.add_argument("--headcrop",
                          help="Trim n nucleotides from start of read",
                          default=None,
                          type=int)
    trimming.add_argument("--tailcrop",
                          help="Trim n nucleotides from end of read",
                          default=None,
                          type=int)
    inputoptions = parser.add_argument_group(
        title='Input options.')
    inputoptions.add_argument("-s", "--summary",
                              help="Use albacore or guppy summary file for quality scores")
    inputoptions.add_argument("--readtype",
                              help="Which read type to extract information about from summary. \
                              Options are 1D, 2D or 1D2",
                              default="1D",
                              choices=['1D', '2D', "1D2"])
    args = parser.parse_args()
    if args.minGC > args.maxGC:
        sys.exit("NanoFilt: error: argument --minGC should be smaller than --maxGC")
    if args.minGC == 0.0 and args.maxGC == 1.0:
        args.GC_filter = False
    else:
        args.GC_filter = True
    logging.info('NanoFilt {} started with arguments {}'.format(__version__, args))
    return args


class CustomHelpFormatter(HelpFormatter):
    def _fill_text(self, text, width, indent):
        return ''.join(indent + line for line in text.splitlines(keepends=True))

    def _split_lines(self, text, width):
        text = self._whitespace_matcher.sub(' ', text).strip()
        return _textwrap.wrap(text, 80)


def custom_formatter(prog):
    return CustomHelpFormatter(prog)


def valid_GC(x):
    """type function for argparse to check GC values.

    Check if the supplied value for minGC and maxGC is a valid input, being between 0 and 1
    """
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise ArgumentTypeError("{} not in range [0.0, 1.0]".format(x))
    return x


def start_logging(logfile):
    try:
        logging.basicConfig(
            format='%(asctime)s %(message)s',
            filename=logfile,
            level=logging.INFO)
    except PermissionError:
        pass  # indicates that user has no write permission in this directory. No logs then

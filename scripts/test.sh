set -ev

git clone https://github.com/wdecoster/nanotest.git

NanoFilt -h

gunzip -c nanotest/reads.fastq.gz | NanoFilt -q 8 -l 500 --headcrop 50 > /dev/null

gunzip -c nanotest/reads.fastq.gz | NanoFilt -q 7 --maxGC 0.75 --headcrop 75 | gzip > trimmed-reads.fastq.gz

gunzip -c nanotest/reads.fastq.gz | NanoFilt -q 6 -s nanotest/sequencing_summary.txt > /dev/null

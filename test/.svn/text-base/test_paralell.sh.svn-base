
d=$(dirname $(pwd)/${0})
cd $d

perl parallel_trim.pl

mkfifo fastq1 fastq2
zcat -c r1.1.fastq.gz r2.1.fastq.gz r3.1.fastq.gz > fastq1 &
zcat -c r1.2.fastq.gz r2.2.fastq.gz r3.2.fastq.gz > fastq2 &
perl pcat.pl fastq1 fastq2

rm fastq1 fastq2
rm r1.1.fastq.gz r2.1.fastq.gz r3.1.fastq.gz r1.2.fastq.gz r2.2.fastq.gz r3.2.fastq.gz

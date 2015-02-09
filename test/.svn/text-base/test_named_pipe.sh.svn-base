mkfifo fastq1 fastq2
perl ~/perllib/pairTrim/trim.pl --fastq1=/icgc/lsdf/mb/analysis/guz/trim/bsseq_high_qual/bsseq_high_qual.r1.fastq.gz --fastq2=/icgc/lsdf/mb/analysis/guz/trim/bsseq_high_qual/bsseq_high_qual.r2.fastq.gz --output1=fastq1 --output2=fastq2 -s 1 -k 10 &
perl pcat.pl fastq1 fastq2
rm fastq1 fastq2


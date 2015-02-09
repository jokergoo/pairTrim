
d=$(dirname $(pwd)/${0})
cd $d

mkfifo o1 o2
perl ../trim.pl --fastq1=/icgc/lsdf/mb/analysis/guz/trim/bsseq_low_qual/bsseq_low_qual.r1.fastq.gz --fastq2=/icgc/lsdf/mb/analysis/guz/trim/bsseq_low_qual/bsseq_low_qual.r2.fastq.gz --output1=o1 --output2=o2 -k 20  1>/dev/null 2>&1 &

perl pcat.pl o1 o2

rm o1 o2



echo
echo "stucked"

mkfifo o1 o2
perl ../trim.pl --fastq1=/icgc/lsdf/mb/analysis/guz/trim/bsseq_low_qual/bsseq_low_qual.r1.fastq.gz --fastq2=/icgc/lsdf/mb/analysis/guz/trim/bsseq_low_qual/bsseq_low_qual.r2.fastq.gz --output1=o1 --output2=o2 -k 20  1>/dev/null 2>&1 &
# stucked, don't know why
perl cat.pl o1 o2

rm o1 o2


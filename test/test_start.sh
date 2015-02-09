
zcat -c /icgc/lsdf/mb/analysis/guz/trim/bsseq_low_qual/bsseq_low_qual.r2.fastq.gz | head -n 4
mkfifo o2
perl ~/perllib/pairTrim/trim.pl --fastq1=/icgc/lsdf/mb/analysis/guz/trim/bsseq_low_qual/bsseq_low_qual.r1.fastq.gz --fastq2=/icgc/lsdf/mb/analysis/guz/trim/bsseq_low_qual/bsseq_low_qual.r2.fastq.gz --output1=tmp_r1.fastq.gz --output2=o2 -s 1 -k 1 --quality-cutoff=0  1>/dev/null 2>&1 &

cat o2

rm o2
rm tmp_r1.fastq.gz


echo
echo

zcat -c /icgc/lsdf/mb/analysis/guz/trim/bsseq_low_qual/bsseq_low_qual.r2.fastq.gz | head -n 8
mkfifo o2
perl ~/perllib/pairTrim/trim.pl --fastq1=/icgc/lsdf/mb/analysis/guz/trim/bsseq_low_qual/bsseq_low_qual.r1.fastq.gz --fastq2=/icgc/lsdf/mb/analysis/guz/trim/bsseq_low_qual/bsseq_low_qual.r2.fastq.gz --output1=tmp_r1.fastq.gz --output2=o2 -s 1 -k 2 --quality-cutoff=0  1>/dev/null 2>&1 &

cat o2

rm o2
rm tmp_r1.fastq.gz


echo
echo

zcat -c /icgc/lsdf/mb/analysis/guz/trim/bsseq_low_qual/bsseq_low_qual.r2.fastq.gz | awk 'NR==4,NR==8'
mkfifo o2
perl ~/perllib/pairTrim/trim.pl --fastq1=/icgc/lsdf/mb/analysis/guz/trim/bsseq_low_qual/bsseq_low_qual.r1.fastq.gz --fastq2=/icgc/lsdf/mb/analysis/guz/trim/bsseq_low_qual/bsseq_low_qual.r2.fastq.gz --output1=tmp_r1.fastq.gz --output2=o2 -s 2 -k 1 --quality-cutoff=0  1>/dev/null 2>&1 &

cat o2

rm o2
rm tmp_r1.fastq.gz


echo
echo

zcat -c /icgc/lsdf/mb/analysis/guz/trim/bsseq_low_qual/bsseq_low_qual.r2.fastq.gz | awk 'NR==4,NR==12'
mkfifo o2
perl ~/perllib/pairTrim/trim.pl --fastq1=/icgc/lsdf/mb/analysis/guz/trim/bsseq_low_qual/bsseq_low_qual.r1.fastq.gz --fastq2=/icgc/lsdf/mb/analysis/guz/trim/bsseq_low_qual/bsseq_low_qual.r2.fastq.gz --output1=tmp_r1.fastq.gz --output2=o2 -s 2 -k 2 --quality-cutoff=0  1>/dev/null 2>&1 &

cat o2

rm o2
rm tmp_r1.fastq.gz

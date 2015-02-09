
d=$(dirname $(pwd)/${0})
cd $d



time1=`date +%s`
perl /home/guz/ngspipeline/BSTools/bismark_bin/trim_galore --paired --gzip -o /icgc/lsdf/mb/analysis/guz/trim/bsseq_low_qual /icgc/lsdf/mb/analysis/guz/trim/bsseq_low_qual/bsseq_low_qual.r1.fastq.gz /icgc/lsdf/mb/analysis/guz/trim/bsseq_low_qual/bsseq_low_qual.r2.fastq.gz 1>/dev/null 2>&1
time2=`date +%s`
echo Trim Galore used $(( $time2 - $time1 )) s.

time1=`date +%s`
perl ../trim.pl --fastq1=/icgc/lsdf/mb/analysis/guz/trim/bsseq_low_qual/bsseq_low_qual.r1.fastq.gz --fastq2=/icgc/lsdf/mb/analysis/guz/trim/bsseq_low_qual/bsseq_low_qual.r2.fastq.gz --output1=/icgc/lsdf/mb/analysis/guz/trim/bsseq_low_qual/bsseq_low_qual.r1.trimmed.fastq.gz --output2=/icgc/lsdf/mb/analysis/guz/trim/bsseq_low_qual/bsseq_low_qual.r2.trimmed.fastq.gz --overlap=1 1>/dev/null 2>&1
time2=`date +%s`
echo pairTrim used $(( $time2 - $time1 )) s.

check_1=$(zcat -c /icgc/lsdf/mb/analysis/guz/trim/bsseq_low_qual/bsseq_low_qual.r1_val_1.fq.gz | md5sum | awk '{print $1}')
check_2=$(zcat -c /icgc/lsdf/mb/analysis/guz/trim/bsseq_low_qual/bsseq_low_qual.r1.trimmed.fastq.gz | md5sum | awk '{print $1}')
check_3=$(zcat -c /icgc/lsdf/mb/analysis/guz/trim/bsseq_low_qual/bsseq_low_qual.r2_val_2.fq.gz | md5sum | awk '{print $1}')
check_4=$(zcat -c /icgc/lsdf/mb/analysis/guz/trim/bsseq_low_qual/bsseq_low_qual.r2.trimmed.fastq.gz | md5sum | awk '{print $1}')

if [[ $check_1 == $check_2 ]] && [[ $check_3 == $check_4 ]];
then
	echo "output same for TrimGalore and pairTrim.";
else
	echo "ERROR: output differetn form TrimGalore and pairTrim.";
fi


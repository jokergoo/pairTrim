

my $k = 30;

# child process to trim R1
my $pid1 = fork();
if($pid1) {
	#print "pid is $pid1, parent $$\n";
} elsif ($pid1 == 0) {

	run(1, $k);
	exit 0;
	
} else {
	die "cannot fork: $!\n";
}

my $pid2 = fork();
if($pid2) {
	#print "pid is $pid1, parent $$\n";
} elsif ($pid2 == 0) {

	run(2, $k);
	exit 0;
	
} else {
	die "cannot fork: $!\n";
}

my $pid3 = fork();
if($pid3) {
	#print "pid is $pid1, parent $$\n";
} elsif ($pid3 == 0) {

	run(3, $k);
	exit 0;
	
} else {
	die "cannot fork: $!\n";
}


waitpid($pid1, 0);
waitpid($pid2, 0);
waitpid($pid3, 0);

sub run {
	my $i = shift;
	my $k = shift;
	my $start = $k*($i-1)+1;
	`perl ~/perllib/pairTrim/trim.pl --fastq1=/icgc/lsdf/mb/analysis/guz/trim/bsseq_low_qual/bsseq_low_qual.r1.fastq.gz --fastq2=/icgc/lsdf/mb/analysis/guz/trim/bsseq_low_qual/bsseq_low_qual.r2.fastq.gz --output1=r$i.1.fastq.gz --output2=r$i.2.fastq.gz -s $start -k $k 1>/dev/null 2>&1`;
}
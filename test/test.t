use strict;
use Test::More;
use File::Basename;

chdir dirname($0)."/../";

my $f1;
my $f2;
my $o1;
my $o2;

################################################################################
# test trimming ending A
################################################################################
# left
($f1, $f2, $o1, $o2) = generate_paired_end(
	seq  => "AAAAAAAAAAAACATCGTACGTACGTATGCGTACGTAGTCGTATGCGTAGTTCACGTCGTAGTCGTACAGTCGT",
	qual => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII");
`perl trim.pl --fastq1=$f1 --fastq2=$f2 --output1=$o1 --output2=$o2 --adapter=CCCCCCC --quality-cutoff=0 --trim-polyA 1>/dev/null 2>&1`;
my $seq = get_seq($o1);
cmp_ok("TCGTACGTACGTATGCGTACGTAGTCGTATGCGTAGTTCACGTCGTAGTCGTACAGTCGT", "eq", $seq, "testing trimming ending As on the left");
unlink($f1); unlink($f2); unlink($o1); unlink($o2);

# left, less than 5
($f1, $f2, $o1, $o2) = generate_paired_end(
	seq  => "AAACTTCGTACGTACGTATGCGTACGTAGTCGTATGCGTAGTTCACGTCGTAGTCGTACAGTCGT",
	qual => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII");
`perl trim.pl --fastq1=$f1 --fastq2=$f2 --output1=$o1 --output2=$o2 --adapter=CCCCCCC --quality-cutoff=0 --trim-polyA 1>/dev/null 2>&1`;
my $seq = get_seq($o1);
cmp_ok("AAACTTCGTACGTACGTATGCGTACGTAGTCGTATGCGTAGTTCACGTCGTAGTCGTACAGTCGT", "eq", $seq, "testing trimming ending As on the left, less than 5");
unlink($f1); unlink($f2); unlink($o1); unlink($o2);

# left, with mismatch, passed
($f1, $f2, $o1, $o2) = generate_paired_end(
	seq  => "AAAAAAATAAAACATCGTACGTACGTATGCGTACGTAGTCGTATGCGTAGTTCACGTCGTAGTCGTACAGTCGT",
	qual => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII");
`perl trim.pl --fastq1=$f1 --fastq2=$f2 --output1=$o1 --output2=$o2 --adapter=CCCCCCC --quality-cutoff=0 --trim-polyA 1>/dev/null 2>&1`;
my $seq = get_seq($o1);
cmp_ok("CATCGTACGTACGTATGCGTACGTAGTCGTATGCGTAGTTCACGTCGTAGTCGTACAGTCGT", "eq", $seq, "testing trimming ending As on the left, with mismatch");
unlink($f1); unlink($f2); unlink($o1); unlink($o2);

# left, with mismatch, not passed
($f1, $f2, $o1, $o2) = generate_paired_end(
	seq  => "AAAAAAATTAAACATCGTACGTACGTATGCGTACGTAGTCGTATGCGTAGTTCACGTCGTAGTCGTACAGTCGT",
	qual => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII");
`perl trim.pl --fastq1=$f1 --fastq2=$f2 --output1=$o1 --output2=$o2 --adapter=CCCCCCC --quality-cutoff=0 --trim-polyA 1>/dev/null 2>&1`;
my $seq = get_seq($o1);
cmp_ok("TTAAACATCGTACGTACGTATGCGTACGTAGTCGTATGCGTAGTTCACGTCGTAGTCGTACAGTCGT", "eq", $seq, "testing trimming ending As on the left, with more mismatch");
unlink($f1); unlink($f2); unlink($o1); unlink($o2);

# right
($f1, $f2, $o1, $o2) = generate_paired_end(
	seq  => "CATCGTACGTACGTATGCGTACGTAGTCGTATGCGTAGTTCACGTCGTAGTCGTACAGTCGTAAAAAAAAAAAA",
	qual => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII");
`perl trim.pl --fastq1=$f1 --fastq2=$f2 --output1=$o1 --output2=$o2 --adapter=CCCCCCC --quality-cutoff=0 --trim-polyA 1>/dev/null 2>&1`;
my $seq = get_seq($o1);
cmp_ok("CATCGTACGTACGTATGCGTACGTAGTCGTATGCGTAGTTCACGTCGTAGTCGTACAGTCGT", "eq", $seq, "testing trimming ending As on the right");
unlink($f1); unlink($f2); unlink($o1); unlink($o2);

# right, less than 5
($f1, $f2, $o1, $o2) = generate_paired_end(
	seq  => "CATCGTACGTACGTATGCGTACGTAGTCGTATGCGTAGTTCACGTCGTAGTCGTACAGTCGTAAA",
	qual => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII");
`perl trim.pl --fastq1=$f1 --fastq2=$f2 --output1=$o1 --output2=$o2 --adapter=CCCCCCC --quality-cutoff=0 --trim-polyA 1>/dev/null 2>&1`;
my $seq = get_seq($o1);
cmp_ok("CATCGTACGTACGTATGCGTACGTAGTCGTATGCGTAGTTCACGTCGTAGTCGTACAGTCGTAAA", "eq", $seq, "testing trimming ending As on the right, less than 5");
unlink($f1); unlink($f2); unlink($o1); unlink($o2);

# right, with mismatch, passed
($f1, $f2, $o1, $o2) = generate_paired_end(
	seq  => "CATCGTACGTACGTATGCGTACGTAGTCGTATGCGTAGTTCACGTCGTAGTCGTACAGTCGTAAAATAAAAAAA",
	qual => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII");
`perl trim.pl --fastq1=$f1 --fastq2=$f2 --output1=$o1 --output2=$o2 --adapter=CCCCCCC --quality-cutoff=0 --trim-polyA 1>/dev/null 2>&1`;
my $seq = get_seq($o1);
cmp_ok("CATCGTACGTACGTATGCGTACGTAGTCGTATGCGTAGTTCACGTCGTAGTCGTACAGTCGT", "eq", $seq, "testing trimming ending As on the right, with mismatch");
unlink($f1); unlink($f2); unlink($o1); unlink($o2);

# right, with mismatch, not passed
($f1, $f2, $o1, $o2) = generate_paired_end(
	seq  => "CATCGTACGTACGTATGCGTACGTAGTCGTATGCGTAGTTCACGTCGTAGTCGTACAGTCGTAAATTAAAAAAA",
	qual => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII");
`perl trim.pl --fastq1=$f1 --fastq2=$f2 --output1=$o1 --output2=$o2 --adapter=CCCCCCC --quality-cutoff=0 --trim-polyA 1>/dev/null 2>&1`;
my $seq = get_seq($o1);
cmp_ok("CATCGTACGTACGTATGCGTACGTAGTCGTATGCGTAGTTCACGTCGTAGTCGTACAGTCGTAAATT", "eq", $seq, "testing trimming ending As on the right, with more mismatch");
unlink($f1); unlink($f2); unlink($o1); unlink($o2);


done_testing();

sub generate_paired_end {
	my %param = ( seq => undef,
	              qual => undef,
                     @_);
	
	my $r = $$.int(rand(100000));
	my $f1 = "$r.r1.fastq";
	my $f2 = "$r.r2.fastq";
	my $o1 = "$r.r1_trimmed.fastq";
	my $o2 = "$r.r2_trimmed.fastq";
	
	open OUT, ">$f1" or die "cannot create $f1\n";;
	print OUT "\@read 1\n";
	print OUT "$param{seq}\n";
	print OUT "+\n";
	print OUT "$param{qual}\n";
	close OUT;
	
	open OUT, ">$f2" or die "cannot create $f2\n";;
	print OUT "\@read 2\n";
	print OUT "$param{seq}\n";
	print OUT "+\n";
	print OUT "$param{qual}\n";
	close OUT;
				  
	return ($f1, $f2, $o1, $o2);
}

sub get_seq {
	my $f = shift;
	
	open F, $f or die "cannot open $f\n";
	<F>;
	my $seq = <F>;
	close F;
	
	chomp $seq;
	return $seq;
}
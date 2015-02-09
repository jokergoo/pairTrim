use strict;

# get the script's directory. It needs to find where CO::... modules are.
BEGIN {
	use File::Basename;
	unshift(@INC, dirname($0)."");
}

use CO::FastQ;
use Data::Dumper;
use File::Basename;
use Getopt::Long;
use Fcntl;

# print help message
if(scalar(@ARGV) == 0 or grep {/^(-h|--help)$/i} @ARGV) {
	print_msg();
	exit;
}


# read the parameters, these parameters are global variables which will
# be assigned in parse_arg()
our $FASTQ1;
our $FASTQ2;
our $OUTPUT1;
our $OUTPUT2;
our $TMP;
our $START;
our $RECORDS;
our $MINIMUM_LENGTH;
our $QUALITY_CUTOFF;
our $QUALITY_BASE;
our $TRIM_POLYA;
our $LEFT;
our $RIGHT;

# construct arguments for cutadapt
# if user specified multiple adapters, the script will trim all adapters one by
# one through pipes. So for each adapter, there is a cutadapt command.
our @CUTADAPT_ARGS;
parse_arg();


# print parameters
print STDERR <<CONF;

fastq1         = $FASTQ1
fastq2         = $FASTQ2
output1        = $OUTPUT1
output2        = $OUTPUT2
tmp            = $TMP
start          = $START
records        = $RECORDS
minimum-length = $MINIMUM_LENGTH
quality-cutoff = $QUALITY_CUTOFF
quality-base   = $QUALITY_BASE
trim-polyA     = $TRIM_POLYA

CONF

# create named pipe
# If there are multiple lanes, each lane will evoke this script and will get
# differetn pid.
my $r = time().int(rand(9999999));
my $fifo1 = "$TMP/$r.$$.output1";
my $fifo2 = "$TMP/$r.$$.output2";
system("mkfifo $fifo1");
system("mkfifo $fifo2");


################################################################################
# create children processes

# child process to trim R1
my $pid1 = fork();
if($pid1) {
	#print "pid is $pid1, parent $$\n";
} elsif ($pid1 == 0) {

	# for the child process, call cutadapt and exit
	cutadapt($FASTQ1, $START, $RECORDS, $fifo1);
	exit 0;
	
} else {
	die "cannot fork: $!\n";
}

# child process to trim R2
my $pid2 = fork();
if($pid2) {
	#print "pid is $pid2, parent $$\n";
} elsif ($pid2 == 0) {

	# for the child process, call cutadapt and exit
	cutadapt($FASTQ2, $START, $RECORDS, $fifo2);
	exit 0;
	
} else {
	die "cannot fork: $!\n";
}

# child process to filter trimmed reads and keep paired reads
my $pid3 = fork();
if($pid3) {
	#print "pid is $pid3, parent $$\n";
} elsif ($pid3 == 0) {

	# receive trimmed fastq from cutadapt and post process
	post($fifo1, $fifo2, $OUTPUT1, $OUTPUT2);
	exit 0;
	
} else {
	die "cannot fork: $!\n";
}

# recycle child processes
waitpid($pid1, 0);
waitpid($pid2, 0);
waitpid($pid3, 0);

# delete named pipes
system("rm $fifo1");
system("rm $fifo2");

# run cutadapt
# note only value for --minimum_length is not sent to cutadapt, So cutadapt will
# output all reads even if the read length is 0. It ensures numbers and orders
# of reads in paired-end fastq files are the same.
#
# Too short reads will be removed and synchronized in `post` step
sub cutadapt {
	my $file = shift;  # fastq or gzipped fastq
	my $start = shift; # which read to start
	my $n_records = shift;
	my $output = shift;  # can be a named pipe
	
	# through a series of pipes, multiple adapters will be trimmed
	# finally write to the named pipe
	my $cutadapt_command = join " | ", map {"cutadapt --format=fastq $_ -"} @CUTADAPT_ARGS;
	$cutadapt_command = "$cutadapt_command --output=$output";
	
	print STDERR "$cutadapt_command\n";
	
	# transfer fastq to cutadapt through pipe
 	open my $out, "| $cutadapt_command" or die "cannot create pipe for cutadapt.\n";
	
	my $fastq = CO::FastQ->new(file => $file);

	my $i_processed = 0;
	
	# each time read into one read
	# only send read from `start` to `start + n_records` to cutadapt
	while(my $read = $fastq->next) {
	
		if($fastq->i < $start) {
			next;
		}
		
		$read->trim_n(left => $LEFT,
		              right => $RIGHT);
		
		print $out $read->record;
		
		$i_processed ++;
		if($i_processed % 1e5 == 0) {
			print STDERR "[".basename($file)."] $i_processed reads sent to cutadapt.\n";
		}
		
		if($i_processed >= $n_records) {
			last;
		}
	}
}

# read output of two reads and post processing
sub post {
	my $f1 = shift;    # named pipe 1
	my $f2 = shift;    # named pipe 2
	my $gz1 = shift;   # output 1, fastq or gzipped fastq
	my $gz2 = shift;   # output 2, fastq or gzipped fastq
	
	my $fastq1 = CO::FastQ->new(file => $f1);  # trimmed fastq stream 1
	my $fastq2 = CO::FastQ->new(file => $f2);  # trimmed fastq stream 2
	
	my ($fh1, $fh2);
	if(-p $gz1) {
		sysopen($fh1, "$gz1", O_WRONLY);
	} elsif($gz1 =~/\.gz$/) {
		open $fh1, " | gzip -c > $gz1" or die "cannot create pipe to gzip\n";
	} else {
		open $fh1, ">$gz1" or die "cannot create $gz1\n";
	}
	
	if(-p $gz2) {
		sysopen($fh2, "$gz2", O_WRONLY);
	} elsif($gz2 =~/\.gz$/) {
		open $fh2, " | gzip -c > $gz2" or die "cannot create pipe to gzip\n";
	} else {
		open $fh2, ">$gz2" or die "cannot create $gz2\n";
	}
	
	my $i_processed = 0;
	my $i_accepted = 0;
	my $i_drop_cutadapt = 0;
	my $i_drop_poly_A = 0;
	my $read1;
	my $read2;
	my $l1;
	my $l2;
	
	# $read1 and $read2 are with the same name
	while($read1 = $fastq1->next and $read2 = $fastq2->next) {
		
		$i_processed ++;

		if($i_processed % 1e5 == 0) {
		
			# in fact, we do not need to check read name in every cycle.
			# If reads in #i pair have different name, it's almost true that
			# reads in #i + 1e5 will also have different name.
			# Checking read name in every cycle will make efficiency low.
			if($read1->name ne $read2->name) {
				die "Read pair have different names (${[$read1->name]}, ${[$read2->name]}), before #$i_processed.\n";
			}
			print STDERR "$i_processed pairs processed.\n";
		}
		
		if($read1->length < $MINIMUM_LENGTH or $read2->length < $MINIMUM_LENGTH) {
			$i_drop_cutadapt ++;
			next;
		}
		
		# then do a list of post trimming
			
		if($TRIM_POLYA) {

			$read1->trim_ploy_A;
			if($read1->length < $MINIMUM_LENGTH) {
				$i_drop_poly_A ++;
				next;
			}

			$read2->trim_ploy_A;
			if($read2->length < $MINIMUM_LENGTH) {
				$i_drop_poly_A ++;
				next;
			}
		}
		
		$i_accepted ++;
		print $fh1 $read1->record;
		print $fh2 $read2->record;
			
		$l1->{$read1->length} ++;
		$l2->{$read2->length} ++;

	}
	
	close $fh1;
	close $fh2;
	
	my $rate = sprintf '%.2f', $i_accepted/$i_processed*100;
	print STDERR "All reads pairs: $i_processed\n";
	print STDERR "Accepted reads pairs: $i_accepted ($rate\%)\n";
	print STDERR "Short reads pairs (< $MINIMUM_LENGTH) dropped after trimming by cutadapt: $i_drop_cutadapt\n";
	print STDERR "Short reads pairs (< $MINIMUM_LENGTH) dropped after removing poly-A tails: $i_drop_poly_A\n" if ($TRIM_POLYA);
	
	print STDERR "\n";
	print STDERR "Read length in R1:\n";
	foreach (sort {$a <=> $b} keys %$l1) {
		print STDERR "$_\t$l1->{$_}\n";
	}
	
	print STDERR "\n";
	print STDERR "Read length in R2:\n";
	foreach (sort {$a <=> $b} keys %$l2) {
		print STDERR "$_\t$l2->{$_}\n";
	}
	
}

sub parse_arg {

	# parameters started with '#' will be sent to trimPair
	# parameters not started with '#' will be sent to cutadapt
	my $param = { 'adapter' => [],
	        #      'anywhere' => [],      # currently, do not support --anywhere
	        #      'front' => [],         # currently, do not support --front
	              'error-rate' => 0.1,
	              'times' => 1,
	              'overlap' => 3,
	              'match-read-wildcards!' => undef,
	              'no-match-adapter-wildcards!' => undef,
	              'discard!' => undef,
	              'minimum-length' => 0,
	              'maximum-length' => undef,
	              'quality-cutoff' => 20,
	              'quality-base' => 33,
	              'prefix' => undef,
	              'suffix' => undef, 
	              '#fastq1' => undef,
	              '#fastq2' => undef,
	              '#output1' => undef,
	              '#output2' => undef,
	              '#tmp' => './',
	              '#start' => 1,
	              '#records' => 1e20,
	              '#minimum_length' => 20,
	              '#trim-polyA' => 0,
	              '#left' => 0,
	              '#right' => 0};
	              
	
	GetOptions("adapter|a=s@"                 => $param->{'adapter'},
	         #  "anywhere|b=s@"                => $param->{'anywhere'},
	         #  "front|g=s@"                   => $param->{'front'},
	           "error-rate|e=f"               => \$param->{'error-rate'},
	           "times|n=i"                    => \$param->{'times'},
	           "overlap|O=i"                  => \$param->{'overlap'},
	           "match-read-wildcards"         => \$param->{'match-read-wildcards!'},
	           "no-match-adapter-wildcards|N" => \$param->{'no-match-adapter-wildcards!'},
	           "discard-trimmed|discard"      => \$param->{'discard!'},
	           "minimum-length|m=i"           => \$param->{'#minimum_length'},
	           "maximum-length|M=i"           => \$param->{'maximum-length'},
	           "quality-cutoff|q=f"           => \$param->{'quality-cutoff'},
	           "quality-base=i"               => \$param->{'quality-base'},      
	           "prefix|x=s"                   => \$param->{'prefix'},
	           "suffix|y=s"                   => \$param->{'suffix'},
	           "fastq1|f1|1=s"                => \$param->{'#fastq1'},
	           "fastq2|f2|2=s"                => \$param->{'#fastq2'},
	           "output1|o1|3=s"               => \$param->{'#output1'},
	           "output2|o2|4=s"               => \$param->{'#output2'},
	           "tmp-dir|tmp|t=s"              => \$param->{'#tmp'},
	           "start|s=i"                    => \$param->{'#start'},
	           "records|k=i"                  => \$param->{'#records'},
	           "trim-polyA|A"                 => \$param->{'#trim_polyA'},
	           "left=i"                         => \$param->{'#left'},
	           "right=i"                        => \$param->{'#right'},
	           )
	or die;
	
	# maybe useful in future
	$param->{'#quality_cutoff'} = $param->{'quality-cutoff'};
	$param->{'#quality_base'} = $param->{'quality-base'};
	
	# default adapter
	if(scalar(@{$param->{adapter}}) ==0) {
		$param->{adapter} = ["AGATCGGAAGAGC"],
	}
	
	check_mandatory("fastq1", $param->{'#fastq1'});
	check_mandatory("fastq2", $param->{'#fastq2'});
	check_mandatory("output1", $param->{'#output1'});
	check_mandatory("output2", $param->{'#output2'});
	
	my @str;
	for(my $i = 0; $i < scalar(@{$param->{adapter}}); $i ++) {
		$str[$i] = "";
	}
	
	# assign parameters
	foreach my $opt (keys %$param) {
		
		# global variables for pairTrim
		if($opt =~/^#/) {
			my $v = $param->{$opt};
			$opt =~s/^#//;
			{
				no strict;
				${uc($opt)} = $v;
			}
			
		} elsif(ref($param->{$opt}) eq "ARRAY") {
		
			for(my $i = 0; $i < scalar(@str); $i ++) {
				$str[$i] .= "--$opt $param->{$opt}->[$i] ";
			}
			
		} elsif($opt =~/!$/) {  # end with ! means this option has no value (just a flag)
		
			if(defined($param->{$opt})) {
				$opt =~s/!$//;
				for(my $i = 0; $i < scalar(@str); $i ++) {
					$str[$i] .= "--$opt ";
				}
			}
			
		} elsif(defined($param->{$opt})) {
		
			for(my $i = 0; $i < scalar(@str); $i ++) {
				$str[$i] .= "--$opt $param->{$opt} ";
			}
			
		}
		
	}
	
	@CUTADAPT_ARGS = @str;
	
	return @str;
}

sub print_msg {
	my $fh = shift || *STDERR;
	
	my $script_name = basename($0);
	
	print $fh <<MSG;

Trimming paired end FastQ files.

    perl $script_name --fastq1=f1 --fastq2=f2 --output1=o1 --output2=o2
 
trimPair is used to trim paired-end FastQ files. It wraps `cutadapt` 
to trim adapters and uses self-defined codes to remove poly-A/T/N on both ends 
of the reads. Here, post trimming allow mismatches. It also support additional
features such as trim first N or last N based before going to `cutadpat`.

Currently it does not support trimming single-end FastQ files.

trimPair initializes three child processes in which two processes trim reads in
each pair, and the last processes do synchronization from two pairs and do some
post processing. Data is transferred between processes by named pipes which 
can get avoid of I/O error.

The output can be send to downstreaming programs by named pipe (just don't named
the pipe as *.gz).

If --overlap set to 1, it generates exactly the same output as Trim Galore.

NOTE:

  It assumes your paired-end FastQ files are valid (Same number and order of 
  reads in two pairs, four line per read).

Contact:

  Zuguang Gu <z.gu\@dkfz.de>

Parameters:

  --fastq1, --f1, -1   Path of FastQ file, pair 1. Either gzipped or not.
  --fastq2, --f2, -2   Path of FastQ file, pair 2. Either gzipped or not.
  --output1, --o1, -3  Path of output trimmed FastQ file, pair 1. Either gzipped
                       or not. It can also be a named pipe
  --output2, --o2, -4  Path of output trimmed FastQ file, pair 2. Either gzipped
                       or not. It can also be a named pipe
  --tmp, -t            Path of temporary directory.
  --start, -s          Which read in the FastQ files starts for trimming. 
                       Default value is 1 which means start form the beginning.
  --records, -k        How many reads should be trimmed. Default value is 1e20
                       which means to the end.
  --trim-polyA, -A     Whether to trim poly-A. Default FALSE. If it is set to
                       TRUE, The script will trim poly A or poly T from two end
                       of the read with mismatch.
  --process            It is a disabled parameter. It just tells you it will use
                       3 CPUs for a job ( in fact a little bit less than 3).
  --minimum-length, -m Trimmed reads with length less than this value will be
                       dropped. Default is 20.
  --left               How many bases should be firstly trimmed from the left.
  --right              How many bases should be firstly trimmed from the right.
  --help, -h           Help message.
  
  Together with --start and --records, it will be easy to make trimming parallel.
  For example, we spliced all reads into 4 pieces in which each piece contains 
  N reads, then run following commands:
  
    perl $script_name -f1 r1.gz -f2 r2.gz -o1 t1_1.gz -o2 t2_1.gz -s    1 -k N &
    perl $script_name -f1 r1.gz -f2 r2.gz -o1 t1_2.gz -o2 t2_2.gz -s  N+1 -k N &
    perl $script_name -f1 r1.gz -f2 r2.gz -o1 t1_3.gz -o2 t2_3.gz -s 2N+1 -k N &
    perl $script_name -f1 r1.gz -f2 r2.gz -o1 t1_4.gz -o2 t2_4.gz -s 3N+1 -k N &
  
  You can send trimmed FastQ to downstreaming programs:
  
    mkfifo fastq1 fastq2
    perl $script_name -f1 r1.gz -f2 r2.gz -o1 fastq1 -o2 fastq2 &
    cmd fastq1 fastq1
    rm fastq1 fastq2

Parameters that are sent to cutadapt:

  -a, --adapter  Different to cutadapt, if multiple adapters are specified, 
                 all adapters will be trimmed. Examples are:
                 --adapter=A1 --adapter=A2
  -e, --error-rate
  -q, --quality-cutoff
  -i, --quality-base
  -M, --maximum-length
  -n, --times
  -O, --overlap
      --match-read-wildcards
  -N, --no-match-adapter-wildcards
      --discard-trimmed, --discard
  -x, --prefix
  -y, --suffix

MSG
}

sub check_mandatory {
	my $field = shift;
	my $path = shift;
	
	if(!defined($path)) {
		print STDERR "--$field should be specified.\n";
		exit;
	}

}

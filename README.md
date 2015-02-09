Trim paired-end FastQ files
===========================

## Overview:

`trimPair` is used to trim paired-end FastQ files. It wraps `cutadapt` to trim 
adapters and uses self-defined codes to remove poly-A/T/N on both ends of the 
reads.

`trimPair` initializes three child processes in which two processes trim reads in
both pairs, one processes do synchronization from two pairs and post processing.

The output can be send to downstreaming programs by named pipe (just don't named
the pipe as `*.gz`).

If `--overlap` is set to 1, it generates exactly the same output as `Trim Galore`.

## NOTE:

  It assumes your paired-end FastQ files are valid (Same number and order of 
  reads in two pairs, four line per read).

## Parameters:

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
  
Together with `--start` and `--records`, it will be easy to make trimming parallel.
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

Parameters that are sent to `cutadapt`:

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


## mismatch model

This is a simple mismatch model, for each base, this is a value (1|0) returned
by a mismatch subroutine. For example, if we want to trim poly-N, the mismatch
subroutine will be defined as 

    sub {$_[0] eq "N"}
    
It returns 1 if current base is 'N'.

Start from each end, we frist look at the first 5 bases. If the first 5 bases
contain no mismatch or 1 mismatch, we will keep on looking at following bases.
When going through bases from end to center, we calculate mismatch rate for each
base at store the last matched position. If the mismatch rate is larger than a 
cutoff. We continue reading next 10 bases. If mismatch rate never drops, we will
remove bases to the 'last matched position'. If mismatch rate drops lower than
cutoff, we continue reading following bases and repeat above steps.

In case the first five bases which have more than 1 mismatches, for poly-A/T, 
we do not trim anything since the number of ending A/T is too small that we can
not take them as poly-A tail.

## Files:

    trim.pl:
      main script
    
    test/:
      testing script
    
    CO/FastQ.pm
      CO::FastQ class, for a whole FastQ file
        element:
          fh: file handle
          i : index of current read
        methods:
          next: get next read, returns CO::FastQ::Read object
          i: index of current read
    
    CO/FastQ/Read.pm
      CO::FastQ::Read class, for a single read
        element:
          line1: reference to line1 of a read record
          line2: reference to line2 of a read record
          line3: reference to line3 of a read record
          line4: reference to line4 of a read record
          begin: position of the beginning base after trimming
          end  : position of the ending base after trimming
        methods:
          header: first line of the read record
          name: read name, extract from header
          seq: read sequence
          qual_str: read quality string
          qual: read quanlity, transformed into numbers, array reference
          length: length of read after trimming
    
    CO::FastQ::Trim.pm
      Super class of CO::FastQ::Read, provide trimming methods for CO::FastQ::Read class
        methods:
          set: set some parameters for trimming
          trim_poly_A: trim ending A/Ts on both sides, returns CO::FastQ::Read object

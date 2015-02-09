Trim paired-end FastQ files
===========================

## Overview:

`trimPair` is used to trim paired-end FastQ files. It wraps `cutadapt` to trim 
adapters and uses self-defined codes to remove poly-A/T/N on both ends of the 
reads.

`trimPair` initializes three child processes in which two processes trim reads in
both pairs, one processes do synchronization from two pairs and post processing.

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

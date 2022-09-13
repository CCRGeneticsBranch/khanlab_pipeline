#!/usr/bin/perl -w
use strict;
use Getopt::Long qw(GetOptions);
use Time::Piece;
use File::Basename;
use Cwd 'abs_path';


my $fastq_file;
my $usage = <<__EOUSAGE__;

Usage:

$0 [options]

required options:

  -i	FASTQ file
  
__EOUSAGE__



GetOptions (
  'i=s' => \$fastq_file
);


unless ($fastq_file) {
    die "Please input FASTQ file\n$usage";
}


&main();

sub main {
	open(INPUT_FILE, "$fastq_file") or die "Cannot open file $fastq_file";
	my $line_no = 0;
	while (<INPUT_FILE>) {
		if ($line_no % 4 == 1) {
			s/\#/N/g;
		}
		print;
		$line_no++;
	}
	close(INPUT_FILE)
}


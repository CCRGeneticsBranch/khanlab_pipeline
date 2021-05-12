#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;
use Cwd 'abs_path';

my $a;
my $b;
my $i = 1;
my $j = 1;
my $v = 0;
my $c = 0;

my $usage = <<__EOUSAGE__;

Usage:

$0 [options]

Options:

  -a  <string>  The file to be filtered
  -b  <string>  The file that contains the key column
  -i  <string>  The ith in a (default: $i)
  -j  <string>  The jth in b (default: $j)
  -c            Include columns in b
  -v            Remove data in b
  
  
__EOUSAGE__



GetOptions (
  'a=s' => \$a,
  'b=s' => \$b,
  'i=i' => \$i,
  'j=i' => \$j,
  'c' => \$c,
  'v' => \$v,
);

if (!$a || !$b) {
    die "Please specifiy the -a and -b!\n$usage";
}

open(B_FILE, "$b") or die "cannot open file $b"; 
my $header = <B_FILE>;
chomp $header;
my %keys = ();
while(<B_FILE>) {
	chomp;
	my @fields = split(/\t/);
	my $key = $fields[$j-1];
	$keys{$key} = $_;
}
close(B_FILE);


open(MAIN_FILE, "$a") or die "cannot open file $a"; 
my $main_header = <MAIN_FILE>;
chomp $main_header;
if ($c && !$v) {
	print "$main_header\t$header\n";
} else {
	print "$main_header\n";
}
while(<MAIN_FILE>) {
	chomp;
	my @fields = split(/\t/);
	my $key = $fields[$i-1];
	if (exists $keys{$key}) {
		if (!$v) {
			if ($c) {
				print "$_\t$keys{$key}\n";
			} else {
				print "$_\n";
			}
		}
	} else {
		if ($v) {
			print "$_\n";
		}
	}
}
close(MAIN_FILE);

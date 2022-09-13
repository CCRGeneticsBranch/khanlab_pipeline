#!/usr/bin/env perl

use strict;
use File::Basename;

if ($#ARGV < 2) {
	print("usage: perl mergeAquaStat.pl <sample_name> <human_stat.txt> <mouse_stat.txt> \n\n");
	print("Example:\nusage: perl mergeAquaStat.pl RH4_DMSO RH4_DMSO/mergeStats.txt mm10/RH4_DMSO/mergeStats.txt\n");
	exit(0);
}
my $sample = $ARGV[0];
my %data = ();
my @cols = ();
my @samples = ();

if ($#ARGV == 2) {
	push @samples, "${sample}.hg19";
	push @samples, "${sample}.mm10";
}

for (my $i=1;$i<=$#ARGV;$i++) {	
	my $file = $ARGV[$i];
	open(IN_FILE, $file);
	my $header = <IN_FILE>;
	if ($#ARGV > 2) {
		my ($empty, $s) = split(/\t/, $header);
		push @samples, $s;
	}
	while(<IN_FILE>) {
		chomp;
		my ($key, $value) = split(/\t/);
		push @{$data{$key}}, $value;
		push @cols, $key if ($i==1);
	}
	close(IN_FILE);
}

print "\t". join("\t", @samples)."\n";
foreach my $col(@cols) {
	print "$col\t".join("\t", @{$data{"$col"}})."\n";
}


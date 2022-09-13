#!/usr/bin/env perl

use strict;
use File::Basename;

if ($#ARGV < 2) {
	print("usage: perl mergeStatReplicates.pl <merged_sample_name> <rep1_mergeStats.txt rep2_mergeStats.txt...> \n\n");
	print("Example:\nusage: perl mergeStatReplicates.pl RH4_DMSO RH4_DMSO_rep1/mergeStats.txt RH4_DMSO_rep2/mergeStats.txt\n");
	exit(0);
}

my %data = ();
for (my $i=1;$i<=$#ARGV;$i++) {
	my $file = $ARGV[$i];
	open(IN_FILE, $file);
	<IN_FILE>;
	while(<IN_FILE>) {
		chomp;
		my ($key, $value) = split(/\t/);
		$data{$key} += $value;
	}
	close(IN_FILE);
}

my $sample = $ARGV[0];
print "\t$sample\n";
print "valid_interaction\t".$data{"valid_interaction"}."\n";
print "valid_interaction_rmdup\t".$data{"valid_interaction_rmdup"}."\n";
print "percent_valid_interaction_rmdup\t".round($data{"valid_interaction_rmdup"}/$data{"valid_interaction"}*100)."%\n";
print "trans_interaction\t".$data{"trans_interaction"}."\n";
print "cis_interaction\t".$data{"cis_interaction"}."\n";
print "percent_cis_interaction\t".round($data{"cis_interaction"}/$data{"valid_interaction_rmdup"}*100)."%\n";
print "cis_shortRange\t".$data{"cis_shortRange"}."\n";
print "cis_longRange\t".$data{"cis_longRange"}."\n";
print "percent_longRange_interaction\t".round($data{"cis_longRange"}/$data{"cis_interaction"}*100)."%\n";
close(IN_FILE);


sub round() {
	my ($value) = @_;
	return int($value * 1000)/1000;
}
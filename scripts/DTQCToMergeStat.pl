#!/usr/bin/env perl

use strict;
use File::Basename;
#example: for fn in */*PT.qc.txt;do bn=$(basename $fn);dn=$(dirname $fn);s=`echo $bn | sed 's/-PT\.qc\.txt//'`;perl /data/khanlab/projects/pipeline_production/khanlab_pipeline/scripts/DTQCToMergeStat.pl $fn $s $dn;done
open(IN_FILE, $ARGV[0]);
my $sample=$ARGV[1];
my $outdir=$ARGV[2];
if (!$sample) {
	$sample=basename($ARGV[0]);
	$sample =~ s/-PT\.qc\.txt//;
}
my %data = ();
open(OUT_DT, ">$outdir/${sample}-QC.txt");
print OUT_DT "\t$sample\tPercent\n";
while(<IN_FILE>) {
	chomp;
	my $key = substr($_, 0, 46);
	my $value = substr($_, 46, 11);
	my $pct = substr($_, 57);
	$key =~ s/\s+$//;
	$value =~ s/,//g;
	$data{$key} = $value;
	print OUT_DT "$key\t$value\t$pct\n";
}
close(OUT_DT);
open(OUT_STAT, ">$outdir/mergeStats.txt");
print OUT_STAT "\t$sample\n";
print OUT_STAT "valid_interaction\t".$data{"Mapped Read Pairs"}."\n";
print OUT_STAT "valid_interaction_rmdup\t".$data{"No-Dup Read Pairs"}."\n";
print OUT_STAT "percent_valid_interaction_rmdup\t".round($data{"No-Dup Read Pairs"}/$data{"Mapped Read Pairs"}*100)."%\n";
print OUT_STAT "trans_interaction\t".$data{"No-Dup Trans Read Pairs"}."\n";
print OUT_STAT "cis_interaction\t".$data{"No-Dup Cis Read Pairs"}."\n";
print OUT_STAT "percent_cis_interaction\t".round($data{"No-Dup Cis Read Pairs"}/$data{"No-Dup Read Pairs"}*100)."%\n";
print OUT_STAT "cis_shortRange\t".$data{"No-Dup Cis Read Pairs < 1kb"}."\n";
print OUT_STAT "cis_longRange\t".$data{"No-Dup Cis Read Pairs >= 1kb"}."\n";
print OUT_STAT "percent_longRange_interaction\t".round($data{"No-Dup Cis Read Pairs >= 1kb"}/$data{"No-Dup Cis Read Pairs"}*100)."%\n";
close(IN_FILE);
close(OUT_STAT);
sub round() {
	my ($value) = @_;
	return int($value * 1000)/1000;
}
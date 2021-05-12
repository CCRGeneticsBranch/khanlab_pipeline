#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;
use Cwd 'abs_path';

my $a;
my $b = "/data/khanlab/projects/pipeline_production/khanlab_pipeline/ref/gencode_v32_annotation.genes.txt";
my $i = 1;
my $j = 1;
my $o = 1;
my $s = 0;

my $usage = <<__EOUSAGE__;

Usage:

$0 [options]

Options:

  -a  <string>  The main annotation file
  -b  <string>  The file to be appended (default: $b)
  -i  <string>  The ith in a (default: $i)
  -j  <string>  The jth in b (default: $j)
  -o  <string>  The output column in b (default: $o)
  -s            a is sqanti output
  
__EOUSAGE__



GetOptions (
  'a=s' => \$a,
  'b=s' => \$b,
  'i=i' => \$i,
  'j=i' => \$j,
  'o=i' => \$o,
  's' => \$s,
);

if (!$a) {
    die "Please specifiy the main file!\n$usage";
}

open(SEC_FILE, "$b") or die "cannot open file $b"; 
my $header = <SEC_FILE>;
chomp $header;
my $empty_line = "";
my @headers = split(/\t/, $header);
@headers = @headers[($o-1) .. $#headers];
$header = join("\t", @headers);
foreach my $h (@headers) {
	$empty_line = $empty_line."\t.";
}
my %sec = ();
while(<SEC_FILE>) {
	chomp;
	my @fields = split(/\t/);
	my $out_line = join("\t", @fields[($o-1) .. $#fields]);
	my $key = $fields[$j-1];
	if ($s) {
		($key) = $key =~ /(PB\..*)\..*/
	}
	$sec{$key} = $out_line;
}
close(SEC_FILE);


open(MAIN_FILE, "$a") or die "cannot open file $a"; 
my $main_header = <MAIN_FILE>;
chomp $main_header;
print "$main_header\t$header\n";
while(<MAIN_FILE>) {
	chomp;
	my @fields = split(/\t/);
	my $ori_gene_id = $fields[$i-1];
	$ori_gene_id =~ s/"//g;
	my $gene_id = $ori_gene_id;
	if ($ori_gene_id =~ /^ENSG.*_.*/) {
		($gene_id) = $ori_gene_id =~ /^(ENSG.*)_.*/
	}
	if (exists $sec{$gene_id}) {
		print "$_\t$sec{$gene_id}\n";
	} else {
		print "$_\t$ori_gene_id\t$ori_gene_id\tnovel\t.\n";
	}	
}
close(MAIN_FILE);

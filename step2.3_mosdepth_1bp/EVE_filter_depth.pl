#!/usr/bin/perl -w
use strict;

die "perl $0 [ .regions.bed.gz ] [output file]\n" unless @ARGV == 2;

my $total = 0;
my $eve = 0;

open (IN, "gzip -dc $ARGV[0] | ") or die "$ARGV[0] $!\n";
open (OUT, ">$ARGV[1]") or die "$ARGV[1] $!\n";

while(<IN>){
	chomp;
	my @line = split;
	$total ++;

	if($line[-1] >= 5){
		$eve ++;
	}
}
close IN;

my $perc = $eve/$total;
print "Depth percentage: $perc\n";
print OUT "Depth percentage: $perc\n";

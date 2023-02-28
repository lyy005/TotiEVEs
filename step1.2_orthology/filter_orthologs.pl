#!/usr/bin/perl -w
use strict;

die "perl $0 [ Orthogroups.GeneCount.tsv ] [ threshold ] [ output file ]\n" unless @ARGV == 3;

open (LS, "$ARGV[0]") or die "$ARGV[0] $!\n";
open (OUT, ">$ARGV[2]") or die "$ARGV[2] $!\n";

my $threshold = $ARGV[1];

my $header = <LS>;
#print "$header\n";

my @header = split /\s+/, $header;
shift @header; 
pop @header;
my $sum_of_spp = @header;

print "#spp in total: $sum_of_spp\n"; 

while(<LS>){
	chomp;
	my @line = split;
	my $og = shift @line; 
	my $total = pop @line;

	my $flag = 1;

	foreach my $k (@line){
		if ($k > 1){
			$flag = 0;
		}
	}

	my $percentage = $total / $sum_of_spp;
	if($percentage >= $threshold){
	}else{
		$flag = 0;
	}

	if($flag){
		print OUT "$og\t$total\t$percentage\n"; 
	}
}
close LS;


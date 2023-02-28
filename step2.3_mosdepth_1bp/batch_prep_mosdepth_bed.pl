#!/usr/bin/perl -w
use strict;

die "perl $0 [sup5_EVEs_ls.blast.modified] [LsHau_GCA_014465815.1_genomic.fa.bed.1] [output prefix]\n" unless @ARGV == 3;

my $eve = 1;

open (LS, "$ARGV[0]") or die "$ARGV[0] $!\n";
while(<LS>){
	my @line = split;
	my $ref = $line[3];

	my $start;
	my $end;

	if($line[10] < $line[11]){
		$start = $line[10];
		$end = $line[11];
	}else{
		$start = $line[11];
		$end = $line[10];
	}

	print "EVE$eve\t$ref\t$start\t$end\n";

	open (IN, "$ARGV[1]") or die "$ARGV[1] $!\n";
	open (OUT, ">$ARGV[2].EVE$eve") or die "$ARGV[2].EVE$eve $!\n";
	while(<IN>){
		my @bed = split;
		if($bed[0] eq $ref){
			if(($bed[2] >= $start)&&($bed[2] <= $end)){
				print OUT;
			}
		}
	}
	close IN;
	close OUT;
	$eve ++;
}
close LS;


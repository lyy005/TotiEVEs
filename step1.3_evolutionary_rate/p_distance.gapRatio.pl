#!/usr/bin/perl -w
use strict;

die "perl $0 [aligned fasta] [output file] [gap ratio]\n" unless @ARGV == 3;

my %hash;

my $total_distance = 0;
my $count_distance = 0;
my $average_distance = 0;

my $gaps;
my $len;
my $gap_ratio = 0;

# read sequences
open (IN, "$ARGV[0]") or die "$ARGV[0] $!\n";
$/ = ">";
<IN>;
while(<IN>){
	chomp;
	my @line = split /\n+/;
	my $name = shift @line;
	my $seq = join "", @line;
	
	$hash{$name} = $seq;

	$len = length($seq);
	$gaps = ($seq =~ s/\-/\-/g);
	my $ratio_cal = $gaps/$len;

	if($ratio_cal >= $gap_ratio){
		$gap_ratio = $ratio_cal;
	}
}
close IN;

print "Gap ratio: $gap_ratio\n";
if($gap_ratio >= $ARGV[2]){
	print "Gap ratio too high, quit\n";
	open (OT, ">$ARGV[1].gappy") or die "$ARGV[1].gappy $!\n";
	print OT "Gap ratio = $gap_ratio\n";
	close OT;
	exit;
}else{

open (OT, ">$ARGV[1]") or die "$ARGV[1] $!\n";
my @species = sort keys %hash;

# calculate pairwise distance
for my $i (0..$#species){
	for my $j ($i+1..$#species){

		my $seq1 = $hash{$species[$i]};
		my $seq2 = $hash{$species[$j]};

		my $count = ( $seq1 ^ $seq2 ) =~ tr/\0//c;
		my $len = length($hash{$species[$i]}); 
		my $percent = $count/$len;

		$total_distance += $percent;
		$count_distance ++;

		print OT "$species[$i]\t$species[$j]\t$percent\n";
	}
}

$average_distance = $total_distance/$count_distance;
print OT "Average distance: $average_distance\n";
print OT "Gap ratio = $gap_ratio\n";

print "output file: $ARGV[1]\n";
}

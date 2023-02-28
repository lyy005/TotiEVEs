#!/usr/bin/perl -w
use strict;

die "Usage: perl $0 [ file list (Orthogroups.GeneCount.tsv.0.95) ] [ dir to fasta files]\n" unless (@ARGV == 2);
open LST,$ARGV[0] or die "$!\n";
#open OUT,">$ARGV[2]" or die "$!\n";

my $dir = $ARGV[1];

my %hash;

while(<LST>){
	chomp;
	my $line = $_;
	my @line = split /\s+/, $line;
	my $i = $line[0];

# calculate p-distance and average distance of pairwise comparisions of the three species for each OG
	#print "Calculate p-distance for $i\.GBlocks.fas\n";
	#print "perl p_distance.pl $i\.GBlocks.fas $i\.pairwise_distance\n";
	#`perl /share/home/liyiyuan/projects/collab_0_EVE/sample2_Soga/step4_orthoFinder/run1_iqtree/p_distance.pl $dir\/$i\.GBlocks.fas $dir\/$i\.pairwise_distance`;
	
	print "Calculate p-distance for $i\.fa.align\n";
	print "perl p_distance.gapRatio.pl $i\.fa.align $i\.pairwise_distance 0.1\n";
	`perl /share/home/liyiyuan/projects/collab_0_EVE/sample2_Soga/step4b_orthoFinder_FeiLi/run2_distance/p_distance.gapRatio.pl $dir\/$i\.fa.align $dir\/$i\.pairwise_distance 0.1`;
}
close LST;

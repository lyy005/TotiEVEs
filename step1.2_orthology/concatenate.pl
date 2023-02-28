#!/usr/bin/perl -w
use strict;

die "Usage: perl $0 [ file list (Orthogroups.GeneCount.tsv.0.95) ] [ dir to fasta files] [ output ]\n" unless (@ARGV == 3);
open LST,$ARGV[0] or die "$!\n";
open OUT,">$ARGV[2]" or die "$!\n";

my $dir = $ARGV[1];

my %hash;

while(<LST>){
	chomp;
	my $line = $_;
	my @line = split /\s+/, $line;
	my $i = $line[0];
# run alignment
	print "Step 1: running alignment\n";

	if (-e "$dir\/$i\.GBlocks.fas") {
		print "$i\.GBlocks.fas exists\n";
	}else{
		print "mafft-linsi $dir\/$i\.fa 1> $dir\/$i\.fa\.align 2> $dir\/$i\.fa\.align.log\n";
		`mafft-linsi $dir\/$i\.fa 1> $dir\/$i\.fa\.align 2> $dir\/$i\.fa\.align.log`;

		print "Gblocks $dir\/$i\.fa\.align -t=p -b4=10 -b5=n -e=-gb > $dir\/$i\.GBlocks.log\n";
		`Gblocks $dir\/$i\.fa\.align -t=p -b4=10 -b5=n -e=-gb > $dir\/$i\.GBlocks.log`;

		print "perl gblock_convert.pl $dir\/$i\.fa\.align-gb > $dir\/$i\.GBlocks.fas\n";
		`perl gblock_convert.pl $dir\/$i\.fa\.align-gb > $dir\/$i\.GBlocks.fas`;
	}
	print "Step 2: concatenate\n";
	open FA, "$dir\/$i\.GBlocks.fas" or die "$!\n";
	$/=">";
	<FA>;
	while(<FA>){
		chomp;
		my $line = $_;
		my @line = split /\n+/, $line;
		my $name = shift @line;

		my @name = split /\_/, $name;
#		pop @name; 
#		my $id = $name[0]."_".$name[1]."_".$name[2];
		my $id = $name[0];
		my $seq = join "", @line;

		$hash{$id} .= $seq;
	}
	close FA;
	$/="\n";
		
}
close LST;

foreach my $k (sort keys %hash){
	print OUT ">$k\n$hash{$k}\n";
}

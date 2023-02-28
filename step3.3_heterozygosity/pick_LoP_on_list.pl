#!/usr/bin/perl -w
use strict;

die "perl $0 [LoP] [EVE1.summary] [coverage threshold] [output LoP]\n" unless @ARGV == 4;

my %hash;
my $threshold = $ARGV[2];

open (LS, "$ARGV[1]") or die "$ARGV[0] $!\n";
while(<LS>){
	chomp;
	s/\>//;
	my @line = split;

	if($line[1] >= $threshold){
		$hash{$line[0]} = 1;
	}
}
close LS;

open (IN, "$ARGV[0]") or die "$ARGV[0] $!\n";
open (OUT, ">$ARGV[3]") or die "$ARGV[3] $!\n";
while(<IN>){
	chomp;
	my $all = $_;
	my @line = split /\s+/;
	my $id = shift @line;
	my $seq = join "", @line;
#	my @names = split /\_/, $line[0];
#	if($hash{$names}){
	my @id = split /\s+/, $id;

	if($hash{$id[0]}){
#		print OUT ">$id\n$seq\n";
		print OUT "$all\n";
	}

}
close IN;
print "Done!\nOutput: $ARGV[3]\n";

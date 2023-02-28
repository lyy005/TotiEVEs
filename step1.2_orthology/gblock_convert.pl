#!/usr/bin/perl -w
use strict;
 
die "Usage: perl $0 [ -gb file ]\n" unless (@ARGV == 1);
open FA,$ARGV[0] or die "$!\n";

$/ = ">"; 
<FA>; 
while(<FA>){
	chomp;
	my @line = split/\n+/; 
	my $name = shift @line; 
	my @name=split/\s+/, $name; 
	my $seq = join "",@line; 
	$seq =~s/\s+//g; 
	print ">$name[0]\n$seq\n";
}
close FA;

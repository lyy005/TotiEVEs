#!/usr/bin/perl -w
use strict;

die "usage: perl $0 [braker.pep.fasta] [output]\n" unless (@ARGV == 2);
open IN, $ARGV[0] or die "$!\n";
open OUT, ">$ARGV[1]" or die "$!\n";

# input format: 
# g1.p1
# g23.p1
# g23.p2

# >Sfur000016.1
# >Sfur000016.2
#
#
my %gene;
my %found_genes;

my (%len, %pick);
$/=">";
<IN>;
while(<IN>) {
	chomp;
	my $line = $_;
	my @line = split/\n+/, $line;

	my $name = shift @line;
	my @name = split /\.+/, $name;
	#my @gene_id = split /\-/, $name[0];
	
	# calculate the length of the protein
	my $seq = join "", @line;
	my $len = length ($seq);

	# save the length of the protein under the gene ID; only save the longest protein ID for each gene
	# save the protein ID in a hash called pick
	my $gene_id;

	$gene_id = $name[0];

	if(defined $len{$gene_id}){
		if($len > $len{$gene_id}){
			$len{$gene_id} = $len;
			$pick{$gene_id} = $name;
		}
	}else{
		$len{$gene_id} = $len;
		$pick{$gene_id} = $name;
	}
}	
close IN;

my @picks = values %pick;
open IN,$ARGV[0] or die "$!\n";
$/=">";
while(<IN>) {
	chomp;
	my $line = $_;
	my @line = split/\n+/, $line;
	
	my $name = shift @line;
	my $seq = join "", @line;

#	if($name ~~ [values %pick]){
	if($name ~~ @picks){
		print OUT ">$name\n$seq\n";
	}
}
close IN;

print "Done\n";

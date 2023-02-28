#!/usr/bin/perl -w
use strict;

die "perl $0 [BLAST output] [vcf file] [LoP output file]\n" unless @ARGV == 3;

my (%hash_len, %hash_snp);

open (IN, "$ARGV[0]") or die "$ARGV[0] $!\n";
open (VCF, "$ARGV[1]") or die "$ARGV[1] $!\n";
open (OT, ">$ARGV[2]") or die "$ARGV[2] $!\n";
open (TMP, ">$ARGV[2].tmp") or die "$ARGV[2].tmp $!\n";

my ($start, $end);
my $total_len = 0;
my $total_hetero = 0;
my %hetero;
my @sample_header;

# BLAST input example: 
# 975	68309712	Sfur012369.1	CM025290.1	99.803	508	1	0	335	842	19030459	19031982	0.0	614
while(<IN>){
	chomp;
	my @blast_line = split /\s+/;
	
	if($blast_line[10] > $blast_line[11]){
		$start = $blast_line[11];
		$end = $blast_line[10];
	}else{
		$start = $blast_line[10];
		$end = $blast_line[11];
	}

	my $len = $end - $start + 1;
	$total_len += $len;
 
	# my %hetero;
print TMP "BLAST Query $blast_line[2] BLAST db chromosome: $blast_line[3] From $start To $end\n";
	
	open (VCF, "$ARGV[1]") or die "$ARGV[1] $!\n";
	while(<VCF>){
		if(/\#CHROM/){
			my $header = $_;
			my @header = split /\s+/, $header;
			@sample_header = @header[9..$#header];
		}elsif(/^\#/){
			next;
		}

		my @vcf_line = split /\s+/;
		#print "### $vcf_line[0]\n";

		if($vcf_line[0] eq $blast_line[3]){
			if(($vcf_line[1] >= $start)&&($vcf_line[1] <= $end)){

print TMP "@vcf_line\n";

				my @samples = @vcf_line[9..$#vcf_line];

				foreach my $k (0..$#samples){
					my @sample_line = split/\:/, $samples[$k];
					$hetero{$sample_header[$k]} ++ if($sample_line[0] eq "0/1");
					$hetero{$sample_header[$k]} ++ if($sample_line[0] eq "0|1");
					$hetero{$sample_header[$k]} ++ if($sample_line[0] eq "1|0");
					$hetero{$sample_header[$k]} ++ if($sample_line[0] eq "1/0");
				}
			}
		}
		
	}
	close VCF;
}
close IN;

print OT "samples\thetero-sites\ttotal-sites\thetero-ratio\n";
foreach my $l (@sample_header){
	if($hetero{$l}){
		my $hetero_ratio = $hetero{$l} / $total_len;
		#print OT "$l: hetero-sites: $hetero{$l}\t# total sites: $total_len\t# hetero-ratio: $hetero_ratio\n";
		print OT "$l\t$hetero{$l}\t$total_len\t$hetero_ratio\n";
	}else{
		#print OT "$l: hetero-sites: 0\t# total sites: $total_len\t# hetero-ratio: 0\n";
		print OT "$l\t0\t$total_len\t0\n";
	}
}
print "Output: $ARGV[2]\n";

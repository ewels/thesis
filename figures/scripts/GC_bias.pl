#!/usr/bin/perl
use warnings;
use strict;
use Math::Round;
#############################################################
# Name: GC Bias                                             #
# Author: Phil Ewels                                        #
# Version 1.0 - 04/11/2011                                  #
#############################################################

# Set input file by command line
my ($input) = @ARGV;
if (!defined $input) {
	die "Usage is GC_bias.pl [input file]\n";
}

# Get genome into a hash
warn "Loading Genome...\n";
my $fn_base = 'D:\Genome Sequences\Human\chr%s.fa'; 
my @chromosomes = (1..21,'X','Y');
my %chrom;
foreach my $chromosome (@chromosomes) {
	my $filename = sprintf($fn_base, $chromosome);
	open (IN,$filename) or die "Can't read file: $!";
	$_ = <IN>; # Remove fasta header
	while (my $line = <IN>) {
		chomp ($line);
		$chrom{$chromosome} .= uc($line); # Make everything upper case
	}
	warn "Chromosome $chromosome loaded...\n";
}

# Load annotated report
open (IN,$input) or die "Can't read file: $!";
$_ = <IN>; # File Header

my %all_frags;
my %hit_frags;

while (my $line = <IN>) {
	chomp ($line);
	my @field = split(/\t/, $line);
	# $field[1] = Chromosome
	# $field[2] = Start
	# $field[3] = End
	if (exists $chrom{$field[1]}) { # Check that we have the chromosome for this fragment
		my $seq = substr($chrom{$field[1]},$field[2],($field[3] - $field[2] + 1));
		my $at = ($seq =~ tr/AaTt//);
		my $gc = ($seq =~ tr/GgCc//);
		$at = (int($at)) ? $at : 0; # make sure that vars are numeric
		$gc = (int($gc)) ? $gc : 0;
#		print "$seq\n$at\t$gc\n\n"; sleep(1);
		if(($gc + $at) > 0){ # kick out empty strings
			my $ratio = $gc / ($gc + $at);
#			print "$at\t$gc\t$ratio\n";
			$ratio = nearest(0.05, $ratio); # rounds to the nearest 5%
#			print "$at\t$gc\t$ratio\n\n"; sleep(1);
			$all_frags{$ratio} += 1; # count for all fragments
			if($field[11] eq '1.0') { # count only for hit hits
				$hit_frags{$ratio} += 1;
			}
		}
	}
}

my $all_num_frags = 0;
my $hit_num_frags = 0;
while ((my $key, my $value) = each(%all_frags)){
	$all_num_frags += $value;
}
while ((my $key, my $value) = each(%hit_frags)){
	$hit_num_frags += $value;
}
print "\n\n$all_num_frags\t$hit_num_frags\n\n";

# Open output file
open (OUT, '>',$input."_GCbias_output.txt") or die "Can't read file: $!";

print OUT "% GC\tAll Frags\tHit Frags\t% All \t% Hit \t% Hit / % All\n";
foreach my $key (sort keys %all_frags) {
	my $percent_all = 0;
	my $percent_hit = 0;
	my $num_hit = 0;
	my $gc_percent = $key * 100;
	$percent_all = $all_frags{$key}/$all_num_frags;
	if (exists $hit_frags{$key}) {
		$num_hit = $hit_frags{$key};
		$percent_hit = $hit_frags{$key}/$hit_num_frags;
	}
	my $hit_all = $percent_hit / $percent_all;
	$hit_all = ($hit_all == 0) ? 1 : $hit_all;
	$percent_all = sprintf("%.5f", $p
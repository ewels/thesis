#!/usr/bin/perl
use warnings;
use strict;
#############################################################
# Name: Single RE Fragments                                 #
# Author: Phil Ewels                                        #
# Version 1.0 - 05/05/2011                                  #
# --------------------------------------------------------- #
#############################################################

#===============================
#==== CONFIGURATION OPTIONS ====
#===============================

# Set by command line (optional)
my ($output,$re_search) = @ARGV;
if (defined $output && !defined $re_search) {
	die "Usage is single_RE_fragments.pl [output file] [search string]\nLeave blank to use defaults\n";
} elsif(!defined $output) {
	# Path to output file
	$output = 'AseI_fragments.txt';
	# Restriction site to use
	$re_search = 'ATTAAT';
	warn "Using file defaults: search string = $re_search, output = $output\n";
} else {
	warn "Using command line variables: search string = $re_search, output = $output\n";
}

# Path to chromosome fasta files. Replace chromosome number with %s
my $fn_base = 'D:\Genome Sequences\Human\chr%s.fa'; 
#$fn_base = 'D:\Genome Sequences\Human\GRCh37\Homo_sapiens.GRCh37.55.dna.chromosome.%s.fa';
warn "Looking for Genome Sequences in $fn_base\n\n";

# Chromosomes to use. Default is (1..22,'X','Y') - other options are 'MT' etc.
my @chromosomes = (1..22,'X','Y');


#==================================
#== END OF CONFIGURATION OPTIONS ==
#==================================

open (OUT,'>',$output) or die $!;

# go through each chromosome
foreach my $chromosome (@chromosomes) {
	my $filename = sprintf($fn_base, $chromosome);
	open (IN,$filename) or die "Can't read file: $!";
	warn "Starting Chromosome $chromosome ($filename)\n";
	my $sequence = '';
	$_ = <IN>; # Remove fasta header
	while (my $line = <IN>) {
		chomp ($line);
		$sequence .= uc($line); # Make everything upper case
	}
	
	my ($offset, $lastpos, $pm) = (0, 1, "+");
	my $pos = index($sequence, $re_search, $offset);
	while ($pos != -1) {
		# Kick out fragments with lots of N's (centromeres and telomeres)
		#next if ( index( substr($sequence,$lastpos,($pos-$lastpos)),'NNNNN') >=0);
		# Too slow...
	
		print OUT $chromosome."\t". # Chromosome Name
				 $lastpos."\t". # Position Start
				 $pos."\t". # Position Finish
				 $pm."\n"; # Arbitrary orientation to help visulaisation
		if ($pm eq "+") { $pm = "-"; } else { $pm = "+"; }
		$lastpos = $pos + length($re_search);
		$offset = $pos + length($re_search);
		$pos = index($sequence, $re_search, $offset);
	}

	close IN;
}
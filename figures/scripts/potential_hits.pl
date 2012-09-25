#!/usr/bin/perl
use warnings;
use strict;
#############################################################
# Name: Potential Hits (4C)                                 #
# Authors: Phil Ewels                                       #
# Version 1.0 - 08/04/2011                                  #
# --------------------------------------------------------- #
#############################################################

###########################
## CONFIGURATION OPTIONS ##
###########################

# Path to chromosome fasta files. Replace chromosome number with %s
my $fn_base = 'D:\Genome Sequences\Human\chromFa\chr%s.fa'; 
#my $fn_base = 'D:\Genome Sequences\Human\GRCh37\Homo_sapiens.GRCh37.55.dna.chromosome.%s.fa';

# Path to output file
my $output = 'AseI_NlaIII_fragments.txt';

# Chromosomes to use. Default is (1..22,'X','Y') - other options are 'MT' etc.
my @chromosomes = (1..22,'X','Y');

# First restriction site to use (3C digestion - assumes palindromic 6 cutter)
my $re_search1 = 'ATTAAT';

# Second restriction site to use (4C digestion)
my $re_search2 = 'CATG';

# Size selection parameters (in base pairs)
my $minsize = 10;
my $maxsize = 700;

##################################
## END OF CONFIGURATION OPTIONS ##
##################################

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
		$sequence .= $line;
	}

	# Search for restriction enzyme sites (greedy regex)
	while ($sequence =~ /$re_search1(.*?)$re_search1/g) {
		# Kick out fragments with lots of N's (centromeres and telomeres)
		next if (index ($1,'NNNNN') >=0);
		# Add fragment to array, with half of restriction site added back on either end
		# Not using real site so we can revcomp and join without modifying (if AseI = AT^TAAT, fragment is AAT...ATT)
		my $fragment = substr($re_search1,2,4).$1.substr($re_search1,0,2);
		
		# Only proceed if we have a NlaIII site (otherwise ignore fragment)
		if ($fragment =~ /CATG/i) {
			
			# Get sequence from start of string to first NlaIII site
			$fragment =~ /(^.+CATG?).*$/i;
			#print OUT "$chromosome_name\t$start\t$end\t$length_fragment\t$string\n";
			#warn "$chromosome\t".pos($sequence)."\t".(pos($sequence)+length("${1}CATG"))."\t+\t".length("${1}CATG")."\t${1}CATG\n";
			if(length($fragment) > $minsize && length($fragment) < $maxsize) {
				print OUT "$chromosome\t".pos($sequence)."\t".(pos($sequence)+length("${1}CATG"))."\t+\t".length("${1}CATG")."\t${1}CATG\n";
			}
			
			#  Get sequence from last NlaIII site to end of string
			$fragment =~ /^.+CATG(.+$)/i;
			#warn "$chromosome\t".pos($sequence)."\t".(pos($sequence)+length("${1}CATG"))."\t-\t".length("${1}CATG")."\t${1}CATG\n";
			#sleep(1);
			if(length($fragment) > $minsize && length($fragment) < $maxsize) {
				print OUT "$chromosome\t".pos($sequence)."\t".(pos($sequence)+length("${1}CATG"))."\t-\t".length("${1}CATG")."\t${1}CATG\n";
			}		
		}
		# Move array pointer back six so that we don't skip a fragment
		pos($sequence) -= 6;
	}
	close IN;
}
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










#	warn "\n";
#	open (OUT,'>',$output) or die $!;
#	my %concatamer_stats;
#	my %sequenced_concatamers;
#	my $fragment_count = 0;
#	while ($fragment_count < $library_size) {
#	
#	    my $concatamer = '';
#		my $fragment = '';
#		my $concatamer_count = 0;
#		
#		# Select a starting fragment
#		if($fourc) {
#			# For 4C - predetermined starting fragment
#			$concatamer = $fourc_bait;
#			++$concatamer_count;
#		} else {
#			# For Hi-C (everything to everything)
#			while (1) {
#				# pick random starting fragment
#				$fragment = $fragments[int rand(@fragments)];
#				
#				# if we don't have an NlaIII site choose another starting fragment
#				next unless ($fragment =~ /CATG/i);
#				
#				# randomly reverse complement
#				if (int rand(2)) {
#					$fragment = reverse($fragment);
#					$fragment =~ tr/GATCgatc/CTAGctag/;
#				}
#				
#				$concatamer = $fragment;
#				
#				# Strip sequence from start of fragment to first NlaIII site
#				$concatamer =~ s/^\w+CATG//i;
#				die "Nothing removed from $concatamer" if ($concatamer eq $fragment);
#				die "Didn't end with ATT" unless ($concatamer =~ /ATT$/i);
#				++$concatamer_count;
#				last;
#			}
#		}
#		
#		# Keep adding more framents until we hit another NlaIII site
#		while (1) {
#	
#			# Pick second fragment randomly
#			my $fragment = $fragments[int rand(@fragments)];
#	
#			# Randomly reverse complement
#			if (int rand(2)) {
#				$fragment = reverse($fragment);
#				$fragment =~ tr/GATCgatc/CTAGctag/;
#			}
#	
#			# Only proceed if we have a second NlaIII site (otherwise add another fragment)
#			if ($fragment =~ /CATG/i) {
#				
#				# Strip sequence from first NlaIII site to end of fragment
#				$fragment =~ s/(^.*?CATG).*$/$1/i;
#	
#				# Add second fragment to concatamer
#				$concatamer .= $fragment;
#	
#				unless ($concatamer =~ /$re_search1/i) {
#					die "No restriction site $re_search1 from last addition. Fragment:\n$fragment\n\nConcatamer:\n$concatamer\n\n";
#				}
#				
#				# Only use this finished fragment if we satisfy the size selection
#				if(length($concatamer) > $minsize && length($concatamer) < $maxsize) {
#					# Print fake sequencing result from fragment to file 
#					print OUT substr($concatamer,$seqtrim_start,$readlength),"\n";
#					++$concatamer_count;
#					++$concatamer_stats{$concatamer_count};
#					my $seq_concatamers = 0;
#					$seq_concatamers++ while substr($concatamer,$seqtrim_start,$readlength) =~ /$re_search1/g;
#					++$sequenced_concatamers{$seq_concatamers};
#					++$fragment_count;
#					warn "Fragment count = $fragment_count\n" if ($fragment_count % (int($library_size)/10) == 0);
#				}
#	#			else {
#	#				warn "Rejecting concatamer with length ".length($concatamer)."\n";
#	#				sleep 1;
#	#			}
#				last;
#			}
#			
#			# No NlaIII site found, add fragment to concatamer and loop again
#			else {
#				++$concatamer_count;
#				$concatamer .= $fragment;
#			}
#	
#		} # finished this fragment
#	
#	} # finished all fragments
#	
#	# Concatamer Stats
#	warn "\nTotal Concatamer Statistics before trimming for Sequencing:\n";
#	foreach my $key (sort keys %concatamer_stats) {
#		warn "\t$key fragments = $concatamer_stats{$key}\t\t(".sprintf("%.6f",(($concatamer_stats{$key}/$library_size)*100))."%)\n";
#	}
#	
#	warn "\nNumber of Restriction Enzyme sites after trimming for Sequencing:\n";
#	foreach my $key (sort keys %sequenced_concatamers) {
#		warn "\t$key fragments = $sequenced_concatamers{$key}\t\t(".sprintf("%.6f",(($sequenced_concatamers{$key}/$library_size)*100))."%)\n";
#	}
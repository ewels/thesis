#!/usr/bin/perl
use warnings;
use strict;
use Math::Round;
#############################################################
# Name: GC Bias                                             #
# Author: Phil Ewels                                        #
# Version 1.0 - 04/11/2011                                  #
#############################################################

####################
####   SETUP    ####
####################
my $GC_binsize = 5/100;
my $AN_fraglength_binsize = 10;



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
	open (CHROMOSOMES,$filename) or die "Can't read file: $!";
	$_ = <CHROMOSOMES>; # Remove fasta header
	while (my $line = <CHROMOSOMES>) {
		chomp ($line);
		$chrom{$chromosome} .= uc($line); # Make everything upper case
	}
	warn "Chromosome $chromosome loaded...\n";
}

# Load annotated report
open (IN,$input) or die "Can't read file: $!";
$_ = <IN>; # File Header
my %all_frags_percent;
my %hit_frags_percent;
my %all_frags_length;
my %hit_frags_length;

##########################################
####   WORK OUT CORRECTION FACTORS    ####
##########################################
warn "\nCalculating statistics";
my $counter = 0;
while (my $line = <IN>) {
	$counter++;
	chomp ($line);
	my @field = split(/\t/, $line);
	my $chr = $field[1]; # Chromosome
	my $start = $field[2]; # Start
	my $end = $field[3]; # End
	my $hit = ($field[11] eq '1.0') ? 1 : 0;
	if (exists $chrom{$chr}) { # Check that we have the chromosome for this fragment
		########### GC CONTENT
		my $seq = substr($chrom{$chr},$start,($end - $start + 1));
		my $at = ($seq =~ tr/AaTt//);
		my $gc = ($seq =~ tr/GgCc//);
		$at = (int($at)) ? $at : 0; # make sure that vars are numeric
		$gc = (int($gc)) ? $gc : 0;
		if(($gc + $at) > 0){ # kick out empty strings
			my $ratio = $gc / ($gc + $at);
			$ratio = nearest($GC_binsize, $ratio);
			$all_frags_percent{$ratio} += 1; # count for all fragments
			if($hit) { # count only for hit hits
				$hit_frags_percent{$ratio} += 1;
			}
		}
		
		######### FRAGMENT LENGTH
		my $length = $end - $start + 1; # +1 because of genomic co-ordinates, eg. pos 4 to 5 is 2 bases, 5-4 = 1
		$length = nearest($AN_fraglength_binsize, $length);
		$all_frags_length{$length} += 1; # count for all fragments
		if($hit) { # count only for hit hits
			$hit_frags_length{$length} += 1;
		}
	}
	if($counter % 200000 == 0) { warn "$counter lines analysed\n"; }
}
warn "\nStatistics calculated.";

##### Total number of fragments processed
my $all_num_frags = 0;
my $hit_num_frags = 0;
while ((my $key, my $value) = each(%all_frags_percent)){
	$all_num_frags += $value;
}
while ((my $key, my $value) = each(%hit_frags_percent)){
	$hit_num_frags += $value;
}
warn "\n\nNumber of Fragments Processed\nTotal: $all_num_frags\tHit: $hit_num_frags\n\n";

##### Calculate GC correction values
warn "Caculating GC correction values\n";
my %GC_correction;
foreach my $key (sort keys %all_frags_percent) {
	my $percent_all = 0;
	my $percent_hit = 0;
	$GC_correction{$key} = 1;
	if (exists $hit_frags_percent{$key}) {
		$percent_all = $all_frags_percent{$key} / $all_num_frags;
		$percent_hit = $hit_frags_percent{$key} / $hit_num_frags;
		$GC_correction{$key} = $percent_hit / $percent_all;
	}
}
##### Optional - print output GC correction values
warn "Printing GC correction values\n";
open (GC_OUT, '>',$input."_GC_correction_factors.txt") or die "Can't read file: $!";
foreach my $key (sort keys %GC_correction) {
	print GC_OUT ($key*100)."%\t".$GC_correction{$key}."\n";
}

##### Calculate AseI - NlaIII fragment length correction values
warn "Caculating AseI - NlaIII fragment length correction values\n";
my %ANlength_correction;
foreach my $key (sort {$a<=>$b} keys %all_frags_length) {
	my $percent_all = 0;
	my $percent_hit = 0;
	$ANlength_correction{$key} = 1;
	if (exists $hit_frags_length{$key}) {
		$percent_all = $all_frags_length{$key} / $all_num_frags;
		$percent_hit = $hit_frags_length{$key} / $hit_num_frags;
		$ANlength_correction{$key} = $percent_hit / $percent_all;
	}
}

##### Optional - print output frag length correction values
warn "Printing AseI - NlaIII fragment length correction values\n";
open (FL_OUT, '>',$input."_fragLength_correction_factors.txt") or die "Can't read file: $!";
foreach my $key (sort {$a<=>$b} keys %ANlength_correction) {
	print FL_OUT "$key\t".$ANlength_correction{$key}."\n";
}

#########################################
##### OUTPUT
#########################################
##### Print final correction value for each fragment
warn "\nPrinting final correction value file\n\n";
my $starting_chance = $hit_num_frags / $all_num_frags;
open (OUT, '>',$input."_correction_values.txt") or die "Can't read file: $!";
# reset input file pointer
seek(IN, 0, 0);
# loop through input file again
$counter = 0;
while (my $line = <IN>) {
	$counter++;
	chomp ($line);
	my @field = split(/\t/, $line);
	my $chr = $field[1]; # Chromosome
	my $start = $field[2]; # Start
	my $end = $field[3]; # End
	my $chance = $starting_chance;
	if (exists $chrom{$chr}) { # Check that we have the chromosome for this fragment
		########### GC CONTENT
		my $seq = substr($chrom{$chr},$start,($end - $start + 1));
		my $at = ($seq =~ tr/AaTt//);
		my $gc = ($seq =~ tr/GgCc//);
		$at = (int($at)) ? $at : 0; # make sure that vars are numeric
		$gc = (int($gc)) ? $gc : 0;
		if(($gc + $at) > 0){ # kick out empty strings
			my $ratio = $gc / ($gc + $at);
			$ratio = nearest($GC_binsize, $ratio);
			# Multiple the chance by the correction value for this ratio
			$chance = $chance * $GC_correction{$ratio}
		}
		
		######### FRAGMENT LENGTH
		my $length = $end - $start + 1; # +1 because of genomic co-ordinates, eg. pos 4 to 5 is 2 bases, 5-4 = 1
		$length = nearest($AN_fraglength_binsize, $length);
		# Multiple the chance by the correction value for
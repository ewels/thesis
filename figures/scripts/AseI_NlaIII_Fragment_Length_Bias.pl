#!/usr/bin/perl
use warnings;
use strict;
use Math::Round;
#############################################################
# Name: AseI - NlaIII Fragment Length Bias                  #
# Author: Phil Ewels                                        #
# Version 1.0 - 04/11/2011                                  #
#############################################################

# Set input file by command line
my ($input) = @ARGV;
if (!defined $input) {
	die "Usage is AseI_NlaIII_Fragment_Length_Bias.pl [input file]\n";
}

# Load annotated report
open (IN,$input) or die "Can't read file: $!";
$_ = <IN>; # File Header

my %all_frags;
my %hit_frags;

my $count = 0;
while (my $line = <IN>) {
	chomp ($line);
	my @field = split(/\t/, $line);
	# $field[1] = Chromosome
	# $field[2] = Start
	# $field[3] = End
	my $length = $field[3] - $field[2] + 1; # +1 because of genomic co-ordinates, eg. pos 4 to 5 is 2 bases, 5-4 = 1
	$length -= 35; # minumum cutoff was 35, so take this off here (and add on later)
	if($length < 0 ){
		next; # shouldn't have any less than 35??
	}
	$length = nearest(10, $length); # rounds to the nearest 10bp
	$all_frags{($length + 35)} += 1; # count for all fragments
	if($field[11] eq '1.0') { # count only for hit hits
		$hit_frags{$length} += 1;
	}
	$count++;
	if($count % 500000 == 0) { warn "$count lines processed..\n"; }
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
open (OUT, '>',$input."_AseI_NlaIII_frag_length.txt") or die "Can't read file: $!";

print OUT "Fragment Length\tAll Frags\tHit Frags\t% All\t% Hit\t% Hit / % All\n";
foreach my $key (sort {$a<=>$b} keys %all_frags) {
	my $percent_all = 0;
	my $percent_hit = 0;
	my $num_hit = 0;
	my $frag_length = $key;
	$percent_all = $all_frags{$key}/$all_num_frags;
	if (exists $hit_frags{$key}) {
		$num_hit = $hit_frags{$key};
		$percent_hit = $hit_frags{$key}/$hit_num_
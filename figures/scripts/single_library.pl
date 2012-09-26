#!/usr/bin/perl
use warnings;
use strict;
use Math::Random::MT::Perl qw(srand rand);
use IO::File;
#############################################################
# Name: Single Window p Value                               #
# Author: Phil Ewels                                        #
# Version 1.0 - 04/11/2011                                  #
#############################################################

####################
####   SETUP    ####
####################
my $repeats = 5;

# Set input file by command line
my ($corrections_file) = @ARGV;
if (!defined $corrections_file) {
	die "Usage is single_window.pl [corrections file]\n";
}

# Open output filehandles
my @output;
for my $i (0 .. ($repeats-1)) {
	push(@output, IO::File->new('> library_'.($i+1).'.txt'));
}

# Load corrections
warn "Loading correction factors..\n\n";
my $counter = 0;
open (CORRECTIONS, $corrections_file) or die "Can't read file: $!";
while (my $line = <CORRECTIONS>) {
	$counter++;
	chomp ($line);
	my @field = split(/\t/, $line)
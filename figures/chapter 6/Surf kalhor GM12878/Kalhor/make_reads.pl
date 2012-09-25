#!/usr/bin/perl -w
use strict;
use warnings;
#############################################################################
## Quantitative Data to Reads - takes a file with Chromosome co-ordinates  ##
## and an associated value and replicates the line according to its value  ##
## Allows data to be imported into SeqMonk for analysis next to data.      ##
##                                                                         ##
## If left blank, min max and multiplier will use defaults of 0 1 and 1000 ##
## This is designed to run on output files generated by make_windows.pl    ##
#############################################################################
## Written by Phil Ewels, 01/07/10                                         ##
## Last updated 13/12/10                                                   ##
#############################################################################
## Changelog:                                                              ##
##    13/12/10 - removed excessively complicated menu system               ##
#############################################################################

my ($infile,$outfile,$max,$min,$multiplier) = @ARGV;

unless (defined $multiplier) {
  die "Usage is make_reads.pl [input file] [output file] [max] [min] [multiplier]\nTypically, max = 1, min = 0 and multiplier = 1000\n\n";
}

unless ($max > $min) {
  die "Maximum value must be more than Minimum value";
}

unless ($multiplier > 0) {
  die "Multiplier must be greater than 0";
}

open (INPUT,$infile) or die "Can't read $infile: $!";

if (-e $outfile) {
  print "$outfile already exists.  Overwrite [Y/N]?";
  my $answer = <STDIN>;
  if ($answer !~ /^y/i) {
    die "Exiting...\n";
  }
}

open (OUTPUT,'>',$outfile) or die "Can't write to $outfile: $!";




## DEFINE VARIABLES
my $loop = "true";
my $loop2 = "true";
my $i = 0;
my $j = 0;
my $line;
my @linea;
my $num = 0;
my $output;
my $last_chr = '';


while ($line = <INPUT>) {
	if ($i > 0) {
		@linea = split(/\t/, $line);
		$num = int((($linea[5] / ($max - $min)) * $multiplier)+0.5);
		#print "((($linea[5] / ($max - $min)) * $multiplier)+0.5);\n";
		if($num > 0){
			for ($j=0; $j<$num; $j++){
				print OUTPUT $line;
			}
		}
		if ($linea[0] ne $last_chr) {
			print "Processing Chromosome $linea[0]\n";
		}
		$last_chr = $linea[0];
	}
	$i++;
}
close INPUT;
close OUTPUT;
print "\n$outfile file saved\n\n";
exit;
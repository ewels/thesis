#!/usr/bin/perl
use warnings;
use strict;

my ($infile,$outfile,$window,$min_count) = @ARGV;

unless (defined $min_count) {
  die "Usage is make_windows.pl [input file] [output file] [window size] [minimum count]\n";
}

open (IN,$infile) or die "Can't read $infile: $!";

if (-e $outfile) {
  print "$outfile already exists.  Overwrite [Y/N]?";
  my $answer = <STDIN>;
  if ($answer !~ /^y/i) {
    die "Exiting...\n";
  }
}

open (OUT,'>',$outfile) or die "Can't write to $outfile: $!";

# This is the store of lines we're currently dealing with
my @current_lines;

# This is the chromosome we're dealing with
my $last_chromosome = '';

# Record the header
my $header = <IN>;
chomp $header;

my @header = split(/\t/,$header);

for my $l (0 .. $#header) {
	warn "header[$l] = $header[$l]\n";
}

print OUT join("\t",qw(Chr Start End Count),@header[10..$#header]),"\n";

while (<IN>) {

  chomp;
  my @sections = split(/\t/);

  # Check for putting out a status update
  if ($sections[1] ne $last_chromosome) {
    warn "Processing Chromosome $sections[1]\n";
    $last_chromosome = $sections[1];
  }

  if (@current_lines) {
    # Check if we can add this to the existing set

    if ($current_lines[0]->[1] eq $sections[1]
			and int ($current_lines[0]->[2] / $window) == int ($sections[3] / $window)) {
      push @current_lines,\@sections;
    }
    else {
      check_valid_group(@current_lines);
      @current_lines = (\@sections);
    }

  }
  else {
    push @current_lines, \@sections;
  }

}

sub check_valid_group {

  my @sections = @_;

  if (@sections < $min_count) {
    return;
  }

  # If we're keeping this then we need to work out
  # the position of the window and record the 
  # number of hits and the average measure for
  # each sample

  my $window_start = $window * int($sections[0]->[2]/$window);
  my $window_end = $window_start + $window;

  my @line = ($sections[0]->[1],$window_start,$window_end, scalar @sections);

  for my $index (10..$#header) {
    my $average;
    foreach my $probe (@sections) {
      $average += $probe->[$index];
    }
    $average /= @sections;
    push @line, $average;
  }
  print OUT join("\t",@line),"\n";

}


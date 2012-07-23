#!/usr/bin/perl
use warnings;
use strict;

### This script intends to analyse Phil's FastQ files and sort the sequences according to the first 3 bp barcode information. The experiment was a 4C
### analysis with two different baits which were ligated to prey sequences using an AseI digest. Therefore all of the sequences should look like:
### XXX ATTAAT NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN whereby XXX is the bait-specific barcode information, ATTAAT is the AseI site, and NNNN is the prey
### sequence.

my %baits = (
	     MLL => {
	     	     name => 'MLL_bait',
		     sequence => 'TTT',
		    },
 	     BCR => {
 		     name => 'BCR_bait',
		     sequence => 'GGA',
		    }
	    );
my $barcode_detected = 0;
my $count = 0;
my $asei_counter = 0;
my $nla_counter = 0;
my $double_asei_counter =0;
my $nla_too_short_counter = 0;
my $sequence_in_file = shift @ARGV;
my $results;

my @bait=();
foreach my $bait (keys %baits){
  # resetting all counters
  ($double_asei_counter,$nla_too_short_counter,$count,$asei_counter,$nla_counter,$barcode_detected) = (0,0,0,0,0);
  $results = '';
  @bait = ();
  @bait = split (//,$baits{$bait}->{sequence});
  print "length of barcode sequence in $baits{$bait}->{name} is ",scalar@bait," bp\t";
  print join (" ",@bait),"\n";
  $results = create_filehandles($baits{$bait}->{name});

  read_file ();
  warn "Total number of sequences processed: $count\n\n";
  warn "The correct barcode ",@bait," was detected in $barcode_detected cases\n";
  warn "Total number of prey sequences containing an AseI and an NlaIII site which are getting too short (hence being discarded): $nla_too_short_counter\n";
  warn "Total number of prey sequences containing two  AseI sites (hence being discarded): $double_asei_counter\n";
  warn "Total number of prey sequences containing a AseI site as well as an NlaIII site which are still large enough to map: $nla_counter\n";
  warn "Total number of prey sequences with an intact AseI site printed out (including long enough NlaIII sequences): $asei_counter\n\n";
}

sub create_filehandles {
  my $bait_name = shift;
  my %filehandles;
  my $outfile = "${bait_name}_seqs.txt";
  open ($filehandles{isbait},'>',$outfile) or die $!;
  warn "Writing bait-side output to $outfile\n";
  return \%filehandles;
}

sub read_file {
  open (my $SEQ,$sequence_in_file) or die "Can't read sequence: $!";
  while (1) {
    my $seq = read_next_sequence ($SEQ);
    last unless ($seq);
    #   last if ($count > 1000000);
    ++$count;
    if ($count % 2500000 == 0) {
      warn "Processed $count sequences\n";
    }
    process_sequence($seq);
  }
}
sub read_next_sequence {
  my ($fh) = @_;
  my $seq;
  # First line should be the id
  my $id_line = <$fh>;
  return undef unless ($id_line);
  $seq->{id} = $id_line;
  chomp $seq->{id};

  # Next line is the actual sequence
  $seq->{seq} = <$fh>;
  chomp $seq->{seq};

  # next line is like first line just with a + in the start
  $seq->{third_line} = <$fh>;
  chomp $seq->{third_line};

  # finally we get the quality string
  $seq->{qual} = <$fh>;
  chomp $seq->{qual};

  return $seq;
}


sub process_sequence {
  my $seq = shift;
  # $seq is a sequence which contains all 4 lines of the FastQ file. None of the lines contains a \n.
  if (is_bait($seq)) {
    ++$barcode_detected;
    process_bait_end($seq);
  }
}

sub is_bait {
  my $seq = shift;

  my @seq = split(//,$seq->{seq});
  my @quality = split(//,$seq->{qual});

  my $n_count = 0;
  my $match_count = 0;
  my $mismatch_count = 0;

  for my $index (0..$#bait) {
    last if ($index > $#seq);

    if ($seq[$index] eq 'N') {
      ++$n_count;
      next;
    }

    if ($seq[$index] eq $bait[$index]) {
      ++$match_count;
    }
    else {
      ++$mismatch_count;
    }
  }

  # barcode sequence is a perfect bait match
  if ($match_count == 3) {
    return 1;
  }
  else{
    return 0;
  }
}

sub process_bait_end {
  my $seq = shift;

  # Check for the presence of an AseI site
  if ($seq->{seq} =~ /^(\w{3})(ATTAAT\w+)$/ ){
    my $barcode = $1;
    my $downstream_sequence = $2;
    my $offset = length($seq->{seq})-length($downstream_sequence);

    # kicking the sequence out if it contained another AseI site
    if ($downstream_sequence =~ /^ATTAAT\w*ATTAAT/){
      ++$double_asei_counter;
      return;
    }


    # Recognition site for NlaIII is: CATG (CATG #cut)
    if ($downstream_sequence =~ /^(.*CATG)/ ){
      $downstream_sequence = $1;
      # discarding very short fragments
      if (length($downstream_sequence) < 25){
	#	print "$downstream_sequence too short!\n";
	++$nla_too_short_counter;
	return;
      }

      # otherwise we count how many fragments do contain both AseI and NlaIII
      ++$nla_counter;
      ++$asei_counter;
    }

    # does contain an intact AseI site but no NlaIII site
    elsif ($downstream_sequence !~ /CATG/ ){
      ++$asei_counter;
    }
    else {
      die "Either there is an NlaIII site or not!\n";
    }
    # print "$barcode\t$downstream_sequence\t$seq->{seq}\n";
    # Print out the sequence downstream of the to a FastQ file
    print {$results->{isbait}} $seq->{id},"\n";
    print {$results->{isbait}} $downstream_sequence,"\n";
    print {$results->{isbait}} $seq->{third_line},"\n";
    print {$results->{isbait}} substr($seq->{qual},$offset,length $downstream_sequence),"\n";
  }
}

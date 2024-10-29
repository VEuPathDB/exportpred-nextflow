#!/usr/bin/perl

use strict;

use Getopt::Long;

use Bio::SeqFeature::Generic;
use Bio::Tools::GFF;


my ($file, $outputFile);

&GetOptions(
    "inputFile=s"       => \$file,
    "outputFile=s"       => \$outputFile,
    );


open(IN, "<$file") or die "Couldn't open file '$file': $!\n";

open(OUT, ">$outputFile") or die "Couldn't open file $outputFile for writing: $!";
my $source = "veupathdb";

my $count = 0;
while (<IN>) {
    chomp;
    my ($sourceId, $type, $score, $parse) =
      m/^ (\S+)            # sourceId
          .*? \s+
          (KLD | RLE)      # type
          \s+
          (\d+\.\d+)  # score
          \s+
          (\S+)            # parse
        $/x;


    my $offset = 0;


    my @subfeats;

    # example: [a-met:M][a-leader:YSN][a-hydrophobic:LSLC]

    my $uid = 1;
    while ($parse =~ m/\[
                       ( [^ \: ]+ )  # name
                       \:
                       ( [^ \] ]+ )  # seq
                       \]
                      /xg) {
      my ($name, $seq) = ($1, $2);



      my $seqLength = length($seq);

      my $subfeat = Bio::SeqFeature::Generic->new(
          -start        => $offset + 1,
          -end          => $offset + $seqLength,
          -primary      => $name, # -primary_tag is a synonym
          -seq_id       => $sourceId,
          -source_tag   => $source,
          -score        => ".",
          -tag          => { Parent => "${sourceId}_${type}",
                            ID => "${sourceId}_${type}_${uid}",
                            } );


      $subfeat->gff_format(Bio::Tools::GFF->new(-gff_version => 3));

      push @subfeats, $subfeat;

      $offset += $seqLength;
      $uid++;
    }

      my $feat = Bio::SeqFeature::Generic->new(
          -start        => 1,
          -end          => $offset,
          -primary      => $type, # -primary_tag is a synonym
          -seq_id       => $sourceId,
          -source_tag   => $source,
          -score        => $score,
          -frame        => ".",
          -tag          => {ID => "${sourceId}_${type}",
                            } );


    $feat->gff_format(Bio::Tools::GFF->new(-gff_version => 3));
    print OUT $feat->gff_string, "\n" ;
    foreach(@subfeats) {
        print OUT $_->gff_string, "\n" ;
   }
}


close IN;
close OUT;

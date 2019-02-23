#!/usr/bin/perl

use strict;
use Bio::SeqIO;
use Data::Dumper;

if (scalar @ARGV==0){
  die "USAGE: $0 <samtools calmd output> <bowtie2 output> <reference sequence> <output file name> <QualitySymbols> <Quality_cutoff>\n";
}


my $seq_in = Bio::SeqIO -> new(-file => $ARGV[2], -format => 'fasta');
my $seq_obj = $seq_in -> next_seq();

my $seq = $seq_obj-> seq();

my $wt= substr($seq,50,240);   #### protein starts after 50 bp in reference


my @qual=`cat $ARGV[4]`;
my %badqual;
for (my $count=0; $count<$ARGV[5]; $count++){
  chomp $qual[$count];
  $badqual{$qual[$count]}++;
}


my %code = (
'TCA' => 'S', # Serine
'TCC' => 'S', # Serine
'TCG' => 'S', # Serine
'TCT' => 'S', # Serine
'TTC' => 'F', # Phenylalanine
'TTT' => 'F', # Phenylalanine
'TTA' => 'L', # Leucine
'TTG' => 'L', # Leucine
'TAC' => 'Y', # Tyrosine
'TAT' => 'Y', # Tyrosine
'TAA' => '-', # Stop
'TAG' => '-', # Stop
'TGC' => 'C', # Cysteine
'TGT' => 'C', # Cysteine
'TGA' => '-', # Stop
'TGG' => 'W', # Tryptophan
'CTA' => 'L', # Leucine
'CTC' => 'L', # Leucine
'CTG' => 'L', # Leucine
'CTT' => 'L', # Leucine
'CCA' => 'P', # Proline
'CAT' => 'H', # Histidine
'CAA' => 'Q', # Glutamine
'CAG' => 'Q', # Glutamine
'CGA' => 'R', # Arginine
'CGC' => 'R', # Arginine
'CGG' => 'R', # Arginine
'CGT' => 'R', # Arginine
'ATA' => 'I', # Isoleucine
'ATC' => 'I', # Isoleucine
'ATT' => 'I', # Isoleucine
'ATG' => 'M', # Methionine
'ACA' => 'T', # Threonine
'ACC' => 'T', # Threonine
'ACG' => 'T', # Threonine
'ACT' => 'T', # Threonine
'AAC' => 'N', # Asparagine
'AAT' => 'N', # Asparagine
'AAA' => 'K', # Lysine
'AAG' => 'K', # Lysine
'AGC' => 'S', # Serine
'AGT' => 'S', # Serine
'AGA' => 'R', # Arginine
'AGG' => 'R', # Arginine
'CCC' => 'P', # Proline
'CCG' => 'P', # Proline
'CCT' => 'P', # Proline
'CAC' => 'H', # Histidine
'GTA' => 'V', # Valine
'GTC' => 'V', # Valine
'GTG' => 'V', # Valine
'GTT' => 'V', # Valine
'GCA' => 'A', # Alanine
'GCC' => 'A', # Alanine
'GCG' => 'A', # Alanine
'GCT' => 'A', # Alanine
'GAC' => 'D', # Aspartic Acid
'GAT' => 'D', # Aspartic Acid
'GAA' => 'E', # Glutamic Acid
'GAG' => 'E', # Glutamic Acid
'GGA' => 'G', # Glycine
'GGC' => 'G', # Glycine
'GGG' => 'G', # Glycine
'GGT' => 'G'  # Glycine
);

my %counts;
my @wtcodons= unpack("(A3)*", $wt);
open FILE1, $ARGV[0];
open FILE2, $ARGV[1];

my $bad1=0;       #bad alignment
my $bad2=0;       #more than one codon mutated
my $bad3=0;       #low quality wt bases
my $bad4=0;       #low quality mutant bases
my $reads=0;
while (!eof(FILE1) and !eof(FILE2)){
  my $line1=<FILE1>;
  my $line2=<FILE2>;
  
  next if ($line1 =~ /^@/);
  my @parts1 = split /\t/, $line1;
  my $pos = $parts1[3];
  my $n=51-$pos;        ####### location of stop codon
  my $off=$parts1[5];
  my $start=0;
  my $match=0;
  if ($off =~ /(\d+)S(\d+)M/) {
      $start=$1;
      $match=$2;
  }

  
  if ($match-(240+$n)<0){      
    $bad1++;
    next;
  };
  
  my $str = substr($parts1[9],$n+$start,240);
  
  #### QUALITY
  my $qstr = substr($parts1[10],$n+$start,240);
  my @qstr = split //, $qstr; 
  
  my @parts2=split /\t/, $line2;
  my $str2 = substr($parts2[9],$n+$start,240);
  
  my @codons2= unpack("(A3)*", $str2);
  unless ($codons2[$#codons2] eq 'GGC' or $codons2[$#codons2-1] eq 'GAT'){   #### These codons should be at residues 80 and 79
    $bad1++;next;
  }
  
  if ($str !~ /[ACTG]/){
      $counts{'wt'}++;
  } else {
    my @codons= unpack("(A3)*", $str);
    my $count;
    my $changed;
    for (my $codon=0; $codon<@codons; $codon++){
      if ($codons[$codon] =~ /[ACTG]/){
        $count++;
        $changed=$codon;
      }
    }
    if ($count ==1){
      my @nts=split//,$codons2[$changed];
      my @wtnts=split//,$wtcodons[$changed];
      
      my $test=0;
      my $base1=($changed-1)*3;
      for (my $i=$base1; $i<$base1+3; $i++){
        my $B=$i-$base1;
        if (exists $badqual{$qstr[$i]} and $nts[$B] ne $wtnts[$B]){
          $test++; last;
        }
      }
      
      if ($test>0){
        $bad4++;
        next;
      }
      
      
      my $codon=join '', @nts;
      my $revcom=$codon;
      $revcom =~ tr/ACTG/TGAC/;
      $revcom = reverse $revcom;
      my $aa = $code{$revcom};
      my $residue=160-$changed;      ########## Count residues from stop codon
      my $name = "$residue"."_$aa"."_$revcom";
      $counts{$name}++;
    } else {
      $bad2++;
      next;
    }
  }
  $reads++;
  if ($reads %100 ==0){
    print STDERR "Processed reads: $reads\r";
  }
}

print STDERR "\nTOTAL:$reads\nBAD:$bad1\t$bad2\t$bad3\t$bad4\n";

open OUT, ">$ARGV[3]";
foreach my $seq (keys %counts){
  print OUT "$seq\t$counts{$seq}\n";
}




#!/usr/bin/perl
#
# This script parses PacBio sequences and identifies the likely mod_mdg4 variant contained
# It is specific to the experimental design of this article (primers used).
#
# Example usage:
# S=Body_1 ; gzip -dc $S.fa.gz | ./translate.pl $S
#
# (c) N.Guex
# Bioinformatics Competence Center (BICC)
# University of Lausanne
# Switzerland
# https://bix.unil.ch
#

use strict;
my $SAMPLE = $ARGV[0];

my(%codons) = (
'UUU' => 'F',
'UUC' => 'F',
'UUA' => 'L',
'UUG' => 'L',
'CUU' => 'L',
'CUC' => 'L',
'CUA' => 'L',
'CUG' => 'L',
'AUU' => 'I',
'AUC' => 'I',
'AUA' => 'I',
'AUG' => 'M',
'GUU' => 'V',
'GUC' => 'V',
'GUA' => 'V',
'GUG' => 'V',
'UCU' => 'S',
'UCC' => 'S',
'UCA' => 'S',
'UCG' => 'S',
'CCU' => 'P',
'CCC' => 'P',
'CCA' => 'P',
'CCG' => 'P',
'ACU' => 'T',
'ACC' => 'T',
'ACA' => 'T',
'ACG' => 'T',
'GCU' => 'A',
'GCC' => 'A',
'GCA' => 'A',
'GCG' => 'A',
'UAU' => 'Y',
'UAC' => 'Y',
'UAA' => '.',
'UAG' => '.',
'CAU' => 'H',
'CAC' => 'H',
'CAA' => 'Q',
'CAG' => 'Q',
'AAU' => 'N',
'AAC' => 'N',
'AAA' => 'K',
'AAG' => 'K',
'GAU' => 'D',
'GAC' => 'D',
'GAA' => 'E',
'GAG' => 'E',
'UGU' => 'C',
'UGC' => 'C',
'UGA' => '.',
'UGG' => 'W',
'CGU' => 'R',
'CGC' => 'R',
'CGA' => 'R',
'CGG' => 'R',
'AGU' => 'S',
'AGC' => 'S',
'AGA' => 'R',
'AGG' => 'R',
'GGU' => 'G',
'GGC' => 'G',
'GGA' => 'G',
'GGG' => 'G',
);

# read known transcripts database.
open F, "mod_mdg4.infos.txt" or die "oops\n";
my @CTERM;
my @INFOS;
my $annotCnt = 0;
while(<F>)
{
    chomp;
    $INFOS[$annotCnt]   = substr($_,1);
    $CTERM[$annotCnt]   = $INFOS[$annotCnt];
    $CTERM[$annotCnt++] =~ s/.*SFVDTSGDQGNTEAQ//;  # keep only C-terminal isoform differences
}
close F;
# ---------------------------------------------
# annotate sequences
while(<STDIN>)
{
    # read sequence
    chomp;
    my $hdr = substr($_,1);
    my $dna = <STDIN>;
    chomp($dna);

    # identify first start codon (look only at beginning since we use primers and we know the expected start)
    my $ATGpos = index($dna,'ATG');
    if ($ATGpos == -1 || $ATGpos > 56) { print $hdr . "\t$SAMPLE\tNA\tNA\tNA\tNA\n";  next; }

    # get CDS codons
    my $rna = substr($dna,$ATGpos); # trim 5p UTR
    $rna =~ s/T/U/g;
    my @codon;
    my $j = 0;
    for (my $i = 0 ; $i < length($rna); $i = $i+3)
    {
        @codon[$j] = substr( $rna, $i, 3);
        if    (length(@codon[$j]) == 1) { @codon[$j] .= 'AA'; } # assume we have polyA at the end, which was trimmed.
        elsif (length(@codon[$j]) == 2) { @codon[$j] .= 'A';  } # assume we have polyA at the end, which was trimmed.
        $j++;
    }

    # translate CDS
    my $cds = "";
    my $prot = "";
    for (my $j = 0 ; $j < length($rna)/3; $j++)
    {
        if( exists($codons{@codon[$j]}) )
        {
            my $aa = $codons{@codon[$j]};
            if ($aa eq '.') { last; }
            $cds .= $codon[$j];
            $prot .= $aa;
        }
        else
        {
            print "ERROR at $j:@codon[$j]:$rna\n";
            exit(1);
        }
    }
    
    # annotate CDS with protein
    my $hit = "";
    my $annotHits = 0;
    for (my $i = 0; $i < $annotCnt; $i++)
    {
        my ($h,$r,$p) = split("\t",$INFOS[$i]);
        $h =~ s/.gb$//;

        if ($prot eq $p) # exact CDS match ?
        {
            $hit .= "$h;";
            $annotHits++;
        }
        else
        {
            if ($prot =~ /^$p/) # exact CDS match but with C terminal extension ?
            {
                $hit .= "${h}_CtermExt;";
                $annotHits++;
            }
            elsif ($p =~ /^$prot/)  # exact CDS match but trunctated (possibly due to sequencing error)
            {
                $hit .= "${h}_Truncated;";
                $annotHits++;
            }
            elsif ($prot =~ /$CTERM[$i]$/) # exact CDS match of the C-terminal end of a given variant (probably sequencing error upstream)
            {
                $hit .= "${h}_Variant;";
                $annotHits++;
            }
        }
    }
    
    # print result
    if ($hit eq "") { $hit = "mod_mdg4_unknown;"; }
    print $hdr . "\t" . $SAMPLE . "\t" . $annotHits . "\t" . $hit . "\t" . $prot . "\t" . $cds . "\n";
}
# ---------------------------------------------

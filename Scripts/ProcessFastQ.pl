#!/usr/bin/perl
#
# This script parses PacBio sequences and identifies the likely mod_mdg4 variant contained
# It is specific to the experimental design of this article (primers used).
#
# Example usage:
# gzip -dc fastq/m84177_250731_113831_s3.hifi_reads.bc1005_5p--IsoSeqX_3p.hifi_reads.fastq.gz | ./ProcessFastQ.pl Body_1 | gzip>  Body_1.fa.gz &
# gzip -dc fastq/m84177_250731_113831_s3.hifi_reads.bc1007_5p--IsoSeqX_3p.hifi_reads.fastq.gz | ./ProcessFastQ.pl Body_2 | gzip>  Body_2.fa.gz &
# gzip -dc fastq/m84177_250731_113831_s3.hifi_reads.bc1008_5p--IsoSeqX_3p.hifi_reads.fastq.gz | ./ProcessFastQ.pl Head_1 | gzip>  Head_1.fa.gz &
# gzip -dc fastq/m84177_250731_113831_s3.hifi_reads.bc1012_5p--IsoSeqX_3p.hifi_reads.fastq.gz | ./ProcessFastQ.pl Head_2 | gzip>  Head_2.fa.gz &
# wait
#
# (c) N.Guex
# Bioinformatics Competence Center (BICC)
# University of Lausanne
# Switzerland
# https://bix.unil.ch
#

use strict;
my $sample = $ARGV[0];
while(<STDIN>)
{
    my $HDR = $_;
    my $SEQ = <STDIN>;
    <STDIN>; # skip separator
    <STDIN>; # skip QC
    chomp($SEQ);
    my $len = length($SEQ);
    my $orientation = 'FwdKO';
    if ($SEQ =~ /^AGTCAAGAGCCAACAAACGCATAG/) { $orientation = 'FwdOK'; }
    else
    {
        my $hit = index($SEQ,'CTATGCGTTTGTTGGCTCTTGACT');
        if ($hit != -1)
        {
            if (($len - $hit) == 24) { $orientation = 'RevOK'; }
        }
    }
    
    if ($orientation eq 'FwdKO')
    {
        if ($SEQ =~ /^TTTTTTTTTTTTTTTTTT/) { $orientation = 'RevKO'; }
    }

    if ($orientation eq 'RevOK' || $orientation eq 'RevKO')
    {
        $SEQ = reverse($SEQ);
        $SEQ =~ tr/ACGTNacgtn/TGCANtgcan/;
    }
    print ">$sample\_$orientation:$len:$HDR$SEQ\n";
}

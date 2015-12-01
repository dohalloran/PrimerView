#!/usr/bin/perl

use strict;
use warnings;
use SequenceIO;
use PRIMERVIEW;
use Getopt::Std;

#get options from command line
my %opts;
getopt( 'abcdefghijk', \%opts );
my $fasta           = $opts{a};
my $five_prime_end  = $opts{b} || "150";
my $three_prime_end = $opts{c} || "150";
my $kmer_max        = $opts{d} || "28";
my $kmer_min        = $opts{e} || "22";
my $clamp           = $opts{f} || "Y";
my $higher_gc       = $opts{g} || "60";
my $lower_gc        = $opts{h} || "40";
my $upper_tm        = $opts{i} || "68";
my $lower_tm        = $opts{j} || "55";
my $spec            = $opts{k} || "N";

#parse the input file into arrays of sequence and id using SequenceIO.pm
my ( $seq_file, $id_file ) = parser($fasta);
my $specificty = join( ",", @$seq_file );

#input each sequence into PRIMERVIEW.pm 
foreach my $sequence (@$seq_file) {
    my $len_seq  = length $sequence;
    my $shift_id = shift(@$id_file);
    $shift_id =~ tr/a-zA-Z0-9//cd;
    my $id_uniq  = $len_seq . "'" . $shift_id;
    my $len_uniq = $shift_id . "'" . $len_seq;
    primerview(
        $fasta,           $sequence, $len_seq,    $id_uniq,
        $len_uniq,        $shift_id, $specificty, $five_prime_end,
        $three_prime_end, $kmer_max, $kmer_min,   $clamp,
        $higher_gc,       $lower_gc, $upper_tm,   $lower_tm,
        $spec
    );
    align_muscle();
    align_convert();
    graphics();
    graphics_all_primers();
    clean_up();
}


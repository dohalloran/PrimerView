package SequenceIO;

use strict;
use warnings;
use Bio::SeqIO;
use base 'Exporter';

our @EXPORT = qw/ parser /;

sub parser {
    my $fasta = shift;
    my @sequence;
    my @id;

    #create an instance of a Seqio object
    my $seqio = Bio::SeqIO->new(
        -file   => $fasta,
        -format => "fasta",
    );

    while ( my $seqobj = $seqio->next_seq() ) {
        my $seqStream = $seqobj->seq();
        my $idStream  = $seqobj->id();
        push @sequence, $seqStream;
        push @id,       $idStream;
    }
    #send back a reference to the arrays containing seqs and ids
    return ( \@sequence, \@id );

}

1;

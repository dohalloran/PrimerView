package PRIMERVIEW;

#############################################################
#
# PrimerView version 1.0
# Designs forward and reverse primers
# from fastA formatted files of any size
# and generates graphical outputs of primer
# location and distribution
# Written by Damien O'Halloran
# March 2015
# readme.txt file for usage
#
#############################################################

use strict;
use warnings;
use Getopt::Std;
use Cwd;
use Bio::SeqIO;
use Bio::Tools::Run::Alignment::Muscle;
use Bio::AlignIO;
use Bio::Align::Graphics;
use Bio::Graphics;
use Bio::SeqFeature::Generic;
use base 'Exporter';

#########################

our @EXPORT =
  qw(primerview align_muscle align_convert graphics graphics_all_primers clean_up);

#########################

our $VERSION = '1.0';

#########################
sub new {
    my $class = shift;
    return bless {}, $class;
}
#########################

#get options from command line
my %opts;
getopt( 'abcdefghij', \%opts );
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

=pod
 
=head1 INPUT PARAMETERS
 
This script can have 10 parameters. The input fasta file name, the 
5' and 3' search areas, the primer max and min, the GC% max and min,
 the Tm max and min, and whether to include a 3' GC clamp. Default 
 options are in place for all parameters except the input file. To 
 change these defaults, simply provide a parameter in the command-line
 e.g. -c 28 for primer upper length.
 All fasta format description lines must be unique for each gene

=cut

#define strings and hashes
my @array_length;
my @array_name;
my $Tm2;
my $Tm;
my $kmer;
my $id;
my $out_image;
my $kmer_diff  = $kmer_max - $kmer_min;
my $outputfile = $fasta . "_primers.txt";   #good to add PATH to appropriate dir
my $selfie_cuttof = 12;                     #change to make less/more stringent
my %nn_s          = (
    "AA" => 240,
    "AC" => 173,
    "AG" => 208,
    "AT" => 239,
    "AN" => 215,
    "CA" => 129,
    "CC" => 266,
    "CG" => 278,
    "CT" => 208,
    "CN" => 220,
    "GA" => 135,
    "GC" => 267,
    "GG" => 266,
    "GT" => 173,
    "GN" => 210,
    "TA" => 169,
    "TC" => 135,
    "TG" => 129,
    "TT" => 240,
    "TN" => 168,
    "NA" => 168,
    "NC" => 210,
    "NG" => 220,
    "NT" => 215,
    "NN" => 203,
    "aa" => 240,
    "ac" => 173,
    "ag" => 208,
    "at" => 239,
    "an" => 215,
    "ca" => 129,
    "cc" => 266,
    "cg" => 278,
    "ct" => 208,
    "cn" => 220,
    "ga" => 135,
    "gc" => 267,
    "gg" => 266,
    "gt" => 173,
    "gn" => 210,
    "ta" => 169,
    "tc" => 135,
    "tg" => 129,
    "tt" => 240,
    "tn" => 168,
    "na" => 168,
    "nc" => 210,
    "ng" => 220,
    "nt" => 215,
    "nn" => 203
);
my %nn_h = (
    "AA" => 91,
    "AC" => 65,
    "AG" => 78,
    "AT" => 86,
    "AN" => 80,
    "CA" => 58,
    "CC" => 110,
    "CG" => 119,
    "CT" => 78,
    "CN" => 91,
    "GA" => 56,
    "GC" => 111,
    "GG" => 110,
    "GT" => 65,
    "GN" => 85,
    "TA" => 60,
    "TC" => 56,
    "TG" => 58,
    "TT" => 91,
    "TN" => 66,
    "NA" => 66,
    "NC" => 85,
    "NG" => 91,
    "NT" => 80,
    "NN" => 80,
    "aa" => 91,
    "ac" => 65,
    "ag" => 78,
    "at" => 86,
    "an" => 80,
    "ca" => 58,
    "cc" => 110,
    "cg" => 119,
    "ct" => 78,
    "cn" => 91,
    "ga" => 56,
    "gc" => 111,
    "gg" => 110,
    "gt" => 65,
    "gn" => 85,
    "ta" => 60,
    "tc" => 56,
    "tg" => 58,
    "tt" => 91,
    "tn" => 66,
    "na" => 66,
    "nc" => 85,
    "ng" => 91,
    "nt" => 80,
    "nn" => 80
);

#########################
sub primerview {

    #create an instance of a Seqio object
    my $seqio = Bio::SeqIO->new(
        -file   => $fasta,
        -format => "fasta",
    );

    while ( my $seqobj = $seqio->next_seq() ) {
        my $sequence = $seqobj->seq();
        $id = $seqobj->id();
        $id =~ tr/a-zA-Z0-9//cd;
        my $len_seq  = length $sequence;
        my $id_uniq  = $len_seq . "'" . $id;
        my $len_uniq = $id . "'" . $len_seq;

        #declare output file
        $out_image = "GRAPHIC_$id.txt";

    #skip to next if selected 5' or 3' length is longer than the sequence length
        if ( $five_prime_end > length $sequence ) {
            next;
        }
        if ( $three_prime_end > length $sequence ) {
            next;
        }

        #skip to next if sequence length is <100bps
        if ( length $sequence < 100 ) {
            next;
        }

######## FORWARD PRIMER
        #start counting
        my $start = 1;
        for ( my $i = $start - 1 ; $i < $five_prime_end - 1 ; $i += 2 ) {
            $kmer = int( rand($kmer_diff) ) + $kmer_min;
            $_ = substr( $sequence, $i, $kmer );

            #get self complementarity score
            my $revF = reverse($_);
            $revF =~ tr/ATGCatgc/TACGtacg/;
            my $selfie_score = selfie( $_, $revF );

            #Count Gs and Cs
            my $countGC = tr/GCgc//;

            #Calculate percent GC
            my $percentGC = 100 * $countGC / $kmer;
            my $percentGCrounded = sprintf( "%0.1f", $percentGC );

            #calculate Tm
            if ( $kmer <= 36 ) {
                $Tm = calcTm( $_, 100, 50 );
            }
            else {
                $Tm = calclongTm( $_, 100, 50 );
            }
            my $Tmrounded = sprintf( "%0.1f", $Tm );

            my $primer_end = $i + $kmer;

            #define dinucleotide repeats and repetitive sequence
            #and print results if statements are matched
            if (   open( FUSIONFILE, ">>$outputfile" )
                && $Tmrounded ge $lower_tm
                && $Tmrounded le $upper_tm
                && $percentGC ge $lower_gc
                && $percentGC le $higher_gc
                && $selfie_score < $selfie_cuttof
                && $_ !~ /AAAAA/i
                && $_ !~ /TTTTT/i
                && $_ !~ /GGGGG/i
                && $_ !~ /CCCCC/i
                && $_ !~ /ATATATAT/i
                && $_ !~ /TATATATA/i
                && $_ !~ /GCGCGCGC/i
                && $_ !~ /CGCGCGCG/i
                && $_ =~ m/gc$/i
                && $clamp eq "Y" )
            {
                print
"$id\t$i\t$Tmrounded degC\tF:$_\t$selfie_score\t$percentGCrounded%\n";
                push @array_length, $len_uniq;
                push @array_name,   $id_uniq;
                print FUSIONFILE
"$id\t$i\t$Tmrounded degC\tF:$_\t$selfie_score\t$percentGCrounded%\n";
                open( OLIGOS, ">>$out_image" ) or die;
                print OLIGOS "$i\t$selfie_score\t$i\t$primer_end\n";
                close(OLIGOS);

                #declare output files
                my $out = "$i.fa";
                open( ALIGN, ">$out" ) or die;
                print ALIGN ">$id\n$sequence\n>$i\n$_\n";
                close(ALIGN);
                close(FUSIONFILE);
            }
            elsif (open( FUSIONFILE, ">>$outputfile" )
                && $Tmrounded ge $lower_tm
                && $Tmrounded le $upper_tm
                && $percentGC ge $lower_gc
                && $percentGC le $higher_gc
                && $selfie_score < $selfie_cuttof
                && $_ !~ /AAAAA/i
                && $_ !~ /TTTTT/i
                && $_ !~ /GGGGG/i
                && $_ !~ /CCCCC/i
                && $_ !~ /ATATATAT/i
                && $_ !~ /TATATATA/i
                && $_ !~ /GCGCGCGC/i
                && $_ !~ /CGCGCGCG/i
                && $_ =~ m/cg$/i
                && $clamp eq "Y" )
            {
                print
"$id\t$i\t$Tmrounded degC\tF:$_\t$selfie_score\t$percentGCrounded%\n";
                push @array_length, $len_uniq;
                push @array_name,   $id_uniq;
                open( OLIGOS, ">>$out_image" ) or die;
                print OLIGOS "$i\t$selfie_score\t$i\t$primer_end\n";
                close(OLIGOS);

                #declare output files
                my $out = "$i.fa";
                open( ALIGN, ">$out" ) or die;
                print ALIGN ">$id\n$sequence\n>$i\n$_\n";
                close(ALIGN);
                print FUSIONFILE
"$id\t$i\t$Tmrounded degC\tF:$_\t$selfie_score\t$percentGCrounded%\n";
                close(FUSIONFILE);
            }
            elsif (open( FUSIONFILE, ">>$outputfile" )
                && $Tmrounded ge $lower_tm
                && $Tmrounded le $upper_tm
                && $percentGC ge $lower_gc
                && $percentGC le $higher_gc
                && $selfie_score < $selfie_cuttof
                && $_ !~ /AAAAA/i
                && $_ !~ /TTTTT/i
                && $_ !~ /GGGGG/i
                && $_ !~ /CCCCC/i
                && $_ !~ /ATATATAT/i
                && $_ !~ /TATATATA/i
                && $_ !~ /GCGCGCGC/i
                && $_ !~ /CGCGCGCG/i
                && $clamp eq "N" )
            {
                print
"$id\t$i\t$Tmrounded degC\tF:$_\t$selfie_score\t$percentGCrounded%\n";
                push @array_length, $len_uniq;
                push @array_name,   $id_uniq;
                open( OLIGOS, ">>$out_image" ) or die;
                print OLIGOS "$i\t$selfie_score\t$i\t$primer_end\n";
                close(OLIGOS);

                #declare output files
                my $out = "$i.fa";
                open( ALIGN, ">$out" ) or die;
                print ALIGN ">$id\n$sequence\n>$i\n$_\n";
                close(ALIGN);
                print FUSIONFILE
"$id\t$i\t$Tmrounded degC\tF:$_\t$selfie_score\t$percentGCrounded%\n";
                close(FUSIONFILE);
            }
        }

######## REVERSE PRIMER

        #start counting for reverse primer
        for (
            my $j = length($sequence) - $three_prime_end ;
            $j < length($sequence) ;
            $j += 2
          )
        {
            $kmer = int( rand($kmer_diff) ) + $kmer_min;
            $_ = substr( $sequence, $j, $kmer );

            #rev comp
            my $revR = reverse($_);
            $revR =~ tr/ATGCatgc/TACGtacg/;

            #get self complementarity score
            my $selfie_scoreR = selfie( $_, $revR );

            #Count Gs and Cs
            my $count_GC = tr/GCgc//;

            #Calculate percent GC
            my $percent_GC = 100 * $count_GC / $kmer;
            my $percentGC_rounded = sprintf( "%0.1f", $percent_GC );

            #calculate Tm
            if ( $kmer <= 36 ) {
                $Tm2 = calcTm( $_, 100, 50 );
            }
            else {
                $Tm2 = calclongTm( $_, 100, 50 );
            }
            my $Tm_rounded = sprintf( "%0.1f", $Tm2 );
            my $primer_start_R = $j + $kmer;

            #define dinucleotide repeats and repetitive sequence
            #and print results if statements are matched
            if (   open( FUSIONFILE, ">>$outputfile" )
                && $Tm_rounded ge $lower_tm
                && $Tm_rounded le $upper_tm
                && $percent_GC ge $lower_gc
                && $percent_GC le $higher_gc
                && $selfie_scoreR < $selfie_cuttof
                && $revR !~ /AAAAA/i
                && $revR !~ /TTTTT/i
                && $revR !~ /GGGGG/i
                && $revR !~ /CCCCC/i
                && $revR !~ /ATATATAT/i
                && $revR !~ /TATATATA/i
                && $revR !~ /GCGCGCGC/i
                && $revR !~ /CGCGCGCG/i
                && $clamp eq "N" )
            {
                open( OLIGOS, ">>$out_image" ) or die;
                print OLIGOS "$j\t$selfie_scoreR\t$primer_start_R\t$j\n";
                close(OLIGOS);
                push @array_length, $len_uniq;
                push @array_name,   $id_uniq;

                #declare output files
                my $out_r = "$j.fa";
                open( ALIGN_R, ">$out_r" ) or die;
                print ALIGN_R ">$id\n$sequence\n>$j\n$_\n";
                close(ALIGN_R);
                print
"$id\t$j\t$Tm_rounded degC\tR:$revR\t$selfie_scoreR\t$percentGC_rounded%\n";
                print FUSIONFILE
"$id\t$j\t$Tm_rounded degC\tR:$revR\t$selfie_scoreR\t$percentGC_rounded%\n";
                close(FUSIONFILE);
            }
            elsif (open( FUSIONFILE, ">>$outputfile" )
                && $Tm_rounded ge $lower_tm
                && $Tm_rounded le $upper_tm
                && $percent_GC ge $lower_gc
                && $percent_GC le $higher_gc
                && $selfie_scoreR < $selfie_cuttof
                && $revR !~ /AAAAA/i
                && $revR !~ /TTTTT/i
                && $revR !~ /GGGGG/i
                && $revR !~ /CCCCC/i
                && $revR !~ /ATATATAT/i
                && $revR !~ /TATATATA/i
                && $revR !~ /GCGCGCGC/i
                && $revR !~ /CGCGCGCG/i
                && $_ =~ m/^gc/i
                && $clamp eq "Y" )
            {
                open( OLIGOS, ">>$out_image" ) or die;
                print OLIGOS "$j\t$selfie_scoreR\t$primer_start_R\t$j\n";
                close(OLIGOS);
                push @array_length, $len_uniq;
                push @array_name,   $id_uniq;

                #declare output files
                my $out_r = "$j.fa";
                open( ALIGN_R, ">$out_r" ) or die;
                print ALIGN_R ">$id\n$sequence\n>$j\n$_\n";
                close(ALIGN_R);
                print
"$id\t$j\t$Tm_rounded degC\tR:$revR\t$selfie_scoreR\t$percentGC_rounded%\n";
                print FUSIONFILE
"$id\t$j\t$Tm_rounded degC\tR:$revR\t$selfie_scoreR\t$percentGC_rounded%\n";
                close(FUSIONFILE);
            }
            elsif (open( FUSIONFILE, ">>$outputfile" )
                && $Tm_rounded ge $lower_tm
                && $Tm_rounded le $upper_tm
                && $percent_GC ge $lower_gc
                && $percent_GC le $higher_gc
                && $selfie_scoreR < $selfie_cuttof
                && $revR !~ /AAAAA/i
                && $revR !~ /TTTTT/i
                && $revR !~ /GGGGG/i
                && $revR !~ /CCCCC/i
                && $revR !~ /ATATATAT/i
                && $revR !~ /TATATATA/i
                && $revR !~ /GCGCGCGC/i
                && $revR !~ /CGCGCGCG/i
                && $_ =~ m/^cg/i
                && $clamp eq "Y" )
            {
                open( OLIGOS, ">>$out_image" ) or die;
                print OLIGOS "$j\t$selfie_scoreR\t$primer_start_R\t$j\n";
                close(OLIGOS);
                push @array_length, $len_uniq;
                push @array_name,   $id_uniq;

                #declare output files
                my $out_r = "$j.fa";
                open( ALIGN_R, ">$out_r" ) or die;
                print ALIGN_R ">$id\n$sequence\n>$j\n$_\n";
                close(ALIGN_R);
                print
"$id\t$j\t$Tm_rounded degC\tR:$revR\t$selfie_scoreR\t$percentGC_rounded%\n";
                print FUSIONFILE
"$id\t$j\t$Tm_rounded degC\tR:$revR\t$selfie_scoreR\t$percentGC_rounded%\n";
                close(FUSIONFILE);
            }
        }
    }

}

####################################
sub calclongTm {
    my $sequence = shift;
    my $DNA_nM   = shift;
    my $K_mM     = shift;
    return 81.5 + ( 16.6 * ( log( $K_mM / 1000.0 ) / log(10) ) ) +
      ( 41.0 * percentGC($sequence) ) - ( 600.0 / length($sequence) );
}

####################################
sub calcTm {
    my $sequence = uc(shift);
    my $DNA_nM   = shift;
    my $K_mM     = shift;
    my $dH       = 0;
    my $dS       = 108;
    my $seq_len  = length($sequence);
    my $i;

    # Compute dH and dS
    for ( $i = 0 ; $i < ( $seq_len - 1 ) ; $i++ ) {
        my $pair = substr( $_, $i, 2 );
        $dH += $nn_h{$pair};
        $dS += $nn_s{$pair};
    }

    $dH *= -100.0;
    $dS *= -0.1;

    return $dH / ( $dS + 1.987 * log( $DNA_nM / 4000000000.0 ) ) - 273.15 +
      16.6 * ( log( $K_mM / 1000.0 ) / log(10) );
}

####################################
sub selfie {
    my $primer_f  = shift;
    my $primer_rc = shift;               #revcomp primer
    my $FLEN      = length $primer_rc;
    my $RC_LEN    = length $primer_f;
    my $D         = [];
    for ( my $t = 0 ; $t <= $FLEN ; ++$t ) {
        $D->[$t][0] = 0;
    }
    for ( my $p = 0 ; $p <= $RC_LEN ; ++$p ) {
        $D->[0][$p] = $p;
    }
    for ( my $t = 1 ; $t <= $FLEN ; ++$t ) {
        for ( my $p = 1 ; $p <= $RC_LEN ; ++$p ) {
            $D->[$t][$p] =

              min3(
                substr( $primer_rc, $t - 1, 1 ) eq
                  substr( $primer_f, $p - 1, 1 )
                ? $D->[ $t - 1 ][ $p - 1 ]
                : $D->[ $t - 1 ][ $p - 1 ] + 1,

                $D->[ $t - 1 ][$p] + 1,

                $D->[$t][ $p - 1 ] + 1
              );
        }
    }

    my @matches   = ();
    my $bestscore = 10000000;

    for ( my $t = 1 ; $t <= $FLEN ; ++$t ) {
        if ( $D->[$t][$RC_LEN] < $bestscore ) {
            $bestscore = $D->[$t][$RC_LEN];
            @matches   = ($t);
        }
        elsif ( $D->[$t][$RC_LEN] == $bestscore ) {
            push( @matches, $t );
        }
    }

    return $bestscore;
}

####################################
sub min3 {
    my ( $i, $j, $k ) = @_;
    my ($tmp);

    $tmp = ( $i < $j ? $i : $j );
    $tmp < $k ? $tmp : $k;
}

####################################
sub align_muscle {
    my $dir = cwd();

    foreach my $fp ( glob("$dir/*.fa") ) {
        open my $fh, "<", $fp or die;

        my @params = ( 'IN' => "$fp", 'OUT' => "$fp.fasta", 'MAXITERS' => 1 );

        my $factory = Bio::Tools::Run::Alignment::Muscle->new(@params);

        my $aln = $factory->align( my $inputfilename );

    }

}

####################################
sub align_convert {
    my $dir = cwd();

    foreach my $fp ( glob("$dir/*.fasta") ) {
        open my $fh, "<", $fp or die;

        my $in = Bio::AlignIO->new(
            -file   => "$fp",
            -format => 'fasta'
        );
        my $out = Bio::AlignIO->new(
            -file   => ">$fp.aln",
            -format => 'clustalw'
        );

        while ( my $aln = $in->next_aln ) {
            $out->write_aln($aln);
        }

    }
}

####################################
sub graphics {
    my $dir = cwd();

    foreach my $fp ( glob("$dir/*.fa.fasta.aln") ) {
        open my $fh, "<", $fp or die;

        my $output = "$fp.jpeg";

        #Create an AlignI object using AlignIO
        my $in = new Bio::AlignIO( -file => $fp, -format => 'clustalw' );

        #Read the alignment
        my $aln = $in->next_aln();

        my $print_align = new Bio::Align::Graphics(
            align              => $aln,
            output             => $output,
            font               => 5,
            x_label            => "true",
            y_label            => "true",
            bg_color           => "white",
            font_color         => "black",
            x_label_color      => "red",
            y_label_color      => "blue",
            pad_top            => 5,
            pad_bottom         => 5,
            pad_left           => 5,
            pad_right          => 5,
            x_label_space      => 1,
            y_label_space      => 1,
            reference_id       => "First sequence supplied in alignment",
            block_size         => 20,
            block_space        => 2,
            show_nonsynonymous => 0,
            out_format         => "jpeg"
        );

        $print_align->draw();

    }

}

####################################
sub graphics_all_primers {
    my $dir = cwd();
    my $orientation;
    my @gene_length;
    my @gene_id;
    my @unique_length = uniq(@array_length);
    my @unique_name   = uniq(@array_name);
    my %hash;
    my $len;
    my $name_id;

    foreach my $unique_name (@unique_name) {
        if ( $unique_name =~ m/^\d+'(.+)?/i ) {
            push @gene_id, $1;
        }
    }

    foreach my $unique_length (@unique_length) {
        if ( $unique_length =~ m/^.+'(\d+)?/i ) {
            push @gene_length, $1;
        }
    }

    foreach my $fp ( glob("$dir/GRAPHIC_*.txt") ) {
        open my $fh, "<", $fp or die;

        @hash{@gene_id} = @gene_length;
        while ( my ( $key, $val ) = each %hash ) {
            $len = $val if $fp =~ m/$key/i;
        }
        while ( my ( $key, $val ) = each %hash ) {
            $name_id = $key if $fp =~ m/$key/i;
        }

        my $outputfile = "$fp.png";
        open( OUTFILE, ">$outputfile" ) or die;

        my $panel = Bio::Graphics::Panel->new(
            -length     => $len,
            -width      => 800,
            -pad_left   => 30,
            -pad_right  => 30,
            -pad_top    => 20,
            -pad_bottom => 20,
        );

        my $full_length = Bio::SeqFeature::Generic->new(
            -start        => 1,
            -end          => $len,
            -display_name => $name_id,
        );

        $panel->add_track(
            $full_length,
            -glyph        => 'arrow',
            -tick         => 2,
            -fgcolor      => 'black',
            -double       => 1,
            -label        => 1,
            -strand_arrow => 1,
        );

        my $track = $panel->add_track(
            -glyph        => 'transcript2',
            -label        => 1,
            -strand_arrow => 1,
            -bgcolor      => 'blue',
            -bump         => +1,
            -height       => 12,
        );

        while (<$fh>) {
            chomp;

            my ( $name, $score, $start, $end ) = split /\t/;
            if ( $start - $end < 1 ) {
                $orientation = +1;
            }
            elsif ( $start - $end > 1 ) {
                $orientation = -1;
            }
            my $feature = Bio::SeqFeature::Generic->new(
                -display_name => $start,
                -score        => $score,
                -start        => $start,
                -strand       => $orientation,
                -end          => $end
            );
            $track->add_feature($feature);

        }

        binmode OUTFILE;

        print OUTFILE $panel->png;

    }

}

####################################
sub uniq {
    my %seen;
    return grep { !$seen{$_}++ } @_;
}

####################################
sub clean_up {
    my $dir = cwd();
    unlink glob "$dir/*.fa";
    unlink glob "$dir/*.fa.fasta";
    unlink glob "$dir/*.fa.fasta.aln";
}

####################################
####################################
1;




=pod
 
=head1 INPUT PARAMETERS
 
PrimerView
Author: Damien O'Halloran, The George Washington University, 2015

To run simply execute as follows specifying the Getopt arguments: 
>perl primerview_driver.pl [-a filename e.g. test_seqs.fasta] 
   [-b 5' search area, integer] [-c 3' search area, integer] 
[-d primer max, integer] [-e primer min, integer] [-f GC clamp Y or N]
[-g upper GC, integer] [-h lower GC, integer] [-i upper Tm, integer] 
[-j lower Tm, integer] 

example settings: ">perl primerview_driver.pl -a test_seqs.fasta"

default settings are as follows:

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

the defaults will be overwritten if a commandline parameter is added
   
where:  'primerview_driver.pl' is a simple driver for the primerview package, 
and 'test_seqs.fasta' is a fasta file with sample query sequences.

if only the distribution graphic of primers is required
then cancel the following subs:
$tmp->align_muscle();
$tmp->align_convert();
$tmp->graphics();

you must have muscle.exe in your PATH and you must have following BioPerl modules: 
##  Bio::SeqIO 
##  Bio::Tools::Run::Alignment::Muscle 
##  Bio::AlignIO
##  Bio::Align::Graphics
##  Bio::Graphics
##  Bio::SeqFeature::Generic

WARNING: the subroutine 'clean_up' deletes the '.fa', '.fa.fasta', 
and '.fa.fasta.aln' extension files generated from cwd, 

=cut

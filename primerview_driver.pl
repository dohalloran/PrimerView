#!/usr/bin/perl

use strict;
use warnings;
use PRIMERVIEW;

MAIN: {
    #create an instance of PRIMERVIEW
    my $tmp = PRIMERVIEW->new();
    $tmp->primerview();
    $tmp->align_muscle();
    $tmp->align_convert();
    $tmp->graphics();
    $tmp->graphics_all_primers();
    $tmp->clean_up();
}

=pod
 
=head1 DESC
 
if only the distribution graphic of primers is required
then cancel the following subs:

$tmp->align_muscle();
$tmp->align_convert();
$tmp->graphics();

=cut
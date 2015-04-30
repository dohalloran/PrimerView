PrimerView
Author: Damien O'Halloran, The George Washington University, 2015

To start:

1. Download and extract the primerview.zip file
2. >tar -xzvf primerview.zip 
3. The extracted dir will be called PRIMERVIEW  
4. cd PRIMERVIEW
5. insert sequence file (FASTA format) 
[folder contains a sample file called 'test_seqs.fasta']
6. To run simply execute as follows specifying the Getopt arguments 
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


Any questions:
email: stitcher[AT]ohalloranlab.net
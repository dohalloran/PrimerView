##PrimerView
Author: Damien O'Halloran, The George Washington University, 2015

##Installation
1. Download and extract the primerview.zip file  
`tar -xzvf primerview.zip`  
2. The extracted dir will be called PRIMERVIEW  
  `cd PRIMERVIEW`  
  `perl Makefile.PL`  
  `make`  
  `make test`  
  `make install`  
  
##Getting Started  
1. You must have `muscle.exe` in your PATH  
[Click here to get MUSCLE] (http://www.drive5.com/muscle/)  
2. Must have following [BioPerl] (https://github.com/bioperl) modules:  
Bio::SeqIO  
Bio::Tools::Run::Alignment::Muscle  
Bio::AlignIO  
Bio::Align::Graphics  
Bio::Graphics  
Bio::SeqFeature::Generic  
WARNING: the subroutine 'clean_up' deletes the '.fa', '.fa.fasta', and '.fa.fasta.aln' extension files generated from cwd  
3. Start with a sequence file in FASTA format (for example see `test_seqs.fasta`)  

##Usage 
Run as follows:  
  `perl primerview_driver.pl`  
  
##GetOpts 
  `-a test_seqs.fasta` //filename  
   `-b 200` //5' search area    
   `-c 200` // 3' search area    
   `-d 28` //primer max    
   `-e 20` //primer min   
   `-f Y` // GC clamp Y or N  
   `-g 60` //upper GC%  
   `-h 40` //lower GC%  
   `-i 70` //upper Tm  
   `-j 55` // lower Tm  
   `-k Y` //specificty across entire input file, Y or N  

## Contributing
All contributions are welcome.

## Support
If you have any problem or suggestion please open an issue [here](https://github.com/dohalloran/PrimerView/issues).

## License 
GNU GENERAL PUBLIC LICENSE






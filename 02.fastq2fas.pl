#!/bin/perl -w

use strict;
my $filein =$ARGV[0];
my $fileout=$ARGV[1];
my $line;
my @temp;
my %array;
my $i;
my $find;
my $probe;


	open(SRC,"$filein")||die;
	open(TGT,">$fileout")||die "can not WRITE file:";
	

 while ($line=<SRC>){
    chomp $line;
    if($line=~/^\@/){
    	 $line='>'.(substr $line,1);
    	 my @id=split(/\s+/,$line);    	 
    	 my $pop_id =shift @id;
    	 $line=<SRC>;
    	 chomp $line;
    	 my $len=length($line);
    	
    	
    	 if($len>20){
    	 	 print TGT $pop_id,'_',$len,"\n";    	 
    	     print TGT "$line\n";
    	 }
    	
    		
    		  $line=<SRC>;
    		#print TGT $line;
    		  $line=<SRC>;
    		#print TGT $line;
    	

    }

 }
     
   
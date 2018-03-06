#!/usr/bin/perl -w
# run the script: time perl 00.script/00.download.sra.pl 01.data/00.PriorData/SRA174408 zcluster

use strict;
system("echo 'Running 00.download.sra.pl ....' >> job.monitor.txt");

my $reffile = shift @ARGV;
my $mode = shift @ARGV;
my $platform = lc(shift @ARGV);

my @temp = split(/\//, $reffile);
my $sample = pop @temp;
my $dir = "00.script/00.download.script/$sample";

open(SRC, $reffile);
system("mkdir -p $dir");


foreach my $line (<SRC>){
	my @lines = split(/\s+/, $line);
	my $sra = $lines[0];
	my $sub = $lines[1];
	my $head = substr $sub, 0, 6;
	system("mkdir -p 01.data/01.Fastq/$sra");
	
    my $shell = "$dir/download.$sra.$sub.sh";
    open(SHL, ">$shell");
    if($platform eq "sapelo"){
	    print SHL "#PBS -S /bin/bash\n";
	    print SHL "#PBS -q batch\n";
	    print SHL "#PBS -N download.$sra.$sub.sh\n";
	    print SHL "#PBS -l nodes=1:ppn=1:AMD\n";
	    print SHL "#PBS -l walltime=12:00:00\n";
	    print SHL "#PBS -l mem=2gb\n";
	    print SHL "\n";
	    print SHL "module purge\n";
	    print SHL "module load sratoolkit/2.5.7\n";
	    print SHL "cd \$PBS_O_WORKDIR\n\n";
    }elsif($platform eq "zcluster"){
	    print SHL "#!/bin/bash\n\n";
    }else{
	    die "Please provide the platform: 'Sapelo' or 'Zcluster'";
    }
	if($platform eq "sapelo"){
        print SHL "wget -r -nHd ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/SRX/$head/$sub --directory-prefix=01.data/01.Fastq/$sra/$sub\n";
        print SHL "fastq-dump -O 01.data/01.Fastq/$sra/$sub 01.data/01.Fastq/$sra/$sub/* --split-3\n";
        if($mode eq "single-end"){
            print SHL "cat 01.data/01.Fastq/$sra/$sub/*.fastq > 01.data/01.Fastq/$sra/$sub/$sub.fastq\n";
        }elsif($mode eq "paired-end"){
            print SHL "cat 01.data/01.Fastq/$sra/$sub/*_1.fastq > 01.data/01.Fastq/$sra/$sub/$sub.R1.fastq\n";
            print SHL "cat 01.data/01.Fastq/$sra/$sub/*_2.fastq > 01.data/01.Fastq/$sra/$sub/$sub.R2.fastq\n";
        }else{
            die "Please provide the mode: 'single-end' or 'paired-end'";
        }

		print SHL "ls -l 01.data/01.Fastq/$sra/$sub/\n";
        #print SHL "rm 01.data/01.Fastq/$sra/$sub/*.sra\n";
	    print SHL "touch $dir/done.log\n";
	}elsif($platform eq "zcluster"){
        #print SHL "wget -r -nHd ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/SRX/$head/$sub --directory-prefix=01.data/01.Fastq/$sra/$sub\n";
        print SHL "time /usr/local/aspera/latest/bin/ascp -i /usr/local/aspera/latest/etc/asperaweb_id_dsa.putty -L $dir -k 1 -QTrd -l 200m anonftp\@ftp-trace.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByExp/sra/SRX/$head/$sub 01.data/01.Fastq/$sra/\n\n";
        print SHL "for r in 01.data/01.Fastq/$sra/$sub/*\ndo\n";
        print SHL "\tf=\${r##*/}\n";
        print SHL "\ttime /usr/local/sra/2.3.4-2/bin/fastq-dump -O 01.data/01.Fastq/$sra/$sub 01.data/01.Fastq/$sra/$sub/\$f/\$f.sra --split-3\ndone\n\n";
        if($mode eq "single-end"){
            print SHL "cat 01.data/01.Fastq/$sra/$sub/*.fastq > 01.data/01.Fastq/$sra/$sub/$sub.fastq\n";
        }elsif($mode eq "paired-end"){
            print SHL "cat 01.data/01.Fastq/$sra/$sub/*_1.fastq > 01.data/01.Fastq/$sra/$sub/$sub.R1.fastq\n";
            print SHL "cat 01.data/01.Fastq/$sra/$sub/*_2.fastq > 01.data/01.Fastq/$sra/$sub/$sub.R2.fastq\n";
        }else{
            die "Please provide the mode: 'single-end' or 'paired-end'";
        }
        print SHL "ls -l 01.data/01.Fastq/$sra/$sub/\n";
        #print SHL "rm 01.data/01.Fastq/$sra/$sub/*.sra\n";
	    print SHL "touch $dir/done.log\n";
	}else{
		die "Please provide the platform: 'Sapelo' or 'Zcluster'";
	}

    close SHL;
    system("chmod u+x $shell");
    if($platform eq "sapelo"){
	    system("qsub $shell");
    }elsif($platform eq "zcluster"){
        system("qsub -q rcc-30d $shell");
    }else{
	    die "Please provide the platform: 'Sapelo' or 'Zcluster'";
    }
}


close SRC;
system("echo 'Finished 00.download.sra.pl!' >> job.monitor.txt");

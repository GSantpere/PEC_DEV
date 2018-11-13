##AIM:	Run STAR alignment

use Getopt::Std;

my %opt=();
getopts("b:s:",\%opt);
my $brainName=$opt{b} if defined($opt{b}) || die "Please Input brain Name! -b\n";
my $regionName=$opt{s} if defined($opt{b}) || die "Please Input region Name! -s\n";
print "$brainName\t$regionName\n";

##-----set up output folder
my $brainFold="../$brainName";
mkdir $brainFold unless (-e $brainFold);
my $outFold="../$brainName/$regionName";
mkdir $outFold unless (-e $outFold);

##---set path
my $inFold = "./$brainName/$regionName";
my $indexDir="./STAR_2.4.0e/hg38ANDspikeIn";
my $comDir="./STAR_2.4.0e";
my $inFile=join(".",$brainName,$regionName,"fq");
my $sampName=join("",$brainName,".",$regionName,".");

##---STAR hg38
system("$comDir/STAR --runMode alignReads --readFilesIn $inFold/$inFile --outFileNamePrefix $outFold/$sampName --genomeDir $indexDir --runThreadN  8 --outSAMattributes  All --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM  62000000000 --quantMode  TranscriptomeSAM --outFilterMismatchNoverLmax 0.1 --alignSJoverhangMin  8 --alignSJDBoverhangMin  1 --outSAMunmapped  Within --outFilterType  BySJout --alignMatesGapMax  500 --outFilterMultimapNmax  50 --alignEndsType  Local --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical");

my $now=localtime;
print "All Jobs Finished!:\t$now\n";



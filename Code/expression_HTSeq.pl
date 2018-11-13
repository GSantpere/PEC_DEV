##AIM:	get reads count per gene by using HTSeq


use Getopt::Std;

my %opt=();
getopts("b:s:",\%opt);
my $brainName=$opt{b} if defined($opt{b}) || die "Please Input brain Name! -b\n";
my $regionName=$opt{s} if defined($opt{s}) || die "Please Input region Name! -s\n";
print "$brainName\t$regionName\n";

##---set data folder
my $sampFold ="./$brainName/$regionName";
chdir("$sampFold");
	
##------convert bam file to sam file
my $inFile = join(".",$brainName,$regionName,"Aligned.sortedByCoord.out.bam");
my $outFile = "myTmp.sam";
system("samtools view $inFile > $outFile");

##-----remove multiple mapping
$inFile="myTmp.sam";
$outFile=join(".",$brainName,$regionName,"myTmp2.sam");
open INF, $inFile or die $!;
open OUTF2,">$outFile" or die $!;
while(<INF>){
	chomp;
	my @line=split(/\t/,$_);
	my $flag=0;
	for(my $i=11;$i<@line;$i++){
		$flag=1 if($line[$i] eq "NH:i:1");
		
		##----change for STAR
		$line[$i]="" if($line[$i] =~ "jM|jI");
	}
	if($flag==1){
		my $newline = join("\t",@line);
		print OUTF2 "$newline\n";   ##just for STAR 
	}
}
close(INF);
close(OUTF2);
system("rm myTmp.sam");


##-----get reads map to exons
$outFile2 = join(".",$brainName,$regionName,"gene.count.txt");
my $cmd = "./HTSeq-0.6.1/scripts/htseq-count";
my $annot = "./gencode.v21.annotation.forHTSeq.gtf";
system("python $cmd -s no -t exon  $outFile $annot > $outFile2");
system("rm $outFile");

##------back to source code folder
chdir("/home/ml724/program/brainSpan");

	
my $now=localtime;
print "All Jobs Finished!:\t$now\n";





##AIM:	use RSEQTools to get gene RPKM expression


#############################------------
##AIM:	get mrf file from bam file
use Getopt::Std;

my %opt=();
getopts("b:s:",\%opt);
my $brainName=$opt{b} if defined($opt{b}) || die "Please Input brain Name! -b\n";
my $regionName=$opt{s} if defined($opt{s}) || die "Please Input region Name! -b\n";
print "$brainName\t$regionName\n";

##---recording multiple reads
open OUTF5,">>./brainSpan.multRead.txt" or die $!;

##---set data folder
my $sampFold ="./$brainName/$regionName";
chdir("$sampFold");
	
##------convert bam file to sam file
my $inFile = join(".",$brainName,$regionName,"Aligned.sortedByCoord.out.bam");
my $outFile="myTmp.sam";
system("samtools view $inFile > $outFile");

##-----remove multiple mapping
my %multRead=();
$inFile="myTmp.sam";
$outFile="myTmp2.sam";
open INF, $inFile or die $!;
open OUTF2,">$outFile" or die $!;
while(<INF>){
	chomp;
	my @line=split(/\t/,$_);
	my $flag=0;
	for(my $i=11;$i<@line;$i++){
		$flag=1 if($line[$i] eq "NH:i:1");
	}
	if($flag==1){
		print OUTF2 "$_\n"; 
	}
	else{
		$multRead{$line[0]}=1;
	}	
}
close(INF);
close(OUTF2);
system("rm myTmp.sam");

##---count multiple hitting reads
my $ct_multRead=keys %multRead;
print OUTF5 "$brainName\t$regionName\t$ct_multRead\n";

##-----convert sam file to mrf file
$outFile=join(".",$brainName,$regionName,"mrf");
system("./python2.7/bin/python sam2mrf.py  myTmp2.sam  myTmp3.mrf");
system("cut -f 1 myTmp3.mrf > $outFile");
system("rm myTmp2.sam myTmp3.mrf");

##------back to source code folder
chdir("/home/ml724/program/brainSpan");

	
my $now=localtime;
print "All Jobs Finished!:\t$now\n";



############################-----------
##AIM:	generate coverage,junction, insertion and deletion file in BigWig and BigBed format
use Getopt::Std;

my %opt=();
getopts("b:",\%opt);
my $brainName=$opt{b} if defined($opt{b}) || die "Please Input brain Name! -b\n";

##---set data folder
my $brainFold="./$brainName";
my @regionList=`ls $brainFold`;
my $chromFile="./hg38.chromInfo.txt";
open OUTF, ">>./brainSpan.bigwig.track.txt" or die $!;

##----read files
for(my $i=0;$i<@regionList;$i++){
	my $regionName=$regionList[$i];
	chomp($regionName);
	print "$brainName\t$regionName\n";
	
	#---change work folder
	chdir("$brainFold/$regionName");
	
	##---remove chrM reads
	my $inFile=join(".",$brainName,$regionName,"mrf");
	system("grep -v chrM $inFile > myTmp_Bigfile.mrf");
	
	##------generate wig files
	my $sampName=join(".",$brainName,$regionName);
	$inFile="myTmp_Bigfile.mrf";
	print "$inFile\n";
	my $trackName=$sampName;
	system("mrf2bedGraph_MFL $trackName < $inFile");
	system("rm myTmp_Bigfile.mrf");


	##------generate bigwig files for coverage
	$inFile=join(".",$sampName,"wig");
	print "$inFile\n";
	my $outFile=join(".",$sampName,"bw");
	system("wigToBigWig $inFile $chromFile $outFile");
	my $name=join(".",$sampName,"coverage");
	my $output="track type=bigWig name=\"$name\" description=\"$name\" bigDataUrl=./$outFile visibility=full";
	print OUTF "$output\n";
	#system("rm $inFile");
		
	##---change back home folder
	chdir("/home/ml724/program/brainSpan/");
}

my $now=localtime;
print "All Jobs Finished!:\t$now\n";


#########
##AIM:	get gene/exon RPKM and reads counts
use Getopt::Std;

my %opt=();
getopts("b:",\%opt);
my $brainName=$opt{b} if defined($opt{b}) || die "Please Input brain Name! -b\n";

##---set data folder
my $annotFold="./";
my $brainFold="./$brainName";
my @regionList=`ls $brainFold`;
my $stat="singleOverlap";

##----reading files
for(my $i=0;$i<@regionList;$i++){
	my $regionName=$regionList[$i];
	chomp($regionName);
	print "$brainName\t$regionName\n";
	
	#---change work folder
	chdir("$brainFold/$regionName");
	
	my $inFile;
	my $outFile;
	my $annotFile;	
	my $cmd;

	##-----change back
	$cmd="mrfQuantifier_MFL75";

	##------remove chrM and ERCC
	$inFile=join(".",$brainName,$regionName,"mrf");
	system("grep -v chrM $inFile | grep -v ERCC  > myTmp_getExpression.mrf");
	

	##------wholeGene.exonComposite
	$annotFile="gencode.v21.wholeGene.exonComposite.interval";
	$outFile=join(".",$brainName,$regionName,"wholeGene.exonComposite.exon.expression");
	print "$annotFile\n";
	system("$cmd $annotFold/$annotFile $stat < myTmp_getExpression.mrf > $outFile");		
	
	##------wholeGene.geneComposite
	$annotFile="gencode.v21.wholeGene.geneComposite.interval";
	$outFile=join(".",$brainName,$regionName,"wholeGene.geneComposite.gene.expression");
	print "$annotFile\n";
	system("$cmd $annotFold/$annotFile $stat < myTmp_getExpression.mrf > $outFile");	
	
	##------delete temporal file
	system("rm myTmp_getExpression.mrf");
	
	##------back to source code folder
	chdir("/home/ml724/program/brainSpan");
}	
my $now=localtime;
print "All Jobs Finished!:\t$now\n";



##################
##AIM:	get expression of spike-in




use Getopt::Std;

my %opt=();
getopts("b:",\%opt);
my $brainName=$opt{b} if defined($opt{b}) || die "Please Input brain Name! -b\n";

##---set data folder
my $annotFold="./";
my $brainFold="./$brainName";
my @regionList=`ls $brainFold`;
my $stat="singleOverlap";

##----reading files
for(my $i=0;$i<@regionList;$i++){
	my $regionName=$regionList[$i];
	chomp($regionName);
	print "$brainName\t$regionName\n";
	
	##---debug
	#next if($regionName ne "AMY");

	#---change work folder
	chdir("$brainFold/$regionName");
	
	my $inFile;
	my $outFile;
	my $annotFile;	
	my $cmd;

	##-----change back
	$cmd="mrfQuantifier_MFL75";

	##------remove chrM and ERCC
	$inFile=join(".",$brainName,$regionName,"mrf");
	system("grep -v chrM $inFile > myTmp_getExpression.mrf");
	

	##------get spike-in expression files
	$outFile=join(".",$brainName,$regionName,"spike.expression");
	$annotFile="spike_in_subset.interval";
	print "$annotFile\n";
	system("$cmd $annotFold/$annotFile $stat < myTmp_getExpression.mrf > $outFile");
	
	##------delete temporal file
	system("rm myTmp_getExpression.mrf");
	
	##------back to source code folder
	chdir("/home/ml724/program/brainSpan");
}	
my $now=localtime;
print "All Jobs Finished!:\t$now\n";











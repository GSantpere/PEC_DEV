##AIM: conduct gene splicing analysis


##----generate annotation
python ~/package/R/R-3.1.0/lib64/R/library/DEXSeq/python_scripts/dexseq_prepare_annotation.py gencode.v21.annotation.forHTSeq.gtf gencode.v21.annotation_reduce.gtf
awk '{OFS="\t"}{if ($3 == "exonic_part") print $1,$2,$3,$4,$5,$6,$7,$8,$14":"$12}' gencode.vM1.annotation_reduce.gtf | sed 's=[";]==g' > gencode.v21.annotation.forDEXSeq.collapse.exonUniqID.gff
rm gencode.v21.annotation_reduce.gtf


#####################--------------------
##AIM:	calculate PSI for brainSpan
use Getopt::Std;

my %opt=();
getopts("b:s:",\%opt);
my $brainName = $opt{b} if defined($opt{b}) || die "Please Input brain Name! -b\n";
my $regionName = $opt{s} if defined($opt{b}) || die "Please Input region Name! -s\n";
print "$brainName\t$regionName\n";

##-----set up output folder
my $outFold="./$brainName/$regionName/PSI";
mkdir $outFold unless (-e $outFold);

##---set path
my $inFold = "./$brainName/$regionName";
my $annotFile ="./gencode.v21.annotation.forDEXSeq.collapse.exonUniqID.gff";
my $bamFile=join(".",$brainName,$regionName,"Aligned.sortedByCoord.out.bam");
my $sjFile=join(".",$brainName,$regionName,"SJ.out.tab");
my $comDir = "/home/ml724/soft/bedtools-2.17.0/bin";


###-----run com part_1
chdir("$inFold");
system("$comDir/coverageBed -split -abam $bamFile -b $annotFile | awk \'BEGIN{OFS=\"\\t\"} {print \$1,\$4,\$5,\$5-\$4+1,\$9,\$10}\' | sort -k 5 > PSI/exonic_parts.inclusion");
system("awk \'BEGIN{OFS=\"\\t\"}{print \$1, \$2-20-1, \$3+20,\"JUNCBJ\"NR, \$7, (\$4 == 1)? \"+\":\"-\",\$2-20-1, \$3+20,\"255,0,0\", 2, \"20,20\", \"0,300\" }\' $sjFile > PSI/junctions.bed");

##----run com part_2
chdir("$outFold");
system("sed \'s/,/\t/g\' junctions.bed | awk \'BEGIN{OFS=\"\\t\"} {print \$1,\$2,\$2+\$13,\$4,\$5,\$6}\' > left.bed");
system("sed \'s/,/\t/g\' junctions.bed | awk \'BEGIN{OFS=\"\\t\"} {print \$1,\$3-\$14,\$3,\$4,\$5,\$6}\' > right.bed");
system("$comDir/intersectBed -u -s -a left.bed -b $annotFile > left.overlap");
system("$comDir/intersectBed -u -s -a right.bed -b $annotFile > right.overlap");
system("cat left.overlap right.overlap | cut -f4 | sort |uniq -c | awk \'{ if(\$1 == 2) print \$2 }\' > filtered_junctions.txt");
system("grep -F -f filtered_junctions.txt junctions.bed > filtered_junctions.bed");
system("sed \'s/,/\t/g\' filtered_junctions.bed | grep -v description | awk \'{OFS=\"\\t\"}{print \$1,\$2+\$13,\$3-\$14,\$4,\$5,\$6}\' > intron.bed");
system("rm filtered_junctions.bed");
system("$comDir/intersectBed -wao -f 1.0 -s -a $annotFile -b intron.bed | awk \'BEGIN{OFS=\"\\t\"}{\$16== 0? s[\$9] += 0:s[\$9] += \$14}END{for (i in s) {print i,s[i]}}\' | sort -k 1 > exonic_parts.exclusion");
system("paste exonic_parts.inclusion exonic_parts.exclusion | awk -v \"len=75\" \'BEGIN{OFS=\"\\t\"; print \"exon_ID\" , \"length\" , \"inclusion\" , \"exclusion\" ,\"PSI\"}{NIR=\$6/(\$4+len-1) ; NER=\$8/(len-1)}{print \$5,\$4,\$6,\$8,(NIR+NER<=0)? \"NA\":NIR / (NIR + NER)}\' > exonic_parts.psi");

##----delet files
system("rm left.bed right.bed left.overlap right.overlap filtered_junctions.txt");
system("rm intron.bed");
system("rm exonic_parts.inclusion");
system("rm exonic_parts.exclusion");

##------change folder
chdir("/home/ml724/program/geneSplicing");

####---------------
my $now=localtime;
print "All Jobs Finished!:\t$now\n";



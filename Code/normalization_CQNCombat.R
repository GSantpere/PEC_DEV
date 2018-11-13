##AIM:	do conditional quantile normalization and combat correction for brainspan


library("cqn");
source("ComBat.R");
source("cqn2.R");
options(stringsAsFactors = F);


##------gene expression
genexp = read.table("./brainSpan.hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.count.txt",header=T);
geneList = as.character(genexp[,1]);
genexp = as.matrix(genexp[, -1]);
sampList = colnames(genexp);


##------gene annotation
geneAnnot = read.table("./gencode.v21.gene.annot.txt");
myindex = which(as.character(geneAnnot[,3]) == "chrM");
geneList_chrM = as.character(geneAnnot[myindex,1]);

##------assigne chrM gene to 0
myindex = match(geneList_chrM, geneList);
genexp[myindex,] = 0;


##---cqn conditional quantile normalization
##---gene info
geneInfo = read.table("./gencode.v21.wholeGene.geneComposite.geneGCcontent.txt",header=T,sep="\t");
myindex = match(geneList, as.character(geneInfo[,1]));
geneInfo = geneInfo[myindex,];

##---sequencing depth info(removed chrM reads)
mydata = read.table("./brainSpan.chrom.read.xls",header=T,sep="\t");
seqDepth = as.numeric(mydata[,3])*(1-mydata[,ncol(mydata)]/100);
seqDepth = round(seqDepth);
myindex = match(sampList,paste(as.character(mydata[,1]),as.character(mydata[,2]),sep="."));
seqDepth = seqDepth[myindex];

##---cqn normlization
backupGenexp = genexp;
cqn.genexp = cqn2(genexp,lengths = as.numeric(geneInfo[,2]), x= as.numeric(geneInfo[,3]),
					sizeFactors=seqDepth,lengthMethod="smooth",sqn=T);					
RPKM.cqn = cqn.genexp$y + cqn.genexp$offset;
genexp = RPKM.cqn;
##Notes: it is log2 and quantile normalized

##---save CQN'ed RPKM
res=cbind(geneList,genexp);
colnames(res)=c("Geneid",sampList);
write.table(res,file="./brainSpan.hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.RPKM.log2.normalized.CQN.txt",sep="\t",quote=F,row.name=F,col.names=T,append=F);

##----save CQN'ed count
res=cbind(geneList,genexp);
colnames(res)=c("Geneid",sampList);
write.table(res,file="./brainSpan.hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.count.log2.normalized.CQN.txt",sep="\t",quote=F,row.name=F,col.names=T,append=F);


##---save CQN'ed RPM
res=cbind(geneList,genexp);
colnames(res)=c("Geneid",sampList);
write.table(res,file="./brainSpan.hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.RPM.log2.normalized.CQN.txt",sep="\t",quote=F,row.name=F,col.names=T,append=F);


######################---batch effect correction by combat correction
if(file.exists("./brainSpan.info.site.txt")) file.remove("./brainSpan.info.site.txt");
cat("array name","Sample name", "Batch",file="./brainSpan.info.site.txt",sep="\t",append=T);
cat("\n",file="./brainSpan.info.site.txt",sep="\t",append=T);
for(i in 1:length(sampList)){
	sampName=sampList[i];
	
	##--- set site info
	##---yale brains
	if(length(grep("HSB123",sampName,fixed=T)) >0) site=1;
	if(length(grep("HSB126",sampName,fixed=T)) >0) site=1;
	if(length(grep("HSB130",sampName,fixed=T)) >0) site=1;
	if(length(grep("HSB145",sampName,fixed=T)) >0) site=1;
	if(length(grep("HSB114",sampName,fixed=T)) >0) site=1;
	if(length(grep("HSB122",sampName,fixed=T)) >0) site=1;
	if(length(grep("HSB124",sampName,fixed=T)) >0) site=1;
	if(length(grep("HSB141",sampName,fixed=T)) >0) site=1;
	if(length(grep("HSB149",sampName,fixed=T)) >0) site=1;
	if(length(grep("HSB153",sampName,fixed=T)) >0) site=1;
	if(length(grep("HSB159",sampName,fixed=T)) >0) site=1;
	if(length(grep("HSB175",sampName,fixed=T)) >0) site=1;
	if(length(grep("HSB194",sampName,fixed=T)) >0) site=1;
	if(length(grep("HSB121",sampName,fixed=T)) >0) site=1;
	if(length(grep("HSB127",sampName,fixed=T)) >0) site=1;
	if(length(grep("HSB139",sampName,fixed=T)) >0) site=1;
	if(length(grep("HSB148",sampName,fixed=T)) >0) site=1;
	if(length(grep("HSB150",sampName,fixed=T)) >0) site=1;
	if(length(grep("HSB155",sampName,fixed=T)) >0) site=1;
	if(length(grep("HSB173",sampName,fixed=T)) >0) site=1;
	if(length(grep("HSB178",sampName,fixed=T)) >0) site=1;
	if(length(grep("HSB96",sampName,fixed=T))  >0) site=1;
	                                         
	##---use brains                          
	if(length(grep("HSB135",sampName,fixed=T)) >0) site=2;
	if(length(grep("HSB136",sampName,fixed=T)) >0) site=2;
	if(length(grep("HSB103",sampName,fixed=T)) >0) site=2;
	if(length(grep("HSB105",sampName,fixed=T)) >0) site=2;
	if(length(grep("HSB107",sampName,fixed=T)) >0) site=2;
	if(length(grep("HSB112",sampName,fixed=T)) >0) site=2;
	if(length(grep("HSB113",sampName,fixed=T)) >0) site=2;
	if(length(grep("HSB118",sampName,fixed=T)) >0) site=2;
	if(length(grep("HSB119",sampName,fixed=T)) >0) site=2;
	if(length(grep("HSB131",sampName,fixed=T)) >0) site=2;
	if(length(grep("HSB132",sampName,fixed=T)) >0) site=2;
	if(length(grep("HSB143",sampName,fixed=T)) >0) site=2;
	if(length(grep("HSB154",sampName,fixed=T)) >0) site=2;
	if(length(grep("HSB168",sampName,fixed=T)) >0) site=2;
	if(length(grep("HSB169",sampName,fixed=T)) >0) site=2;
	if(length(grep("HSB171",sampName,fixed=T)) >0) site=2;
	if(length(grep("HSB172",sampName,fixed=T)) >0) site=2;
	if(length(grep("HSB174",sampName,fixed=T)) >0) site=2;
	if(length(grep("HSB92",sampName,fixed=T)) >0) site=2;
	if(length(grep("HSB97",sampName,fixed=T)) >0) site=2;
	if(length(grep("HSB98",sampName,fixed=T)) >0) site=2;

	##---output
	cat(sampName,sampName,site,file="./brainSpan.info.site.txt",sep="\t",append=T);
	cat("\n",file="./brainSpan.info.site.txt",sep="\t",append=T);
}

write.table(genexp,file="tmp.genexp.txt",col.names=T,row.names=T,quote=F,sep="\t");
combatGenexp = ComBat("tmp.genexp.txt","./brainSpan.info.site.txt",skip=1,write=F,par.prior=T,prior.plots=F);
if(file.exists("tmp.genexp.txt")) file.remove("tmp.genexp.txt");
combatGenexp  = data.matrix(combatGenexp);

##----back to RPKM scale
##Notes: very important
combatGenexp = 2**combatGenexp ;
combatGenexp = round(combatGenexp,6); ##slim data

##----output normlized gene expression
res = cbind(geneList,combatGenexp);
colnames(res) = c("Geneid",sampList);
write.table(res,file="./brainSpan.hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.RPKM.normalized.CQNCombat.txt",sep="\t",quote=F,row.name=F,col.names=T,append=F);









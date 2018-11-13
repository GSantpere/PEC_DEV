##AIM: Conduct spatial or temporal differentially expression analysis for bulk tissue mRNA-seq using DESeq2 package


###################------------------------
#####temporal differentially expressed genes
library("DESeq2");
library("cqn");
source("cqn2.R");

##-----------get gene annotation
geneAnnot = read.table("data/gencode.v21.gene.annot.txt");

##---gene info
geneInfo = read.table("data/gencode.v21.wholeGene.geneComposite.geneGCcontent.txt",header=T,sep="\t");
 
##---get gene RPKM
geneRPKM = read.table("result/genexp/brainSpan.hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.RPKM.txt",header=T,sep="\t");

##---get gene count
geneCount = read.table("result/genexp/brainSpan.hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.count.txt",header=T,sep="\t");
sampList = colnames(geneCount);
geneList = as.character(geneCount[,1]);

##----copy some samples
myindex = grep("MSC", sampList,fixed=T);
addS1C = gsub("MSC","S1C",sampList[myindex]);
geneRPKM = cbind(geneRPKM, geneRPKM[,myindex]);
geneCount = cbind(geneCount, geneCount[,myindex]);
geneList = as.character(geneCount[,1]);  ###
sampList = c(sampList,addS1C) ;          ###
sampList = gsub("MSC","M1C", sampList);
colnames(geneRPKM) = sampList;
colnames(geneCount) = sampList;

####----------------seq
##-----re-name region name
sampList2 = sampList[-1];
rec = unlist(strsplit(sampList2,split=".",fixed=T))                                              
dim(rec) = c(2, length(rec)/2);   
x0 = which(rec[2,] == "CGE"); rec[2,x0] = "STR";
x0 = which(rec[2,] == "DTH"); rec[2,x0] = "MD"; 
x0 = which(rec[2,] == "LGE"); rec[2,x0] ="STR";
x0 = which(rec[2,] == "MGE"); rec[2,x0] = "STR";
x0 = which(rec[2,] == "OC"); rec[2,x0] = "V1C";
x0 = which(rec[2,] == "PC");rec[2,x0] = "IPC";
x0 = which(rec[2,] == "URL"); rec[2,x0] = "CBC";
x0 = which(rec[2,] == "TC"); rec[2,x0] = "ITC";
sampList2 = paste(rec[1,],rec[2,],sep=".");
sampList = c("Geneid",sampList2);  


##-----sequencing depth info
mydata = read.table("result/annot/brainSpan.chrom.read.xls",header=T,sep="\t");
sampSeqdepth = as.numeric(mydata[,3])*(1-mydata[,ncol(mydata)]/100);
sampSeqdepth = round(sampSeqdepth);
sampList_tmp = paste(as.character(mydata[,1]),as.character(mydata[,2]),sep=".");
myindex = grep("MSC", sampList_tmp,fixed=T);
addS1C = gsub("MSC","S1C",sampList_tmp[myindex]);
sampSeqdepth = c(sampSeqdepth, sampSeqdepth[myindex]);
sampList_tmp = c(sampList_tmp, addS1C);
sampList_tmp = gsub("MSC","M1C", sampList_tmp);
names(sampSeqdepth) = sampList_tmp;
##-----re-name region name
sampList2 = sampList_tmp;
rec = unlist(strsplit(sampList2,split=".",fixed=T))                                              
dim(rec) = c(2, length(rec)/2);   
x0 = which(rec[2,] == "CGE"); rec[2,x0] = "STR";
x0 = which(rec[2,] == "DTH"); rec[2,x0] = "MD"; 
x0 = which(rec[2,] == "LGE"); rec[2,x0] ="STR";
x0 = which(rec[2,] == "MGE"); rec[2,x0] = "STR";
x0 = which(rec[2,] == "OC"); rec[2,x0] = "V1C";
x0 = which(rec[2,] == "PC");rec[2,x0] = "IPC";
x0 = which(rec[2,] == "URL"); rec[2,x0] = "CBC";
x0 = which(rec[2,] == "TC"); rec[2,x0] = "ITC";
sampList2 = paste(rec[1,],rec[2,],sep=".");
names(sampSeqdepth) = sampList2;



##----brain information
brainInfo = read.csv("data/brainseq.info.txt",header=T,sep="\t");

##------regions
regionList = c("OFC","DFC","VFC","MFC","M1C","S1C","IPC","A1C","STC","ITC","V1C","HIP","AMY","STR","MD","CBC");
stageList = c(
			"HSB112-HSB148",
			"HSB153-HSB150-HSB113-HSB103-HSB149-HSB114",
			"HSB178-HSB154-HSB96-HSB97",
			"HSB98-HSB107-HSB92-HSB159",
			"HSB155-HSB194-HSB121-HSB132-HSB139",
			"HSB131-HSB171-HSB122-HSB143-HSB173",
			"HSB172-HSB118-HSB141-HSB174-HSB175",
			"HSB124-HSB119-HSB105-HSB127",
			"HSB130-HSB136-HSB126-HSB145-HSB123-HSB135"
			);
			
			
##----pair wise comparison
for(i in 1:length(regionList)){
	regionName = regionList[i];
	
	for(j in 1:(length(stageList)-1)){
		for(k in (j+1):length(stageList)){
			stageName1 = stageList[j];
			stageName2 = stageList[k];

			##----get the expression for two groups
			##Group 1
			##----split sliding window
			rec = unlist(strsplit(as.character(stageName1),split="-",fixed=T));
			myindex = c();
			for(kk in 1:length(rec)){
				stageNameTmp = rec[kk];
				myindex = c(myindex,which(as.character(brainInfo[,2]) == stageNameTmp));
			}
			##--
			brainStage = as.character(brainInfo[myindex,2]);
			brainSite = as.character(brainInfo[myindex,9]);
			oneIndex = match(paste(brainStage,regionName,sep="."), sampList);
			oneSite = brainSite[!is.na(oneIndex)];
			oneIndex = oneIndex[!is.na(oneIndex)];
			oneSample = sampList[oneIndex];
		
			##Group 2
			##----split sliding window
			rec = unlist(strsplit(as.character(stageName2),split="-",fixed=T));
			myindex = c();
			for(kk in 1:length(rec)){
				stageNameTmp = rec[kk];
				myindex = c(myindex,which(as.character(brainInfo[,2]) == stageNameTmp));
			}
			##--
			brainStage = as.character(brainInfo[myindex,2]);
			brainSite = as.character(brainInfo[myindex,9]);
			twoIndex = match(paste(brainStage,regionName,sep="."), sampList);
			twoSite = brainSite[!is.na(twoIndex)];
			twoIndex = twoIndex[!is.na(twoIndex)];
			twoSample = sampList[twoIndex];
			
			##------skip comparison having duplcated samples
			if(length(intersect(oneSample, twoSample)) > 0){
				cat("Duplicated samples in compared pairs",regionName,stageName1,stageName2,length(oneIndex), length(twoIndex),"\n",sep="\t");
				next;
			}
			
			##-----screen for debug
			if(length(oneIndex) < 2 || length(twoIndex) < 2){
				cat("No replicates",regionName,stageName1,stageName2,length(oneIndex), length(twoIndex),"\n",sep="\t");
				next;
			}


			##----print screen
			cat(regionName,stageName1,stageName2,length(oneIndex), length(twoIndex),"\n",sep="\t");
			
			##-----get gene expression
			geneCount1 = geneCount[,oneIndex];
			geneCount2 = geneCount[,twoIndex];
			geneRPKM1 = geneRPKM[,oneIndex];
			geneRPKM2 = geneRPKM[,twoIndex];
			geneCountTable = cbind(geneCount1,geneCount2);
			rownames(geneCountTable) = geneList;
			colnames(geneCountTable) = c(oneSample,twoSample);
			
			##------mean of expression
			RPKM1=apply(geneRPKM1,1,mean);
			RPKM2=apply(geneRPKM2,1,mean);
			count1=apply(geneCount1,1,mean);
			count2=apply(geneCount2,1,mean);
			
			##-----meta data
			sampDesign = data.frame(row.names = colnames(geneCountTable),
									condition = c(rep("W",ncol(geneCount1)),rep("K",ncol(geneCount2))),
									type = c(oneSite, twoSite));
			
	
	
			##-----input from count matrix
			dds = DESeqDataSetFromMatrix(countData = geneCountTable,
										colData = sampDesign,
										design = ~ condition);
			dds$condition = factor(dds$condition,levels=c("W","K"));
			
			##--------------------------choose one normalization method ---------------------------------------##
			##---1-----cqn normlization
			dds1= dds;
			seqDepth = sampSeqdepth[match(colnames(geneCountTable), names(sampSeqdepth))];
			cqnObject = cqn2(geneCountTable,lengths=as.numeric(geneInfo[,2]), x= as.numeric(geneInfo[,3]),
							sizeFactors=seqDepth,lengthMethod="smooth",sqn=T);
			cqnOffset = cqnObject$glm.offset;
			cqnNormFactors = exp(cqnOffset);
			cqnNormFactors = cqnNormFactors/mean(cqnNormFactors);
			normalizationFactors(dds1) = cqnNormFactors;
			dds1 = estimateDispersions(dds1);
			
			##---2------only sequencing depth normalization
			#dds2 = dds;
			#dds2 <- estimateSizeFactors(dds2);
			#dds2 = estimateDispersions(dds2);
			
			##!!!!!!!!!choose Normalization method!!!!!!!!!!!!!!!!!!
			ddsMF = dds1;
			
			##-----DESeq2 do comparison for multiple factor
			siteNumber = length(unique(c(oneSite,twoSite)));
			if(siteNumber > 1){
				design(ddsMF) = formula(~ type + condition);
				ddsMF = nbinomWaldTest(ddsMF);
				resMF = results(ddsMF);
			}
			else{
				ddsMF = nbinomWaldTest(ddsMF);
				resMF = results(ddsMF);
			}

			res = cbind(rownames(resMF), 
						RPKM1, RPKM2,
						count1, count2,
						2**resMF[,2],
						resMF[,2],
						resMF[,5],
						resMF[,6]
						)
			
			##---add gene type and assign change type
			myindex=match(as.character(res[,1]),as.character(geneAnnot[,1]));
			geneType=as.character(geneAnnot[myindex,2]);
			res=cbind(res[,1],geneType,res[,2:ncol(res)]);
			colnames(res)=c("Geneid","geneType","baseMeanA(RPKM)","baseMeanB(RPKM)","baseMeanA(count)","baseMeanB(count)",
							"foldChange","log2FoldChange", "pval","padj");			
			##-----output
			outFile = paste("result/DEX/dexDESeq2_developWindow/","W",j,"-vs-", "W",k, ".", regionName,".all.DEX.xls",sep="");
			write.table(res,file=outFile,quote=F,row.names=F, col.names=T,sep="\t");
		}
	}
}	



#########################################
##AIM: get spatial differentially expressed genes
library("DESeq2");
library("cqn");
source("cqn2.R");

##-----------get gene annotation
geneAnnot = read.table("data/gencode.v21.gene.annot.txt");

##---gene info
geneInfo = read.table("data/gencode.v21.wholeGene.geneComposite.geneGCcontent.txt",header=T,sep="\t");
 
##---get gene RPKM
geneRPKM = read.table("result/genexp/brainSpan.hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.RPKM.txt",header=T,sep="\t");

##---get gene count
geneCount = read.table("result/genexp/brainSpan.hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.count.txt",header=T,sep="\t");
sampList = colnames(geneCount);
geneList = as.character(geneCount[,1]);

##----copy some samples
myindex = grep("MSC", sampList,fixed=T);
addS1C = gsub("MSC","S1C",sampList[myindex]);
geneRPKM = cbind(geneRPKM, geneRPKM[,myindex]);
geneCount = cbind(geneCount, geneCount[,myindex]);
geneList = as.character(geneCount[,1]);  ###
sampList = c(sampList,addS1C) ;          ###
sampList = gsub("MSC","M1C", sampList);
colnames(geneRPKM) = sampList;
colnames(geneCount) = sampList;

####----------------seq
##-----re-name region name
sampList2 = sampList[-1];
rec = unlist(strsplit(sampList2,split=".",fixed=T))                                              
dim(rec) = c(2, length(rec)/2);   
x0 = which(rec[2,] == "CGE"); rec[2,x0] = "STR";
x0 = which(rec[2,] == "DTH"); rec[2,x0] = "MD"; 
x0 = which(rec[2,] == "LGE"); rec[2,x0] ="STR";
x0 = which(rec[2,] == "MGE"); rec[2,x0] = "STR";
x0 = which(rec[2,] == "OC"); rec[2,x0] = "V1C";
x0 = which(rec[2,] == "PC");rec[2,x0] = "IPC";
x0 = which(rec[2,] == "URL"); rec[2,x0] = "CBC";
x0 = which(rec[2,] == "TC"); rec[2,x0] = "ITC";
sampList2 = paste(rec[1,],rec[2,],sep=".");
sampList = c("Geneid",sampList2);  



##-----sequencing depth info
mydata = read.table("result/annot/brainSpan.chrom.read.xls",header=T,sep="\t")    ;
sampSeqdepth = as.numeric(mydata[,3])*(1-mydata[,ncol(mydata)]/100);
sampSeqdepth = round(sampSeqdepth);
sampList_tmp = paste(as.character(mydata[,1]),as.character(mydata[,2]),sep=".");
myindex = grep("MSC", sampList_tmp,fixed=T);
addS1C = gsub("MSC","S1C",sampList_tmp[myindex]);
sampSeqdepth = c(sampSeqdepth, sampSeqdepth[myindex]);
sampList_tmp = c(sampList_tmp, addS1C);
sampList_tmp = gsub("MSC","M1C", sampList_tmp);
names(sampSeqdepth) = sampList_tmp
##-----re-name region name
sampList2 = sampList_tmp;
rec = unlist(strsplit(sampList2,split=".",fixed=T))                                              
dim(rec) = c(2, length(rec)/2);   
x0 = which(rec[2,] == "CGE"); rec[2,x0] = "STR";
x0 = which(rec[2,] == "DTH"); rec[2,x0] = "MD"; 
x0 = which(rec[2,] == "LGE"); rec[2,x0] ="STR";
x0 = which(rec[2,] == "MGE"); rec[2,x0] = "STR";
x0 = which(rec[2,] == "OC"); rec[2,x0] = "V1C";
x0 = which(rec[2,] == "PC");rec[2,x0] = "IPC";
x0 = which(rec[2,] == "URL"); rec[2,x0] = "CBC";
x0 = which(rec[2,] == "TC"); rec[2,x0] = "ITC";
sampList2 = paste(rec[1,],rec[2,],sep=".");
names(sampSeqdepth) = sampList2;


##----brain information
brainInfo = read.csv("data/brainseq.info.txt",header=T,sep="\t");

##------regions
regionList = c("OFC","DFC","VFC","MFC","M1C","S1C","IPC","A1C","STC","ITC","V1C","HIP","AMY","STR","MD","CBC");
stageList = c(
			"HSB112-HSB148",
			"HSB153-HSB150-HSB113-HSB103-HSB149-HSB114",
			"HSB178-HSB154-HSB96-HSB97",
			"HSB98-HSB107-HSB92-HSB159",
			"HSB155-HSB194-HSB121-HSB132-HSB139",
			"HSB131-HSB171-HSB122-HSB143-HSB173",
			"HSB172-HSB118-HSB141-HSB174-HSB175",
			"HSB124-HSB119-HSB105-HSB127",
			"HSB130-HSB136-HSB126-HSB145-HSB123-HSB135"
			);

##----pair wise comparison
for(i in 1:length(stageList)){
	stageName = stageList[i];
	
	##----split sliding window
	rec = unlist(strsplit(as.character(stageName),split="-",fixed=T));
	myindex = c();
	for(kk in 1:length(rec)){
		stageName2 = rec[kk];
		myindex = c(myindex,which(as.character(brainInfo[,2]) == stageName2));
	}	
	
	##----get brain list and site list	
	brainStage = as.character(brainInfo[myindex,2]);
	brainSite = as.character(brainInfo[myindex,9]);
	
	for(j in 1:(length(regionList)-1)){
		for(k in (j+1):length(regionList)){
			regionName1 = regionList[j];
			regionName2 = regionList[k];
		
			##----get the expression for two groups
			##Group 1
			oneIndex = match(paste(brainStage,regionName1,sep="."), sampList);
			oneSite = brainSite[!is.na(oneIndex)];
			oneIndex = oneIndex[!is.na(oneIndex)];
			oneSample = sampList[oneIndex];
			
			##Group 2
			twoIndex = match(paste(brainStage,regionName2,sep="."), sampList);
			twoSite = brainSite[!is.na(twoIndex)];
			twoIndex = twoIndex[!is.na(twoIndex)];
			twoSample = sampList[twoIndex];
	
			##-----screen for debug
			if(length(oneIndex) < 2  || length(twoIndex) < 2 ) {
				cat("No replicates",stageName,regionName1,regionName2,length(oneIndex), length(twoIndex),"\n",sep="\t");
				next;
			}
			
			##-----print screen
			cat(stageName,regionName1,regionName2,length(oneIndex), length(twoIndex),"\n",sep="\t");
			
			
			##-----get gene expression
			geneCount1 = geneCount[,oneIndex];
			geneCount2 = geneCount[,twoIndex];
			geneRPKM1 = geneRPKM[,oneIndex];
			geneRPKM2 = geneRPKM[,twoIndex];
			geneCountTable = cbind(geneCount1,geneCount2);
			rownames(geneCountTable) = geneList;
			colnames(geneCountTable) = c(oneSample,twoSample);
			
			##------mean of expression
			RPKM1=apply(geneRPKM1,1,mean);
			RPKM2=apply(geneRPKM2,1,mean);
			count1=apply(geneCount1,1,mean);
			count2=apply(geneCount2,1,mean);

			##-----meta data
			sampDesign = data.frame(row.names = colnames(geneCountTable),
									condition = c(rep("W",ncol(geneCount1)),rep("K",ncol(geneCount2))),
									type = c(oneSite, twoSite));	
	
			##-----input from count matrix
			dds = DESeqDataSetFromMatrix(countData = geneCountTable,
										colData = sampDesign,
										design = ~ condition);
			dds$condition = factor(dds$condition,levels=c("W","K"));
			
			##--------------------------choose one normalization method ---------------------------------------##
			##---1-----cqn normlization
			dds1= dds;
			seqDepth = sampSeqdepth[match(colnames(geneCountTable), names(sampSeqdepth))];
			cqnObject = cqn2(geneCountTable,lengths=as.numeric(geneInfo[,2]), x= as.numeric(geneInfo[,3]),
							sizeFactors=seqDepth,lengthMethod="smooth",sqn=T);
			cqnOffset = cqnObject$glm.offset;
			cqnNormFactors = exp(cqnOffset);
			cqnNormFactors = cqnNormFactors/mean(cqnNormFactors);
			normalizationFactors(dds1) = cqnNormFactors;
			dds1 = estimateDispersions(dds1);
			
			##---2------only sequencing depth normalization
			#dds2 = dds;
			#dds2 <- estimateSizeFactors(dds2);
			#dds2 = estimateDispersions(dds2);
			
			##!!!!!!!!!choose Normalization method!!!!!!!!!!!!!!!!!!
			ddsMF = dds1;
			
			##-----DESeq2 do comparison for multiple factor
			siteNumber = length(unique(c(oneSite,twoSite)));
			if(siteNumber > 1){
				design(ddsMF) = formula(~ type + condition);
				ddsMF = nbinomWaldTest(ddsMF);
				resMF = results(ddsMF);
			}
			else{
				ddsMF = nbinomWaldTest(ddsMF);
				resMF = results(ddsMF);
			}

			res = cbind(rownames(resMF), 
						RPKM1, RPKM2,
						count1, count2,
						2**resMF[,2],
						resMF[,2],
						resMF[,5],
						resMF[,6]
						)
			
			##---add gene type and assign change type
			myindex=match(as.character(res[,1]),as.character(geneAnnot[,1]));
			geneType=as.character(geneAnnot[myindex,2]);
			res=cbind(res[,1],geneType,res[,2:ncol(res)]);
			colnames(res)=c("Geneid","geneType","baseMeanA(RPKM)","baseMeanB(RPKM)","baseMeanA(count)","baseMeanB(count)",
							"foldChange","log2FoldChange", "pval","padj");
			
			##---output
			outFile = paste("result/DEX/dexDESeq2_regionWindow/",regionName1,"-vs-", regionName2, ".W",i,".all.DEX.xls",sep="");
			write.table(res,file=outFile,quote=F,row.names=F, col.names=T,sep="\t");
		}
	}
}	







##AIM: conduct weighted gene co-expression network analysis


###################-------------
##filter gene expression for WGCNA
##----sample meta data
metaData = read.table("data/brainSpan.metaData.txt",skip=3,header=F,sep="\t");
brainList = unique(as.character(metaData[,3]));

dfcList=c("DFC","FC");
vfcList=c("VFC","FC");
mfcList=c("MFC","FC");
ofcList=c("OFC","FC");
m1cList=c("M1C","FC","MS","MSC","M1CS1C");
s1cList=c("S1C","PC","MS","MSC","M1CS1C");
ipcList=c("IPC","PC");
a1cList=c("A1C","TC");
stcList=c("STC","TC");
itcList=c("ITC","TC");
v1cList=c("V1C","OC");
ncxList=c("DFC","VFC","MFC","OFC","MS","MSC","M1C","S1C","M1CS1C","IPC","A1C","STC","ITC","V1C","FC","PC","TC","OC");
hipList=c("HIP");
amyList=c("AMY");
strList=c("STR","VF","MGE","LGE","CGE","BF");
mdList=c("THM","DIE","DTH","MD");
cbcList=c("URL","CBC");


##----------get gene expression
mydata = read.csv("result/genexp/brainSpan.hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.RPKM.normalized.CQNCombat.txt",header=T,sep="\t");
genexp = mydata[, -1];
rownames(genexp) = as.character(mydata[,1]);

##---------exclude low expressed genes
hiexpindex = which(apply(genexp,1,function (x) {y=length(which(x >=1))}) >= 5);
genexp = genexp[hiexpindex,];
geneList = rownames(genexp);
sampList = colnames(genexp);
##get regionList
rec = unlist(strsplit(sampList,split=".",fixed=T));
dim(rec)=c(2,length(rec)/2);
regionList = rec[2,];


##----re-organize gene expression
recGenexp = geneList;
recSample = "Geneid";
for(i in 1:length(brainList)){
	brainName = brainList[i];
	brainIndex = grep(brainName,sampList,fixed=T);
	cat(brainName,"\n");
	
	##OFC
	regionIndex = c();
	for(k in 1:length(ofcList)){
		y = which(regionList == ofcList[k]);
		regionIndex = c(regionIndex, y);
	}
	x = unique(intersect(brainIndex, regionIndex));	
	if(length(x) > 0){
		recSample = c(recSample, paste(brainName,"OFC",sep="."));
		if(length(x) > 1){
			cat("Error found OFC!\n");
			avexp = apply(genexp[,x],1,mean);
			recGenexp = cbind(recGenexp,avexp);
		}
		else{
			recGenexp = cbind(recGenexp,genexp[,x]);
		}
	}

	##DFC
	regionIndex = c();
	for(k in 1:length(dfcList)){
		y = which(regionList == dfcList[k]);
		regionIndex = c(regionIndex, y);
	}
	x = unique(intersect(brainIndex, regionIndex));	
	if(length(x) > 0){
		recSample = c(recSample, paste(brainName,"DFC",sep="."));
		if(length(x) > 1){
			cat("Error found DFC!\n");
			avexp = apply(genexp[,x],1,mean);
			recGenexp = cbind(recGenexp,avexp);
		}
		else{
			recGenexp = cbind(recGenexp,genexp[,x]);
		}
	}	

	##VFC
	regionIndex = c();
	for(k in 1:length(vfcList)){
		y = which(regionList == vfcList[k]);
		regionIndex = c(regionIndex, y);
	}
	x = unique(intersect(brainIndex, regionIndex));	
	if(length(x) > 0){
		recSample = c(recSample, paste(brainName,"VFC",sep="."));
		if(length(x) > 1){
			cat("Error found VFC!\n");
			avexp = apply(genexp[,x],1,mean);
			recGenexp = cbind(recGenexp,avexp);
		}
		else{
			recGenexp = cbind(recGenexp,genexp[,x]);
		}
	}

	##MFC
	regionIndex = c();
	for(k in 1:length(mfcList)){
		y = which(regionList == mfcList[k]);
		regionIndex = c(regionIndex, y);
	}
	x = unique(intersect(brainIndex, regionIndex));	
	if(length(x) > 0){
		recSample = c(recSample, paste(brainName,"MFC",sep="."));
		if(length(x) > 1){
			cat("Error found MFC!\n");
			avexp = apply(genexp[,x],1,mean);
			recGenexp = cbind(recGenexp,avexp);
		}
		else{
			recGenexp = cbind(recGenexp,genexp[,x]);
		}
	}

	##M1C
	regionIndex = c();
	for(k in 1:length(m1cList)){
		y = which(regionList == m1cList[k]);
		regionIndex = c(regionIndex, y);
	}
	x = unique(intersect(brainIndex, regionIndex));	
	if(length(x) > 0){
		recSample = c(recSample, paste(brainName,"M1C",sep="."));
		if(length(x) > 1){
			cat("Error found M1C!\n");
			avexp = apply(genexp[,x],1,mean);
			recGenexp = cbind(recGenexp,avexp);
		}
		else{
			recGenexp = cbind(recGenexp,genexp[,x]);
		}
	}
	
	##S1C
	regionIndex = c();
	for(k in 1:length(s1cList)){
		y = which(regionList == s1cList[k]);
		regionIndex = c(regionIndex, y);
	}
	x = unique(intersect(brainIndex, regionIndex));	
	if(length(x) > 0){
		recSample = c(recSample, paste(brainName,"S1C",sep="."));
		if(length(x) > 1){
			cat("Error found S1C!\n");
			avexp = apply(genexp[,x],1,mean);
			recGenexp = cbind(recGenexp,avexp);
		}
		else{
			recGenexp = cbind(recGenexp,genexp[,x]);
		}
	}
	
	##IPC
	regionIndex = c();
	for(k in 1:length(ipcList)){
		y = which(regionList == ipcList[k]);
		regionIndex = c(regionIndex, y);
	}
	x = unique(intersect(brainIndex, regionIndex));	
	if(length(x) > 0){
		recSample = c(recSample, paste(brainName,"IPC",sep="."));
		if(length(x) > 1){
			cat("Error found IPC!\n");
			avexp = apply(genexp[,x],1,mean);
			recGenexp = cbind(recGenexp,avexp);
		}
		else{
			recGenexp = cbind(recGenexp,genexp[,x]);
		}
	}
	
	##A1C
	regionIndex = c();
	for(k in 1:length(a1cList)){
		y = which(regionList == a1cList[k]);
		regionIndex = c(regionIndex, y);
	}
	x = unique(intersect(brainIndex, regionIndex));	
	if(length(x) > 0){
		recSample = c(recSample, paste(brainName,"A1C",sep="."));
		if(length(x) > 1){
			cat("Error found A1C!\n");
			avexp = apply(genexp[,x],1,mean);
			recGenexp = cbind(recGenexp,avexp);
		}
		else{
			recGenexp = cbind(recGenexp,genexp[,x]);
		}
	}
	
	##STC
	regionIndex = c();
	for(k in 1:length(stcList)){
		y = which(regionList == stcList[k]);
		regionIndex = c(regionIndex, y);
	}
	x = unique(intersect(brainIndex, regionIndex));	
	if(length(x) > 0){
		recSample = c(recSample, paste(brainName,"STC",sep="."));
		if(length(x) > 1){
			cat("Error found STC!\n");
			avexp = apply(genexp[,x],1,mean);
			recGenexp = cbind(recGenexp,avexp);
		}
		else{
			recGenexp = cbind(recGenexp,genexp[,x]);
		}
	}
	
	##ITC
	regionIndex = c();
	for(k in 1:length(itcList)){
		y = which(regionList == itcList[k]);
		regionIndex = c(regionIndex, y);
	}
	x = unique(intersect(brainIndex, regionIndex));	
	if(length(x) > 0){
		recSample = c(recSample, paste(brainName,"ITC",sep="."));
		if(length(x) > 1){
			cat("Error found ITC!\n");
			avexp = apply(genexp[,x],1,mean);
			recGenexp = cbind(recGenexp,avexp);
		}
		else{
			recGenexp = cbind(recGenexp,genexp[,x]);
		}
	}
	
	##V1C
	regionIndex = c();
	for(k in 1:length(v1cList)){
		y = which(regionList == v1cList[k]);
		regionIndex = c(regionIndex, y);
	}
	x = unique(intersect(brainIndex, regionIndex));	
	if(length(x) > 0){
		recSample = c(recSample, paste(brainName,"V1C",sep="."));
		if(length(x) > 1){
			cat("Error found V1C!\n");
			avexp = apply(genexp[,x],1,mean);
			recGenexp = cbind(recGenexp,avexp);
		}
		else{
			recGenexp = cbind(recGenexp,genexp[,x]);
		}
	}
	
	##NCX
	regionIndex = c();
	for(k in 1:length(ncxList)){
		y = which(regionList == ncxList[k]);
		regionIndex = c(regionIndex, y);
	}
	x = unique(intersect(brainIndex, regionIndex));	
	if(length(x) > 0){
		recSample = c(recSample, paste(brainName,"NCX",sep="."));
		if(length(x) < 2){
			cat("Error found NCX!\n");
		}
		else{
			avexp = apply(genexp[,x],1,mean);
			recGenexp = cbind(recGenexp,avexp);
		}
	}
	
	##HIP
	regionIndex = c();
	for(k in 1:length(hipList)){
		y = which(regionList == hipList[k]);
		regionIndex = c(regionIndex, y);
	}
	x = unique(intersect(brainIndex, regionIndex));	
	if(length(x) > 0){
		recSample = c(recSample, paste(brainName,"HIP",sep="."));
		if(length(x) > 1){
			cat("Error found HIP!\n");
			avexp = apply(genexp[,x],1,mean);
			recGenexp = cbind(recGenexp,avexp);
		}
		else{
			recGenexp = cbind(recGenexp,genexp[,x]);
		}
	}
	
	##AMY
	regionIndex = c();
	for(k in 1:length(amyList)){
		y = which(regionList == amyList[k]);
		regionIndex = c(regionIndex, y);
	}
	x = unique(intersect(brainIndex, regionIndex));	
	if(length(x) > 0){
		recSample = c(recSample, paste(brainName,"AMY",sep="."));
		if(length(x) > 1){
			cat("Error found AMY!\n");
			avexp = apply(genexp[,x],1,mean);
			recGenexp = cbind(recGenexp,avexp);
		}
		else{
			recGenexp = cbind(recGenexp,genexp[,x]);
		}
	}
	
	##STR
	regionIndex = c();
	for(k in 1:length(strList)){
		y = which(regionList == strList[k]);
		regionIndex = c(regionIndex, y);
	}
	x = unique(intersect(brainIndex, regionIndex));	
	if(length(x) > 0){
		recSample = c(recSample, paste(brainName,"STR",sep="."));
		if(length(x) > 1){
			cat("Error found STR!\n");
			avexp = apply(genexp[,x],1,mean);
			recGenexp = cbind(recGenexp,avexp);
		}
		else{
			recGenexp = cbind(recGenexp,genexp[,x]);
		}
	}
	
	##MD
	regionIndex = c();
	for(k in 1:length(mdList)){
		y = which(regionList == mdList[k]);
		regionIndex = c(regionIndex, y);
	}
	x = unique(intersect(brainIndex, regionIndex));	
	if(length(x) > 0){
		recSample = c(recSample, paste(brainName,"MD",sep="."));
		if(length(x) > 1){
			cat("Error found MD!\n");
			avexp = apply(genexp[,x],1,mean);
			recGenexp = cbind(recGenexp,avexp);
		}
		else{
			recGenexp = cbind(recGenexp,genexp[,x]);
		}
	}
	
	##CBC
	regionIndex = c();
	for(k in 1:length(cbcList)){
		y = which(regionList == cbcList[k]);
		regionIndex = c(regionIndex, y);
	}
	x = unique(intersect(brainIndex, regionIndex));	
	if(length(x) > 0){
		recSample = c(recSample, paste(brainName,"CBC",sep="."));
		if(length(x) > 1){
			cat("Error found CBC!\n");
			avexp = apply(genexp[,x],1,mean);
			recGenexp = cbind(recGenexp,avexp);
		}
		else{
			recGenexp = cbind(recGenexp,genexp[,x]);
		}
	}	
}	
colnames(recGenexp) = recSample;
rownames(recGenexp) = geneList;


##----saving data
genexp = recGenexp;
save(genexp, file="result/WGCNA/brainSpan.genexp.WGCNA.filter.RData");	
		

		
#######################-----------------------------		
##########build threshold

library("WGCNA");

##-----
pdf("result/WGCNA/brainSpan.genexp.brainWGCNA.threshold.pdf",12,7);
par(omi=c(0.1,0.1,0.1,0.1));


##-----isoform expression
load("result/WGCNA/brainSpan.genexp.WGCNA.filter.RData");
recGenexp = genexp[,-1];
sampList = colnames(recGenexp);
geneList = as.character(genexp[,1]);

##-----only using 6 regions
x = c(grep("NCX", sampList),
	 grep("HIP", sampList),
	 grep("AMY", sampList),
	 grep("STR", sampList),
	 grep("MD", sampList),
	 grep("CBC", sampList)
	);
recGenexp = recGenexp[,x];
sampList2 = sampList[x];

###############filtering
##low variable genes
backup = recGenexp;
recGenexp = as.numeric(recGenexp);
dim(recGenexp) = dim(backup);

myAve = apply(recGenexp,1, mean);
myStd = apply(recGenexp,1, sd);
myCv = myStd/myAve;
plot(density(myCv),col="deepskyblue",lwd=2,xlab="Coefficient variation",ylab="Density", 
		main="Gene expression");
abline(v=1, lty=2);
abline(v=2, lty=2);
abline(v=3, lty=2);

##----log2 transformed
recGenexp = log2(recGenexp+1);

##---saving filtered data
genexp = recGenexp;
colnames(genexp) = sampList2;
rownames(genexp) = geneList;
save(genexp,file="result/WGCNA/brainSpan.genexp.brainWGCNA.filterAgain.RData");

############-----------simulation for soft threshold
recGenexp = t(recGenexp);

powers = c(c(1:30), seq(from = 31, to=50, by=2));
sft = pickSoftThreshold(recGenexp, powerVector = powers, verbose = 3)


# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence ALL"),ylim=c(-1,1));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.60,col="blue",lty=2)
abline(h=0.70,col="blue",lty=2)
abline(h=0.80,col="blue",lty=2)
abline(h=0.90,col="blue",lty=2)
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity ALL"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,col="red")

dev.off();


###############################---------------------------
##buid modules

library("WGCNA");

##------set figure parameters
pdf("result/WGCNA/brainSpan.genexp.brainWGCNA.module-signed.pdf",7,7);
par(omi=c(0.1,0.1,0.1,0.1));

#############--------------embryonic cluter expression
load("result/WGCNA/brainSpan.genexp.brainWGCNA.filterAgain.RData");
recGenexp = genexp;

recGenexp = t(recGenexp);

##----------get modules
geneNet = blockwiseModules(recGenexp, power=25,minModuleSize=10,mergeCutHeight=0.20,networkType = "signed",TOMType = "signed",numericLabels=T,pamRespectsDendro=F,saveTOMs=F,verbose=0);


##---------plots
mergedColors = labels2colors(geneNet$colors);
plotDendroAndColors(geneNet$dendrograms[[1]],mergedColors[geneNet$blockGenes[[1]]],"Module colors",rowText=paste("M",geneNet$colors,sep=""), addTextGuide=T, dendroLabels = F, hang = 0.03,addGuide = TRUE,guideHang = 0.05,marAll = c(1, 5, 3,1));
##---save data
save(genexp,geneNet, file="result/WGCNA/brainSpan.genexp.brainWGCNA.gene-net-signed.RData");


dev.off()


###########################---------------
##check modules by correlation analysis
library("WGCNA");


##-----data
inFile = "result/WGCNA/brainSpan.genexp.brainWGCNA.gene-net-signed.RData";
outFile="result/WGCNA/brainSpan.genexp.brainWGCNA.geneModule.reassigned.xls";
outFile2="result/WGCNA/brainSpan.genexp.brainWGCNA.geneModule.reassigned.geneid.xls";
outFile3="result/WGCNA/brainSpan.genexp.brainWGCNA.geneModule.reassigned.geneSymbol.xls";

if(file.exists(outFile)) file.remove(outFile);
if(file.exists(outFile2)) file.remove(outFile2);
if(file.exists(outFile3)) file.remove(outFile3);


##loading data
load(inFile);
eigenexp = geneNet$MEs;
moduleList = colnames(eigenexp);
genexp = t(genexp);
geneList = colnames(genexp);
colorList = geneNet$colors;

res=c();
for(i in 0:max(colorList)){
	myindex = which(colorList == i);
	clusterGene = as.character(geneList[myindex]);
	for(j in 1:length(clusterGene)){
		geneName = as.character(clusterGene[j]);
		myindex2 = which(geneList == geneName);
		oneGenexp = as.numeric(genexp[,myindex2]);
		oneGenecor = cor(oneGenexp, eigenexp,method="p");
		
		##---calculate correlation
		moduleName.old = paste("ME",i,sep="");
		myindex3 = which(moduleList == moduleName.old);
		
		moduleCor.old = oneGenecor[myindex3];
		moduleCor.new = max(oneGenecor);
		
		##----reassign module
		output = c(geneName, moduleName.old);
		if((moduleCor.new -moduleCor.old) >0.3){
			maxIndex = which(oneGenecor == max(oneGenecor));
			output = c(geneName, moduleList[maxIndex]);
		}	
		
		##recording	
		if(length(res) == 0){
			res = output;
		}
		else{
			res = rbind(res,output);
		}
	}
}	



###-----output
for(i in 0:max(colorList)){
	moduleName = paste("ME",i,sep="");
	myindex = which(as.character(res[,2]) == moduleName);
	clusterGene = as.character(res[myindex,1]);
	moduleSize = length(clusterGene);
	cat(moduleName,moduleSize,"\n",sep="\t");
	cat(moduleName,clusterGene,"\n",file=outFile,sep="\t",append=T);
	
	##-----split gene id and symbol
	tmpRec = unlist(strsplit(clusterGene,split="|",fixed=T));
	dim(tmpRec) = c(2, length(tmpRec)/2);
	cat(moduleName,tmpRec[1,],"\n",file=outFile2,sep="\t",append=T);
	cat(moduleName,tmpRec[2,],"\n",file=outFile3,sep="\t",append=T);	
}




###################################------------------------
##analyze eigengenes of module

library("WGCNA");
library("gplots");
library("som");



##------lowess function for missing data
lowess.na <- function(x, y = NULL, f = 2/3,...) {  #do lowess with missing data
	x1 <- subset(x,(!is.na(x)) &(!is.na(y)))
	y1 <- subset(y, (!is.na(x)) &(!is.na(y)))
	lowess.na <- lowess(x1,y1,f, ...)
}

##-----color information
mycolor=c("blue","cyan","orange","black","seagreen","red");
mysamp=c("NCX","HIP","AMY","STR","MD","CBC");



##------ gene expression data
load("result/WGCNA/brainSpan.genexp.brainWGCNA.gene-net-signed.RData");
genexp = t(genexp);
sampList = rownames(genexp);
geneList = colnames(genexp);

##-------gene meta data
metaData = read.table("data/brainSpan.metaData.txt",skip=3,header=F,sep="\t");
recAge = c();
recSex = c();
recPeriod = c();
for(i in 1:length(sampList)){
	rec = unlist(strsplit(as.character(sampList[i]),split=".",fixed=T));
	x = which(as.character(metaData[,3]) == rec[1]);
	recAge = c(recAge,as.numeric(metaData[x[1],7]));
	recSex = c(recSex,as.character(metaData[x[1],8]));
	recPeriod = c(recPeriod,as.numeric(metaData[x[1],2]));
}



###########################----------protein-coding genes
pdf("result/WGCNA/brainSpan.genexp.brainWGCNA.eigengenes.v2.pdf");
par(omi=c(0.1,0.1,0.1,0.1));
#layout(matrix(1:48,ncol=6,byrow=T));


##-----calculte eigengene
mydata = readLines("result/WGCNA/brainSpan.genexp.brainWGCNA.geneModule.reassigned.xls");
res = c();
resEigen = c();

for(i in 1:length(mydata)){
	rec = unlist(strsplit(as.character(mydata[i]),split="\t",fixed=T));
	moduleName = rec[1];
	clustGene = rec[2:length(rec)];
	moduleSize = length(clustGene);

	if(moduleName != "ME0"){
		cat(moduleName, "\n");
		
		##----get expression
		myindex = match(clustGene,geneList);
		moduleGenexp = genexp[,myindex];
	
		##----calculate eigengene
		myEigen = moduleEigengenes(moduleGenexp,colors=rep("red",ncol(moduleGenexp)));
		geneigen = unlist(myEigen$eigengenes);
		
		##-----recording eigengenes
		if(length(resEigen) == 0){
			resEigen = geneigen;
		}
		else{
			resEigen = rbind(resEigen,geneigen);
		}	

		##---plot figure
		tmpRes=c();
		for(j in 1:length(mysamp)){
			myindex2 = grep(mysamp[j],sampList,fixed=T);
			xdata = log10(recAge[myindex2]);
			ydata = geneigen[myindex2];
			
			##---sex info
			sexdata=recSex[myindex2];
			maleIndex=which(sexdata == "M");
			femaleIndex=which(sexdata == "F");
			
			##---period info
			periodata = recPeriod[myindex2];
			
			##--set axes
			minXdata=min(log10(as.numeric(recAge)),na.rm=T);
			maxXdata=max(log10(as.numeric(recAge)),na.rm=T);
			minYdata=min(geneigen);				
			maxYdata=max(geneigen)*1.2;
			xTick=c(50,100,200,500,2000,10000);
			#yTick=seq(0,round(maxYdata));
	
			##---plot frames
			titleName = paste(moduleName,", target gene =",moduleSize,sep="");
			if(j==1){
				plot(c(minXdata,maxXdata),c(minYdata,maxYdata),type="n",main= titleName,xlab="Age (Days)",ylab="Eigengenes",xaxt="n",cex.lab=1.2);
				axis(1,at=log10(xTick),labels=xTick);
			#	axis(2,at=yTick,labels=yTick,las=1);
			}
			
			##--draw figure
			x1 <- subset(xdata,(!is.na(xdata)) &(!is.na(ydata)));
			if(length(x1)>2){
				##----1-----
				#points(xdata,ydata,pch=myshap[j],col=mycolor[j],cex=0.8);
				
				##----2----separate sex
				points(xdata[maleIndex],ydata[maleIndex],pch=19,col=mycolor[j],cex=0.8);
				points(xdata[femaleIndex],ydata[femaleIndex],pch=17,col=mycolor[j],cex=0.8);

				##--lowess
				newdata=lowess.na(xdata,ydata,f=.5);
				lines(newdata$x,newdata$y,col=mycolor[j],lwd=2);
				
				##---for heatmap
				tmpRes2=newdata$y;
				for(k in 2:13){
					y = which(periodata == k);
					if(length(y) == 0){
						tmpRes = c(tmpRes,NA);
					}else if(length(y) == 1){
						tmpRes = c(tmpRes,tmpRes2[y])
					}else{
						tmpRes = c(tmpRes,mean(tmpRes2[y]));
					}
				}
			}
		}
		
		rect(minXdata-1,minYdata-1,maxXdata+1,maxYdata+1,col=rgb(0.11,0.16,0.05,alpha=.02),border=NA);
	
		
		##---option 2-----Using window
		abline(v=log10(266)); ##birth line
		abline(v=log10(c(63.5,91.5,119.5,154.5,388.5,1175.5,4184.5,7206.5)),,lty="dashed");
		##--write window number
		tmp1=c(63,91,119,154,388,1175,4184,7206,10**maxXdata);
		tmp2=c(10**minXdata,63,91,119,154,388,1175,4184,7206);
		textLoc=log10((tmp1+tmp2)/2);
		textLoc[1]=textLoc[1]-0.05;
		textLoc[5]=textLoc[5]-0.05;
		textLoc[6]=textLoc[6]-0.1;
		textLoc[7]=textLoc[7]-0.1;
		textLoc[8]=textLoc[8]-0.05;
		
		mtext(c("Window:",1:9),side=3,line=0.3,at=c(1.4,textLoc));
		
		
		###---------------
		legend("topright",mysamp,pch=19,col=mycolor,ncol=6,bg="white");
		legend("topleft",c("Male","Female"),pch=c(19,17),col="black",ncol=1,bg="white");
		
		##-----recording
		if(length(res) ==0){
			res=tmpRes;
		}
		else{
			res=rbind(res,tmpRes);
		}
	}
}	
dev.off();



##-----saving data
colnames(resEigen) = sampList;
rownames(resEigen) = paste("ME",1:(length(mydata)-1),sep="");
save(resEigen, file="result/WGCNA/brainSpan.genexp.brainWGCNA.eigengenes.RData");


####-------------plot heatmap of eigengenes
pdf("result/WGCNA/brainSpan.genexp.brainWGCNA.eigengenes.heatmap.pdf",12,12);
res=normalize(res);
rownames(res)=paste("M",1:nrow(res),sep="");
colnames(res)=rep(paste("P",2:13,sep=""),6);
pairs.breaks <- seq(min(res), max(res), length.out=101);
mycol <- colorpanel(n=100,low="blue",mid="white",high="red");

heatmap.2(res,breaks=pairs.breaks,col=mycol,main="",Rowv=T,Colv=F,dendrogram="row"
			,na.rm=TRUE,trace="none",na.color="gray",density.info="none"
			,colsep=12*1:ncol(res),rowsep=1:nrow(res),sepcolor="black",sepwidth=c(0.05,0.0001)
			,ColSideColors=rep(mycolor,each=12)
            ,RowSideColors=c(rep("gray",nrow(res)))
            ,hclustfun = function(x) hclust(x,method = 'ave')
			);
dev.off();








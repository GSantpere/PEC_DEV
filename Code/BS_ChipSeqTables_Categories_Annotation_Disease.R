
##### Prepare chip-seq tables
##### Annotate peaks
##### Compare peaks with Deseq2
##### Annotate GWAS peaks

#### Import libraries
library(reshape2)
library(DESeq2)
library(dplyr)
library(Rsubread)
source('BEDtoolsR.R')

##### Prepare master table (only fetal and adult)
# for i in $(ls *K27ac_peaks.broadPeak); do cut -f1-4 $i | sed 's/-H3K27ac_peak_[0-9]*//g' | sed 's/-K27ac_peak_[0-9]*//g' >${i}.bed; done
# cat *-K27ac_peaks.broadPeak.bed | sortBed | mergeBed  -c 4 -o collapse -delim ' ' >../TABLES/ENHANCERS.collapsed.bed

### Produce master table

e<-read.table('ENHANCERS.collapsed.bed', header=F, sep='\t')
gsub('-','.',e$V4) -> e$V4
samples<-c("HSB100.CBC","HSB100.DFC","HSB102.CBC","HSB102.DFC","HSB240.CBC","HSB240.DFC",
           "HSB126.CBC","HSB126.DFC","HSB136.CBC","HSB136.DFC","HSB187.CBC","HSB187.DFC")
matrix(nrow=dim(e)[1], ncol=length(samples)) -> bin
bin[is.na(bin)]<-0
colnames(bin)<-samples
for (i in c(1:dim(e)[1])){	
  mask<-which(samples %in% unlist(strsplit(e[i,]$V4, " ")))
  bin[i,mask]<-1
}
cbind(e[,c(1:3)], bin)->master
paste('PK',c(1:dim(bin)[1]), sep='')->master$PeakID

## Reproducible enhancers
master[-which(rowSums(master[,c(4:15)])==1),]->masterRep
write.table(masterRep[,c(1,2,3,16)], file='MasterRep.bed', row.names=F, col.names=F, quote=F, sep='\t')

g<-read.table('gencode.v21.annotation.gtf', sep='\t')
g[which(g$V3=='transcript'),c(1,4,5,7,9)]->genes
genes$TSS<-genes$V4  
genes$TSS[which(genes$V7=='-')]<-genes$V5[which(genes$V7=='-')]	
genes$TSSstart<-genes$TSS-1001
genes$TSSend<-genes$TSS+1000
sapply(strsplit(as.character(genes$V9), split=";", fixed=T),"[", 1)->names.full
genes$ENSMBL<-gsub('gene_id ','',names.full)
sapply(strsplit(as.character(genes$V9), split=";", fixed=T),"[", 3)->names.full
genes$gene_type<-gsub(' gene_type ','',names.full)
sapply(strsplit(as.character(genes$V9), split=";", fixed=T),"[", 5)->names.full
genes$GSymbol<-gsub(' gene_name ','',names.full)
genes[,c(1,7,8,4,9,10,11)]->Genes.TSS
Genes.TSS[,2]<-genes$TSS-1
Genes.TSS[,3]<-genes$TSS
Genes.TSS$TSSstartSpan<-genes$TSSstart
Genes.TSS$TSSendSpan<-genes$TSSend
Genes.TSS[-which(Genes.TSS$V1=='chrM'),]-> Genes.TSS

overlap<-bedTools.2in(bed1= bedTools.2sort(bed1=master), bed2=bedTools.2sort(bed1=Genes.TSS[,c(1,8,9,4,5,6,7,2,3)]), opt.string='')
unique(overlap$V4)->yes

master$YesCorrecte<-'NO'	
master$YesCorrecte[which(master$PeakID %in% yes)]<-'YES'

### ---> Gene annotation

bedTools.2closest(bed1= bedTools.2sort(bed1=master),bed2=bedTools.2sort(bed1=Genes.TSS),opt.string='-wb -d') ->cross.all
colapseTable(cross.all) -> colapsed
bedTools.2sort(bed1= colapsed) -> colapsed
colapsed->colapsedAllGenes	

colapsedAllGenes[c(1:18,33,34,44:46)] ->newtable
colnames(master)[1:18]->colnames(newtable)[1:18]
colnames(newtable)[19:23]<-c("H3K4me3","YesCorrecte","bp2TSS","allOverlapAnyENSEMBL","allOverlapAnyGSymbol")
write.table(newtable, file='MASTER_2017_ANNOT.bed', row.names=FALSE, quote=FALSE, sep='\t')

newtable->master
master[-which(master$TSS=="YES" | master$H3K4me3=='YES'),'PeakID']->ENH     		## OLD TSS per gene      
master[-which(master$YesCorrecte=="YES" | master$H3K4me3=='YES'),'PeakID']->ENH2 	## NEW TSS per transcript

### DESEQ
## Get master table

masterRep->krep
rownames(krep)<-krep$PeakID

FetalID<-c( "HSB100.DFC","HSB102.DFC","HSB240.DFC","HSB100.CBC","HSB102.CBC","HSB240.CBC" )
AdultID<-c( "HSB126.DFC","HSB136.DFC","HSB187.DFC","HSB126.CBC","HSB136.CBC","HSB187.CBC" )

# Classify enhancerns and promoters	
ad.dfc<-which(colnames(krep) %in% c("HSB126.DFC","HSB136.DFC","HSB187.DFC"))
fet.dfc<-which(colnames(krep) %in% c("HSB100.DFC","HSB102.DFC","HSB240.DFC"))		
ad.cbc<-which(colnames(krep) %in% c("HSB126.CBC","HSB136.CBC","HSB187.CBC"))
fet.cbc<-which(colnames(krep) %in% c("HSB100.CBC","HSB102.CBC","HSB240.CBC"))

k.enh<-transform(krep[,which(colnames(krep) %in% c("PeakID", "chr", "start", "end", "GSymbol", "ENSEMBL_PC","bp2TSS"))], 
                 DFC.A=apply(krep[,ad.dfc], 1, sum), CBC.A=apply(krep[,ad.cbc], 1, sum),
                 DFC.F=apply(krep[,fet.dfc], 1, sum), CBC.F=apply(krep[,fet.cbc], 1, sum), 
                 TOTAL.DFC=apply(krep[,c(ad.dfc,fet.dfc)], 1, sum), TOTAL.CBC=apply(krep[,c(ad.cbc,fet.cbc)], 1, sum),
                 TOTAL.FETAL=apply(krep[,c(fet.dfc,fet.cbc)], 1, sum), TOTAL.ADULT= apply(krep[,c(ad.dfc,ad.cbc)], 1, sum))

k.enh[which(k.enh$DFC.A >=2 & k.enh$DFC.F==0),] ->k.enh.dfc.adult            ### 0 0 0 X 1 1
k.enh[which(k.enh$DFC.F >=2 & k.enh$DFC.A==0),] ->k.enh.dfc.fetal  			 ### 1 1 X 0 0 0	
k.enh[which(k.enh$CBC.A >=2 & k.enh$CBC.F==0),] ->k.enh.cbc.adult            ### 0 0 0 X 1 1
k.enh[which(k.enh$CBC.F >=2 & k.enh$CBC.A==0),] ->k.enh.cbc.fetal  			 ### 1 1 X 0 0 0

k.enh[which(k.enh$DFC.A >=2 & k.enh$CBC.A==0),] ->k.enh.adult.dfc            ### 0 0 0 X 1 1
k.enh[which(k.enh$DFC.F >=2 & k.enh$CBC.F==0),] ->k.enh.fetal.dfc  			 ### 1 1 X 0 0 0	
k.enh[which(k.enh$CBC.A >=2 & k.enh$DFC.A==0),] ->k.enh.adult.cbc            ### 0 0 0 X 1 1
k.enh[which(k.enh$CBC.F >=2 & k.enh$DFC.F==0),] ->k.enh.fetal.cbc  			 ### 1 1 X 0 0 0

list(k.enh.dfc.adult, k.enh.dfc.fetal, k.enh.cbc.adult, k.enh.cbc.fetal, k.enh.adult.dfc, k.enh.fetal.dfc, k.enh.adult.cbc, k.enh.fetal.cbc)->l
names(l)<-c('DFC_Adult','DFC_Fetal','CBC_Adult','CBC_Fetal','AdultDFC','FetalDFC','AdultCBC','FetalCBC')
for (i in names(l)){
  write.table(l[[i]][,c(1:4)], file=paste('BEDs/',i,'_ALL_ENH.bed',sep=''), quote=F, row.names=F, col.names=F,sep='\t')
}

## cat *bed | sortBed | uniq >ALL.colapse	

## Run in cluster
# library(Rsubread)
# files<-read.table('FILE.list', header=F)$V1
# counts<-c()
# annD<-read.table('BEDs/ALL.colapse', header=F)
# colnames(annD)<-c('Chr','Start','End','GeneID')
# annD$Strand<-paste('+', c(1:dim(annD)[1]), sep='')		
# 
# fc_SE <- featureCounts(files, annot.ext=annD, isPairedEnd=TRUE, allowMultiOverlap=TRUE)
# counts<-fc_SE$counts
# colnames(counts)<-c("HSB100.CBC","HSB100.DFC","HSB102.CBC","HSB102.DFC","HSB126.CBC","HSB126.DFC",
# 					"HSB136.CBC","HSB136.DFC","HSB187.CBC","HSB187.DFC", "HSB240.CBC", "HSB240.DFC")
# save(counts, file='counts.allENH.LAST.rda')

files<-read.table('FILE.list', header=F)$V1
counts<-c()
annD<-read.table('ENHANCERS.collapsed.simple.bed', header=F)
colnames(annD)<-c('Chr','Start','End')
annD$GeneID<-paste('NormPeak', c(1:dim(annD)[1]), sep='')
annD$Strand<-paste('+', c(1:dim(annD)[1]), sep='')		

fc_SE <- featureCounts(files, annot.ext=annD, isPairedEnd=TRUE)
counts<-fc_SE$counts
colnames(counts)<-c("HSB100.CBC","HSB100.DFC","HSB102.CBC","HSB102.DFC","HSB126.CBC","HSB126.DFC",
                    "HSB136.CBC","HSB136.DFC","HSB187.CBC","HSB187.DFC", "HSB240.CBC", "HSB240.DFC")		
save(counts, file='counts.InAnyPeak.rda')

counts.est = estimateSizeFactorsForMatrix( counts )
rownames(counts)<-paste0("PK",c(1:nrow(counts)))
counts.all<-counts

#load<-load('TABLES/counts.allENH.LAST.rda')
#counts.all<-get(load)
#rm(load)

### Test DFC and CBC independently
DFC.peaks<-c(rownames(l[['DFC_Fetal']]), rownames(l[['DFC_Adult']]))
CBC.peaks<-c(rownames(l[['CBC_Fetal']]), rownames(l[['CBC_Adult']]))

### Test Fetal and Adult independently
Fetal.peaks<-c(rownames(l[['FetalDFC']]), rownames(l[['FetalCBC']]))
Adult.peaks<-c(rownames(l[['AdultDFC']]), rownames(l[['AdultCBC']]))

#rownames(counts)<-rownames(all.cnt)
counts.UnvarDFC<-counts[,grep('DFC', colnames(counts))]
counts.UnvarCBC<-counts[,grep('CBC', colnames(counts))]				
counts.DFC<-counts.all[which(rownames(counts.all) %in% DFC.peaks),grep('DFC', colnames(counts.all))]		
counts.est.DFC<- counts.est[grep('DFC', colnames(counts.all))]
counts.CBC<-counts.all[which(rownames(counts.all) %in% CBC.peaks),grep('CBC', colnames(counts.all))]
counts.est.CBC<- counts.est[grep('CBC', colnames(counts.all))]

counts.Fetal<-counts.all[which(rownames(counts.all) %in% Fetal.peaks),colnames(counts.all) %in% FetalID]		
counts.est.Fetal<- counts.est[colnames(counts.all) %in% FetalID]
counts.Adult<-counts.all[which(rownames(counts.all) %in% Adult.peaks),colnames(counts.all) %in% AdultID]		
counts.est.Adult<- counts.est[colnames(counts.all) %in% AdultID]

colData.DFC<-data.frame(Period=c(rep('Fetal',2), rep('Adult',3), rep('Fetal',1)), row.names=colnames(counts.DFC))
colData.CBC<-data.frame(Period=c(rep('Fetal',2), rep('Adult',3), rep('Fetal',1)), row.names=colnames(counts.CBC))

colData.Fetal<-data.frame(Period=c("CBC","DFC","CBC","DFC","CBC","DFC"), row.names=colnames(counts.Fetal))
colData.Adult<-data.frame(Period=c("CBC","DFC","CBC","DFC","CBC","DFC"), row.names=colnames(counts.Adult))

#dds <- DESeqDataSetFromMatrix(countData = counts.all,  colData = colData, design = ~ Region + Region:ind.n + Region:Period)
#sizeFactors(dds)<-counts.est

dds.DFC <- DESeqDataSetFromMatrix(countData = counts.DFC,  colData = colData.DFC, design = ~Period)
sizeFactors(dds.DFC)<-counts.est.DFC
dds.CBC <- DESeqDataSetFromMatrix(countData = counts.CBC,  colData = colData.CBC, design = ~Period)
sizeFactors(dds.CBC)<-counts.est.CBC

dds.unvar.DFC<-DESeqDataSetFromMatrix(countData = counts.UnvarDFC,  colData = colData.DFC, design = ~Period)
sizeFactors(dds.unvar.DFC)<-counts.est.DFC
dds.unvar.CBC<-DESeqDataSetFromMatrix(countData = counts.UnvarCBC,  colData = colData.CBC, design = ~Period)
sizeFactors(dds.unvar.CBC)<-counts.est.CBC

dds.Fetal <- DESeqDataSetFromMatrix(countData = counts.Fetal,  colData = colData.Fetal, design = ~Period)
sizeFactors(dds.Fetal)<-counts.est.Fetal
dds.Adult <- DESeqDataSetFromMatrix(countData = counts.Adult,  colData = colData.Adult, design = ~Period)
sizeFactors(dds.Adult)<-counts.est.Adult

#res<-DESeq(dds)
res.DFCtest<-DESeq(dds.DFC)
res.CBCtest<-DESeq(dds.CBC)

res.DFCtest.u<-DESeq(dds.unvar.DFC)
res.CBCtest.u<-DESeq(dds.unvar.CBC)

res.FetalTest<-DESeq(dds.Fetal)
res.AdultTest<-DESeq(dds.Adult)

#results(res, contrast=list("RegionCBC.PeriodFetal"))->res.CBC	
#results(res, contrast=list("RegionDFC.PeriodFetal"))->res.DFC	

res.DFC<-results(res.DFCtest)
res.CBC<-results(res.CBCtest)
res.uDFC<-results(res.DFCtest.u)
res.uCBC<-results(res.CBCtest.u)
res.Fetal<-results(res.FetalTest)
res.Adult<-results(res.AdultTest)

res.DFC[which(res.DFC$log2FoldChange <= -1 & res.DFC$padj < 0.01),]-> res.DFC.Adult
res.DFC[which(res.DFC$log2FoldChange >= 1 & res.DFC$padj < 0.01),]-> res.DFC.Fetal	
res.CBC[which(res.CBC$log2FoldChange <= -1 & res.CBC$padj < 0.01),]-> res.CBC.Adult	
res.CBC[which(res.CBC$log2FoldChange >= 1 & res.CBC$padj < 0.01),]-> res.CBC.Fetal	

res.uDFC[which(res.uDFC$pvalue>=0.05),]->res.DFC.unvar
res.uCBC[which(res.uCBC$pvalue>=0.05),]->res.CBC.unvar

res.Fetal[which(res.Fetal$log2FoldChange <= -1 & res.Fetal$padj < 0.01),]-> res.Fetal.CBC
res.Fetal[which(res.Fetal$log2FoldChange >= 1 & res.Fetal$padj < 0.01),]-> res.Fetal.DFC	
res.Adult[which(res.Adult$log2FoldChange <= -1 & res.Adult$padj < 0.01),]-> res.Adult.CBC	
res.Adult[which(res.Adult$log2FoldChange >= 1 & res.Adult$padj < 0.01),]-> res.Adult.DFC	

deseq.list<-list( res.DFC.Adult, res.DFC.Fetal, res.CBC.Adult, res.CBC.Fetal, 
                  res.Fetal.CBC,res.Fetal.DFC,res.Adult.CBC,res.Adult.DFC)
names(deseq.list)<-c('DFC_Adult','DFC_Fetal','CBC_Adult','CBC_Fetal',
                     'FetalCBC','FetalDFC','AdultCBC','AdultDFC')

#### We have the significant peaks.
#### Overlap with qualitative classification
M<-matrix(ncol=8, nrow=8)
rownames(M) <- names(l)
colnames(M) <- names(deseq.list)

definitive<-list()	
for(i in names(l)){
  for(j in names(deseq.list)){
    M[i,j]<-length(which(rownames(l[[i]]) %in% rownames(deseq.list[[j]])))
    if(i==j){ definitive[[i]]<-l[[i]][which(rownames(l[[i]]) %in% rownames(deseq.list[[j]])),] }
  }
}

### Use DESeq2 1.10.1 to obtain the same results in R3.2
### This version is with only temporal enhancers counts

#biasedPeaks -> definitive ## Activate this to run with the subset of temporal enhancers
definitive[['All']]<-krep
definitive[['UnvarDFC']]<-k.enh[which(rownames(k.enh) %in% rownames(res.DFC.unvar)),]
definitive[['UnvarCBC']]<-k.enh[which(rownames(k.enh) %in% rownames(res.CBC.unvar)),]
definitive[['DFC']]<-k.enh[which(k.enh$TOTAL.DFC>=2),]
definitive[['CBC']]<-k.enh[which(k.enh$TOTAL.CBC>=2),]
definitive[['Fetal']]<-k.enh[which(k.enh$TOTAL.FETAL>=2),]
definitive[['Adult']]<-k.enh[which(k.enh$TOTAL.ADULT>=2),]

#### Annotate ENH and PRO
Me<-read.table("1_MASTER_2017_ANNOT.bed", header=T)
filter(Me, TSS=="NO", H3K4me3=="NO")$PeakID -> ENH
filter(Me, TSS=="YES")$PeakID -> PRO

definitive.PRO<-list()
definitive.ENH<-list()
for (i in names(definitive)){
  definitive[[i]]<-krep[which(rownames(krep) %in% rownames(definitive[[i]])),]
  write.table(definitive[[i]][,c(1:4)], file=paste('TABLES/BEDs/BED_COUNTING_2018_Re/',i,'.bed',sep=''), sep='\t', quote=F, row.names=F, col.names=F)
  definitive.ENH[[i]]<-definitive[[i]][which(rownames(definitive[[i]]) %in% ENH),]
  write.table(definitive.ENH[[i]][,c(1:4)], file=paste('TABLES/BEDs/BED_COUNTING_2018_Re/',i,'.ENH.bed',sep=''), sep='\t', quote=F, row.names=F, col.names=F)
  definitive.PRO[[i]]<-definitive[[i]][which(rownames(definitive[[i]]) %in% PRO),]
  write.table(definitive.PRO[[i]][,c(1:4)], file=paste('TABLES/BEDs/BED_COUNTING_2018_Re/',i,'.PRO.bed',sep=''), sep='\t', quote=F, row.names=F, col.names=F)
}		

save(definitive,definitive.ENH,definitive.PRO, file='definitive.2018_edition_AllMatrix.rda') ## Depending on the counts
save(definitive,definitive.ENH,definitive.PRO, file='definitive.2018_edition.rda')	
### Careful, this objects has been updated by reannotating TSS.

### optional remove embryo and infant signal

#cat HSB1296-CBC-K27ac_peaks.broadPeak.bed HSB1472-CBC-K27ac_peaks.broadPeak.bed HSB872-CBC-K27ac_peaks.broadPeak.bed | sortBed | mergeBed >InfantCBC.bed
#cat HSB1296-DFC-K27ac_peaks.broadPeak.bed HSB1472-DFC-K27ac_peaks.broadPeak.bed HSB872-DFC-K27ac_peaks.broadPeak.bed | sortBed | mergeBed >InfantDFC.bed	
#cat CS16.K27ac_peaks.broadPeak.bed CS16b.K27ac_peaks.broadPeak.bed CS23.K27ac_peaks.broadPeak.bed CS23b.K27ac_peaks.broadPeak.bed F2bfrontal.K27ac_peaks.broadPeak.bed F2frontal.K27ac_peaks.broadPeak.bed | sortBed | mergeBed >Embrio.bed

#intersectBed -a MasterRep.bed -b InfantCBC.bed | cut -f4 | sort | uniq > InfantCBC.txt
#intersectBed -a MasterRep.bed -b InfantDFC.bed | cut -f4 | sort | uniq > InfantDFC.txt
#intersectBed -a MasterRep.bed -b Embrio.bed | cut -f4 | sort | uniq > Embrio.txt

InfantCBC<-read.table('TABLES/InfantCBC.txt', header=F)$V1
InfantDFC<-read.table('TABLES/InfantDFC.txt', header=F)$V1
Embrio<-read.table('TABLES/Embrio.txt', header=F)$V1

definitiveRefined<-list()
definitive[['DFC_Fetal']]->temp
definitiveRefined[['OnlyDFC_OnlyFetal']] <- temp[-which( temp$PeakID %in% k.enh[which(k.enh$TOTAL.CBC>=1),'PeakID'] | temp$PeakID %in% InfantDFC), ]		
definitive[['DFC_Adult']]->temp
definitiveRefined[['OnlyDFC_OnlyAdult']] <- temp[-which( temp$PeakID %in% k.enh[which(k.enh$TOTAL.CBC>=1),'PeakID'] | temp$PeakID %in% Embrio), ]
definitive[['CBC_Fetal']]->temp
definitiveRefined[['OnlyCBC_OnlyFetal']] <- temp[-which( temp$PeakID %in% k.enh[which(k.enh$TOTAL.DFC>=1),'PeakID'] | temp$PeakID %in% InfantCBC), ]
definitive[['CBC_Adult']]->temp
definitiveRefined[['OnlyCBC_OnlyAdult']] <- temp[-which( temp$PeakID %in% k.enh[which(k.enh$TOTAL.DFC>=1),'PeakID'] | temp$PeakID %in% Embrio), ]	
definitive[['UnvarDFC']]->temp
definitiveRefined[['OnlyDFC_Unvar']] <- temp[-which( temp$PeakID %in% k.enh[which(k.enh$TOTAL.CBC>=1),'PeakID'] ), ]
definitive[['UnvarCBC']]->temp
definitiveRefined[['OnlyCBC_Unvar']] <- temp[-which( temp$PeakID %in% k.enh[which(k.enh$TOTAL.DFC>=1),'PeakID'] ), ]		

for (i in names(definitiveRefined)){
  write.table(definitiveRefined[[i]][,c(1:4)], file=paste('BEDs/',i,'.bed',sep=''), sep='\t', quote=F, row.names=F, col.names=F)
  write.table(rownames(definitiveRefined[[i]]), file=paste('BEDs/',i,'.list',sep=''), quote=F, row.names=F, col.names=F)
  temp<-definitiveRefined[[i]][which(definitiveRefined[[i]][,'PeakID'] %in% ENH),  ]
  write.table(temp[,c(1:4)], file=paste('BEDs/',i,'.ENH.bed',sep=''), sep='\t', quote=F, row.names=F, col.names=F)
  write.table(rownames(temp), file=paste('BEDs/',i,'.ENH.list',sep=''), sep='\t', quote=F, row.names=F, col.names=F)
  temp<-definitiveRefined[[i]][which(definitiveRefined[[i]][,'PeakID'] %in% PRO),  ]
  write.table(temp[,c(1:4)], file=paste('BEDs/',i,'.PRO.bed',sep=''), sep='\t', quote=F, row.names=F, col.names=F)
  write.table(rownames(temp), file=paste('BEDs/',i,'.PRO.list',sep=''), sep='\t', quote=F, row.names=F, col.names=F)
}	

for (i in names(definitive)){
  write.table(definitive[[i]][,c(1:4)], file=paste('BEDs/',i,'.bed',sep=''), sep='\t', quote=F, row.names=F, col.names=F)
  write.table(rownames(definitive[[i]]), file=paste('BEDs/',i,'.list',sep=''), quote=F, row.names=F, col.names=F)
  temp<-definitive[[i]][which(definitive[[i]][,'PeakID'] %in% ENH),  ]
  write.table(temp[,c(1:4)], file=paste('BEDs/',i,'.ENH.bed',sep=''), sep='\t', quote=F, row.names=F, col.names=F)
  write.table(rownames(temp), file=paste('BEDs/',i,'.ENH.list',sep=''), sep='\t', quote=F, row.names=F, col.names=F)
  temp<-definitive[[i]][which(definitive[[i]][,'PeakID'] %in% PRO),  ]
  write.table(temp[,c(1:4)], file=paste('BEDs/',i,'.PRO.bed',sep=''), sep='\t', quote=F, row.names=F, col.names=F)
  write.table(rownames(temp), file=paste('BEDs/',i,'.PRO.list',sep=''), sep='\t', quote=F, row.names=F, col.names=F)
}	

write.table(file='BEDs/AnyDFC.list', k.enh[which(k.enh$TOTAL.DFC>=1),'PeakID'], sep='\t', quote=F, row.names=F, col.names=F)
write.table(file='BEDs/AnyCBC.list', k.enh[which(k.enh$TOTAL.CBC>=1),'PeakID'], sep='\t', quote=F, row.names=F, col.names=F)

### FUNCTIONS	
colapseTable<-function(t){
  dupPC<-t$V4[duplicated(t$V4)]
  t[which(t$V4 %in% dupPC),]->tableDUP
  t[-which(t$V4 %in% dupPC),]->tableNO_DUP
  tnew<-c()
  for (i in unique(tableDUP$V4)){
    #print(i)
    temp<-t[which(t$V4==i),]
    temp$newE<-paste(unique(temp$V39), collapse=',')
    temp$newG<-paste(unique(temp$V41), collapse=',')
    ## Show preference for protein coding genes
    if('protein_coding' %in% temp$V40){
      tnew<-rbind(tnew,temp[which(temp$V40 %in% 'protein_coding'),][1,])
    } else {
      tnew<-rbind(tnew,temp[1,])
    }
  }
  tableNO_DUP$newE<-tableNO_DUP$V39
  tableNO_DUP$newG<-tableNO_DUP$V41
  tablefinal<-rbind(tnew, tableNO_DUP)
  tablefinal<-tablefinal[order(tablefinal$V4),]
  return(tablefinal)
}

############# PRODUCE TABLE WITH EMBRYO AND INFANT SAMPLES

# for i in $(ls *K27ac_peaks.broadPeak); do cut -f1-4 $i | sed 's/-H3K27ac_peak_[0-9]*//g' | sed 's/-K27ac_peak_[0-9]*//g' >${i}.bed; done
# for i in $(cat InfantEnbrioList.txt); do name=$(basename $i); cut -f1-4 $i | sed 's/\.K27ac_peak_[0-9]*//g' | sed 's/-K27ac_peak_[0-9]*//g' >${name}.bed; done

# cat *broadPeak.bed | sortBed | mergeBed  -c 4 -o collapse -delim ' ' >../TABLES/ENHANCERS.InfEmb.collapsed.bed

### In cluster
### Produce master table
e<-read.table('ENHANCERS.InfEmb.collapsed.bed', header=F, sep='\t')
gsub('-','.',e$V4) -> e$V4
samples<-c(	"F2frontal","F2bfrontal","CS16","CS16b","CS23","CS23b",
            "HSB100.CBC","HSB100.DFC","HSB102.CBC","HSB102.DFC","HSB240.CBC","HSB240.DFC",
            "HSB1296.CBC","HSB1296.DFC","HSB1472.CBC","HSB1472.DFC","HSB872.CBC","HSB872.DFC",
            "HSB126.CBC","HSB126.DFC","HSB136.CBC","HSB136.DFC","HSB187.CBC","HSB187.DFC")

matrix(nrow=dim(e)[1], ncol=length(samples)) -> bin
bin[is.na(bin)]<-0
colnames(bin)<-samples
for (i in c(1:dim(e)[1])){	
  mask<-which(samples %in% unlist(strsplit(e[i,]$V4, " ")))
  bin[i,mask]<-1
}
cbind(e[,c(1:3)], bin)->master
paste('PK',c(1:dim(bin)[1]), sep='')->master$PeakID

## Reproducible enhancers
master[-which(rowSums(master[,c(4:27)])==1),]->masterRep
write.table(masterRep[,c(1,2,3,28)], file='MasterRep.2018.bed', row.names=F, col.names=F, quote=F, sep='\t')

source('~/Desktop/bin/BEDtoolsR.R')

### ---> Overlap with TSS
overlap<-bedTools.2in(bed1= bedTools.2sort(bed1=master), bed2=bedTools.2sort(bed1=Genes.TSS[,c(1,8,9,4,5,6,7,2,3)]), opt.string='-wb')
unique(overlap$V28)->yes

master$YesCorrecte<-'NO'	
master$YesCorrecte[which(master$PeakID %in% yes)]<-'YES'

master[-which(master$YesCorrecte =="YES" | master$H3K4me3=='YES'),'PeakID']->ENH     			
master[which(master$YesCorrecte =="YES"),'PeakID']->PRO

### Gene annotation

bedTools.2closest(bed1= bedTools.2sort(bed1=master),bed2=bedTools.2sort(bed1=Genes.TSS),opt.string='-wb -d') ->cross.all
colapseTable(cross.all) -> colapsed
bedTools.2sort(bed1= colapsed) -> colapsed
colapsed->colapsedAllGenes	

colapsedAllGenes[c(1:3,28,4:27,29,30,40:42)] ->newtable
colnames(master)[4:27]->colnames(newtable)[5:28]
colnames(newtable)[29:33]<-c("H3K4me3","YesCorrecte","bp2TSS","allOverlapAnyENSEMBL","allOverlapAnyGSymbol")
write.table(newtable, file='MASTER.IE_2018_ANNOT.bed', row.names=FALSE, quote=FALSE, sep='\t')
write.table(newtable[,c(1:4)], file='MASTER.IE_2018_ANNOT_Nocols.bed',row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')

### Categories	

master->krep
ad.dfc<-which(colnames(krep) %in% c("HSB126.DFC","HSB136.DFC","HSB187.DFC"))
fet.dfc<-which(colnames(krep) %in% c("HSB100.DFC","HSB102.DFC","HSB240.DFC"))
inf.dfc<-which(colnames(krep) %in% c("HSB1296.DFC","HSB872.DFC","HSB1472.DFC"))		
ad.cbc<-which(colnames(krep) %in% c("HSB126.CBC","HSB136.CBC","HSB187.CBC"))
fet.cbc<-which(colnames(krep) %in% c("HSB100.CBC","HSB102.CBC","HSB240.CBC"))
inf.cbc<-which(colnames(krep) %in% c("HSB1296.CBC","HSB872.CBC","HSB1472.CBC"))
emb<-which(colnames(krep) %in% c("F2frontal","F2bfrontal","CS16","CS16b","CS23","CS23b"))

k<-transform(krep[,which(colnames(krep) %in% c("PeakID", "V1", "V2", "V3"))], 
             DFC.A=apply(krep[,ad.dfc], 1, sum), CBC.A=apply(krep[,ad.cbc], 1, sum),
             DFC.F=apply(krep[,fet.dfc], 1, sum), CBC.F=apply(krep[,fet.cbc], 1, sum),
             DFC.I=apply(krep[,inf.dfc], 1, sum), CBC.I=apply(krep[,inf.cbc], 1, sum),
             EMB=apply(krep[,emb], 1, sum), 
             TOTAL.DFC=apply(krep[,c(ad.dfc,fet.dfc,inf.dfc)], 1, sum), 
             TOTAL.CBC=apply(krep[,c(ad.cbc,fet.cbc,inf.cbc)], 1, sum),
             TOTAL.AD=apply(krep[,c(ad.cbc,ad.dfc)], 1, sum),
             TOTAL.FET=apply(krep[,c(fet.cbc,fet.dfc)], 1, sum),
             TOTAL.INF=apply(krep[,c(inf.cbc,inf.dfc)], 1, sum),
             PRENATAL=apply(krep[,c(emb,fet.dfc,fet.cbc)], 1, sum),
             POSTNATAL=apply(krep[,c(inf.dfc,inf.cbc,ad.dfc,ad.cbc)], 1, sum)
)

save(k, file="ENHANCERS_PRESNECE_SAMPLES.rda")
l<-list()

l[['DFC']]<-k[which(k$TOTAL.DFC >=2),]
l[['CBC']]<-k[which(k$TOTAL.CBC >=2),]

l[['Adult']]<-k[which(k$TOTAL.AD >=2),]
l[['Fetal']]<-k[which(k$TOTAL.FET >=2),]
l[['Infant']]<-k[which(k$TOTAL.INF >=2),]
l[['Embryonic']]<-k[which(k$EMB >=2),]

l[['AdultDFC']]<-k[which(k$DFC.A >=2),]
l[['AdultCBC']]<-k[which(k$CBC.A >=2),]
l[['FetalDFC']]<-k[which(k$DFC.F >=2),]
l[['FetalCBC']]<-k[which(k$CBC.F >=2),]
l[['InfantDFC']]<-k[which(k$DFC.I >=2),]
l[['InfantCBC']]<-k[which(k$CBC.I >=2),]

l[['AdultDFCsp']]<-k[which(k$DFC.A >=2 & k$CBC.A == 0),]
l[['AdultCBCsp']]<-k[which(k$CBC.A >=2 & k$DFC.A== 0),]
l[['FetalDFCsp']]<-k[which(k$DFC.F >=2 & k$CBC.F == 0),]
l[['FetalCBCsp']]<-k[which(k$CBC.F >=2 & k$DFC.F== 0),]
l[['InfantDFCsp']]<-k[which(k$DFC.I >=2 & k$CBC.I == 0),]
l[['InfantCBCsp']]<-k[which(k$CBC.I >=2 & k$DFC.I== 0),]

l[['PostNatal_spDFC']]<-k[which(k$POSTNATAL >=2 & k$PRENATAL == 0 & k$TOTAL.DFC >=2),]
l[['PostNatal_spCBC']]<-k[which(k$POSTNATAL >=2 & k$PRENATAL == 0 & k$TOTAL.CBC >=2),]
l[['PreNatal_spDFC']]<-k[which(k$PRENATAL >=2 & k$POSTNATAL == 0 & k$TOTAL.DFC >=2),]
l[['PreNatal_spCBC']]<-k[which(k$PRENATAL >=2 & k$POSTNATAL == 0 & k$TOTAL.CBC >=2),]

l[['PostNatal_spDFC_sp']]<-k[which(k$POSTNATAL >=2 & k$PRENATAL == 0 & k$TOTAL.DFC >=2 & k$TOTAL.CBC == 0),]
l[['PostNatal_spCBC_sp']]<-k[which(k$POSTNATAL >=2 & k$PRENATAL == 0 & k$TOTAL.CBC >=2 & k$TOTAL.DFC == 0),]
l[['PreNatal_spDFC_sp']]<-k[which(k$PRENATAL >=2 & k$POSTNATAL == 0 & k$TOTAL.DFC >=2 & k$TOTAL.CBC == 0),]
l[['PreNatal_spCBC_sp']]<-k[which(k$PRENATAL >=2 & k$POSTNATAL == 0 & k$TOTAL.CBC >=2 & k$TOTAL.DFC == 0),]

for(i in names(l)){
  l[[i]] -> temp
  write.table(temp[,c(1:4)], file=paste(i,'.ENH_PRO.bed',sep=''), sep='\t', quote=F, row.names=F, col.names=F)
  temp<-l[[i]][which(l[[i]][,'PeakID'] %in% ENH),  ]
  print(c(i, dim(temp)))
  #write.table(temp[,c(1:4)], file=paste(i,'.ENH.bed',sep=''), sep='\t', quote=F, row.names=F, col.names=F)
  temp<-l[[i]][which(l[[i]][,'PeakID'] %in% PRO),  ]
  print(c(i, dim(temp)))
  #write.table(temp[,c(1:4)], file=paste(i,'.PRO.bed',sep=''), sep='\t', quote=F, row.names=F, col.names=F)
}

### Functions
colapseTable<-function(t){
  dupPC<-t$V28[duplicated(t$V28)]
  t[which(t$V28 %in% dupPC),]->tableDUP
  t[-which(t$V28 %in% dupPC),]->tableNO_DUP
  tnew<-c()
  for (i in unique(tableDUP$V28)){
    #print(i)
    temp<-t[which(t$V28==i),]
    temp$newE<-paste(unique(temp$V35), collapse=',')
    temp$newG<-paste(unique(temp$V37), collapse=',')
    ## Show preference for protein coding genes
    if('protein_coding' %in% temp$V36){
      tnew<-rbind(tnew,temp[which(temp$V36 %in% 'protein_coding'),][1,])
    } else {
      tnew<-rbind(tnew,temp[1,])
    }
  }
  tableNO_DUP$newE<-tableNO_DUP$V35
  tableNO_DUP$newG<-tableNO_DUP$V37
  tablefinal<-rbind(tnew, tableNO_DUP)
  tablefinal<-tablefinal[order(tablefinal$V28),]
  return(tablefinal)
}

#### Annotated GWAS hits with Hi-C data 
### This version uses Geschwind loop data

### 1--> Load SNP annotation from Biomart Version 78 (equivalent to Gencode V21)
snp<-read.table("../CLUMPS/ALL_SNPS_CODES.ANNOT.txt", header=FALSE, sep="\t")
colnames(snp)<-c("id","ens","Variant.consequence","chr")
filter(snp, chr %in% c(1:22,"X","Y")) -> snp
unique(snp$Variant.consequence) -> cat
cat[!cat %in% c("upstream_gene_variant","intron_variant","downstream_gene_variant","")] -> cat

### 2--> Load loops from Geschwind

fetal<-read.table("Fetalbrain_TSSinteraction.txt", header=TRUE)
adult<-read.table("Adultbrain_TSSinteraction.txt", header=TRUE)

f<-bedTools.2sort(bed1 =  fetal[,c(1,3,4,5,2)])
a<-bedTools.2sort(bed1 =  adult[,c(1,3,4,5,2)])

f$time<-"Fetal"
a$time<-"Adult"

data<-bedTools.2sort(bed1=rbind(f,a))
write.table(data, file="gerschwind_loops.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

data$ID<-paste0("ID_",c(1:nrow(data)))

### 3--> Input Sestan lab Chip-Seq data
### Liftover coordinates to hg19 to produce MASTER.IE_2018.hg19.bed

load("ENHANCERS_PRESNECE_SAMPLES.rda")
master<-read.table("MASTER.IE_2018.hg19.bed")
Mhg38<-read.table("MASTER.IE_2018_ANNOT.bed", header=TRUE)

dplyr::select(Mhg38, V28, allOverlapAnyENSEMBL,allOverlapAnyGSymbol ) -> Mgene
merge(master, Mgene, by.x='V4', by.y='V28') %>% dplyr::select(V1, V2, V3, V4,allOverlapAnyENSEMBL, allOverlapAnyGSymbol) -> M

filter(Mhg38, H3K4me3=='NO', YesCorrecte=='NO') %>% dplyr::select(V28) -> enh
filter(Mhg38, YesCorrecte=='YES') %>% dplyr::select(V28) %>% unlist -> prom

M$Emb<-"NO"
M$PreDFC<-"NO"
M$PostDFC<-"NO"
M$PreCBC<-"NO"
M$PostCBC<-"NO"

M$Emb[which( M$V4 %in% filter(k, EMB>0)$PeakID )]<-"YES"
M$PreDFC[which( M$V4 %in% filter(k, DFC.F>0)$PeakID )]<-"YES"
M$PostDFC[which( M$V4 %in% filter(k, DFC.I>0 | DFC.A>0)$PeakID )]<-"YES"
M$PreCBC[which( M$V4 %in% filter(k, CBC.F>0)$PeakID )]<-"YES"
M$PostCBC[which( M$V4 %in% filter(k, CBC.I>0 | CBC.A>0)$PeakID )]<-"YES"

write.table(bedTools.2sort(bed1=M[M$Emb=="YES",]), file="M_emb.bed", sep="\t", quote=F, row.names=F, col.names = F)
write.table(bedTools.2sort(bed1=M[M$PreDFC=="YES",]), file="M_PreDFC.bed", sep="\t", quote=F, row.names=F, col.names = F)
write.table(bedTools.2sort(bed1=M[M$PreCBC=="YES",]), file="M_PreCBC.bed", sep="\t", quote=F, row.names=F, col.names = F)
write.table(bedTools.2sort(bed1=M[M$PostDFC=="YES",]), file="M_PostDFC.bed", sep="\t", quote=F, row.names=F, col.names = F)
write.table(bedTools.2sort(bed1=M[M$PostCBC=="YES",]), file="M_PostCBC.bed", sep="\t", quote=F, row.names=F, col.names = F)

### 4--> Overlap tss with open chromatin

dataTSS<-data
dataTSS$start<-data$V5-1000
dataTSS$end<-data$V5+1001
dataTSS<-dataTSS[,c(1,8,9,4,6,7,2,3)]
dataTSS<-bedTools.2sort(bed1=dataTSS)
open_tss<-bedTools.2in( bed1=dataTSS, bed2=M,  opt.string='-wa -wb' )

### 5--> Filter interactions to exclude closed chromatin TSS

data.filt<-open_tss[,c(1,7,8,6,4,5,12,15:19)]

### Filter by HiC specificity
bedTools.2in( bed1=M, bed2=bedTools.2sort(bed1=data.filt),  opt.string='-wa -wb' ) -> loops
loops[,c(4,5,6,7:11,15,12:14,16:23)] -> loops
colnames(loops)<-paste0("V",c(1:ncol(loops)))

#sapply(c(1:nrow(loops)), function(x) { sum(which(loops[x,c(4:8)] == "YES") == which(loops[x,c(16:20)] == "YES"))}) -> concordant
#enh$V28[!enh$V28 %in% unique(loops$V4)] -> peaks.No.interact

### 6--> Annotate enhancer table with TSS only with openChrom
TSS<-read.table("../SULLIVAN/Gencode/Genes.TSS.hg19.bed")
TSS<-bedTools.2sort(bed1=TSS)
TSS$start<-TSS$V2-1000
TSS$end<-TSS$V3+1000

TSS->converter
converter$Ens<-sapply(strsplit(as.character(TSS$V4), split=".", fixed=T),"[", 1)
converter$Sym<-sapply(strsplit(as.character(TSS$V4), split="|", fixed=T),"[", 2)

oTSS<-bedTools.2in( bed1=TSS[,c(1,5,6,4)], bed2=M,  opt.string='-wa -wb' )
oTSS<-oTSS[,c(1:4,8,11:15)]

bedTools.2sort(bed1=oTSS)->oTSS
save(oTSS, file="oTSS.rda")

ClosestTSS<-function(peak){
  filter(M, V4==peak) -> row
  filter(oTSS, as.character(V1) == row$V1) -> possibleTSS
  a<-which(row[c(7:11)] == "YES")
  which(rowSums(as.data.frame(possibleTSS[,c(6:10)][,a]) == "YES")>0) -> m
  bedTools.2closest(bed1= bedTools.2sort(bed1=row[,c(1:3)]),bed2=possibleTSS[m,],opt.string='-wb -d') ->cros
  return(cros)
}

#### Annotate GWAS loci and produce gene lists

for (disease in c("IQ_2018", "ADHD", "ASD", "CLOZUK", "HbA1c", "NEUROTICISM","HEIGHT", "AD","BD","MDD","PD")) { 
  
  print(paste0("Annotating clups for ",disease, "........."))
  clumps<-read.table(paste0("../CLUMPS/Clumps/",disease,".clumps.bed"), header=TRUE)
  bedTools.2merge(bed1=bedTools.2sort(bed1=clumps)) -> clumps
  bedTools.2in(bed1=bedTools.2sort(bed1=clumps), bed2=bedTools.2sort(bed1=M), opt.string='-wb') %>% dplyr::select (V7,V8,V9) -> peaks_snps
  unique(peaks_snps$V7) -> peaks2eval
  
  write.table(bedTools.2sort(bed1=clumps), file=paste0(disease,'.bed'), quote=F, col.names=F, row.names=F)
  
  ### Interacting peaks
  print(paste0("Annotating interactions for ",disease, "........."))
  filter(loops,V1 %in% peaks2eval) -> l
  ## Filter discordant peaks for ChipSeq
  sapply(c(1:nrow(l)), function(x) { sum(which(l[x,c(4,5,7)] == "YES") %in% which(l[x,c(16,17,19)] == "YES"))}) -> concordantFet
  sapply(c(1:nrow(l)), function(x) { sum(which(l[x,c(6,8)] == "YES") %in% which(l[x,c(18,20)] == "YES"))}) -> concordantAd
  ## Filter discordant HiC
  #sapply(c(1:nrow(l)), function(x) { sum(l[x,c(4,5,7,16,17,19)] == "YES")}) -> preLoop
  #sapply(c(1:nrow(l)), function(x) { sum(l[x,c(6,8,18,20)] == "YES")}) -> postLoop
  mask<-c(which(l$V14=="Adult" & concordantAd==0),which(l$V14=="Fetal" & concordantFet==0))
  if(length(mask)>0){
    genes<-unique(l[-mask,]$V13)
    peakInteract<-unique(l[-mask,]$V1)
  } else { 
    genes<-unique(l$V13)
    peakInteract<-unique(l$V1)
  }
  
  HiC_Ensmbl.HiC<-unique(sapply(strsplit(as.character(genes), split="|", fixed=T),"[", 1))
  HiC_Gsymbol.HiC<-unique(sapply(strsplit(as.character(genes), split="|", fixed=T),"[", 2))
  HiC.Ensmbl<-data.frame(gene=HiC_Ensmbl.HiC, case="HiC_Won")
  HiC.Gsymbol<-data.frame(gene=HiC_Gsymbol.HiC, case="HiC_Won")
  
  ### Non-interacting enhancer peaks
  print(paste0("Annotating non-interactions for ",disease, "........."))
  peaks2eval[!peaks2eval %in% peakInteract] -> peaks2evalInNonInteract
  peaks2evalInNonInteract<-peaks2evalInNonInteract[peaks2evalInNonInteract %in% enh$V28]
  NiInt.Ensmbl<-c()
  NiInt.Gsymbol<-c()
  if(length(peaks2evalInNonInteract) > 0) {
    cross.all.table<-c()
    for (i in peaks2evalInNonInteract){
      #print(i)
      cross.all.table<-rbind(cross.all.table, ClosestTSS(i))
    }
    HiC_Ensmbl.NoInt<-unique(sapply(strsplit(as.character(cross.all.table$V7), split="|", fixed=T),"[", 1))
    HiC_Gsymbol.NoInt<-unique(sapply(strsplit(as.character(cross.all.table$V7), split="|", fixed=T),"[", 2))
    NiInt.Ensmbl<-data.frame(gene=HiC_Ensmbl.NoInt, case="NoInt_Won")
    NiInt.Gsymbol<-data.frame(gene=HiC_Gsymbol.NoInt, case="NoInt_Won")
  }  
  ### Peaks in promotor
  peaks2eval[peaks2eval %in% prom] -> peaks2evalInProm
  filter(Mhg38, V28 %in% peaks2evalInProm) %>% dplyr::select(allOverlapAnyENSEMBL,allOverlapAnyGSymbol) -> peaks2evalInPromGenes
  HiC_Ensmbl.Prom<-unique(unlist(strsplit(as.character(peaks2evalInPromGenes$allOverlapAnyENSEMBL), ',')))
  HiC_Gsymbol.Prom<-unique(unlist(strsplit(as.character(peaks2evalInPromGenes$allOverlapAnyGSymbol), ',')))
  Prom.Ensmbl<-data.frame(gene=HiC_Ensmbl.Prom, case="Promoter")
  Prom.Gsymbol<-data.frame(gene=HiC_Gsymbol.Prom, case="Promoter")
  
  ### Exonic variants
  s<-read.table(paste0("../CLUMPS/",disease,".toClump.output.clumped.rsid"), header=FALSE)$V1
  #filter(snp, Biotype=="protein_coding",Variant.consequence %in% c("missense_variant", "stop_gained", "splice_acceptor_variant","splice_region_variant","splice_donor_variant","3_prime_UTR_variant","5_prime_UTR_variant")) -> temp
  filter(snp, id %in% s, Variant.consequence %in% cat) -> temp
  
  Exo.Ensmbl<-c()
  Exo.Gsymbol<-c()
  if(dim(temp)[1]>0){
    unique(temp$ens) -> ge
    unique(filter(converter, Ens %in% ge)$V4) -> exonicVariants
    print(paste0(disease, "...." ,mean(ge %in% converter$Ens), "...", ge[!ge %in% converter$Ens]))
    #print(paste0(disease, "...." ,mean(exonicVariantsSym %in% converter$Sym), "...", exonicVariantsSym[!exonicVariantsSym %in% converter$Sym]))
    HiC_Ensmbl.Exo<-unique(sapply(strsplit(as.character(exonicVariants), split="|", fixed=T),"[", 1))
    HiC_Gsymbol.Exo<-unique(sapply(strsplit(as.character(exonicVariants), split="|", fixed=T),"[", 2))
    Exo.Ensmbl<-data.frame(gene=HiC_Ensmbl.Exo, case="Genic")
    Exo.Gsymbol<-data.frame(gene=HiC_Gsymbol.Exo, case="Genic")
  }
  
  HiC_Ensmbl<-ddply( rbind(Exo.Ensmbl,  Prom.Ensmbl, NiInt.Ensmbl, HiC.Ensmbl)  ,'gene',summarize,X=paste(case,collapse=",") )
  HiC_Gsymbol<-ddply( rbind(Exo.Gsymbol,  Prom.Gsymbol, NiInt.Gsymbol, HiC.Gsymbol)  ,'gene',summarize,X=paste(case,collapse=",") )
  
  write.table(HiC_Gsymbol, file=paste0("GENES/",disease,"_NC2.Gsymbol.2018.Geschwind.PromOpen.E78.ClosestOpen.Concord.CORRECT2.txt"), quote=F, col.names=F, row.names=F)
  write.table(HiC_Ensmbl, file=paste0("GENES/",disease,"_NC2.Ensmbl.2018.Geschwind.PromOpen.ClosestOpen.E78.Concord.CORRECT2.txt"), quote=F, col.names=F, row.names=F)
}

### This version uses Sullivan's loops

library(dplyr)
library(plyr)
source('~/Desktop/bin/BEDtoolsR.R')
snp<-read.table("../CLUMPS/ALL_SNPS_CODES.ANNOT.txt", header=FALSE, sep="\t")
colnames(snp)<-c("id","ens","Variant.consequence","chr")
filter(snp, chr %in% c(1:22,"X","Y")) -> snp
unique(snp$Variant.consequence) -> cat
cat[!cat %in% c("upstream_gene_variant","intron_variant","downstream_gene_variant","")] -> cat

data<-read.table('sullivan.hic.loops.txt', header=TRUE)

### Unique identified
data$id<-paste0("ID_", 1:nrow(data))

### Separate left and right

dataL<-data[,c(1,2,8,4:7)]
dataR<-data[,c(1,3,8,4:7)]

### Expand bins to 10kb

dataL$start <- dataL$l1-5000
dataL$end   <- dataL$l1+5000
dataLbed    <- dplyr::select(dataL, l0c, start, end, id, Adult, Fetal, CP, GZ)
dataLbed$dir<-"LEFT"

dataR$start <- dataR$l2-5000
dataR$end   <- dataR$l2+5000
dataRbed    <- dplyr::select(dataR, l0c, start, end, id, Adult, Fetal, CP, GZ)
dataRbed$dir<-"RIGHT"

dataLwrite<-dataLbed
dataLwrite$ID<-paste0(dataLwrite$id,"--",dataLwrite$dir)
dataLwrite<-dataLwrite[,c(1:3,6)]
write.table(bedTools.2sort(bed1=dataLbed), file='dataL.bed', quote=F, col.names=F, row.names=F, sep="\t")
write.table(bedTools.2sort(bed1=dataRbed), file='dataR.bed', quote=F, col.names=F, row.names=F, sep="\t")

### Input Sestan lab Chip-Seq data

master<-read.table("MASTER.IE_2018.hg19.bed")
Mhg38<-read.table("MASTER.IE_2018_ANNOT.bed", header=TRUE)

dplyr::select(Mhg38, V28, allOverlapAnyENSEMBL,allOverlapAnyGSymbol ) -> Mgene
merge(master, Mgene, by.x='V4', by.y='V28') %>% dplyr::select(V1, V2, V3, V4,allOverlapAnyENSEMBL, allOverlapAnyGSymbol) -> M

filter(Mhg38, H3K4me3=='NO', YesCorrecte=='NO') %>% dplyr::select(V28) -> enh
filter(Mhg38, YesCorrecte=='YES') %>% dplyr::select(V28) %>% unlist -> prom

load("ENHANCERS_PRESNECE_SAMPLES.rda")
M$Emb<-"NO"
M$PreDFC<-"NO"
M$PostDFC<-"NO"
M$PreCBC<-"NO"
M$PostCBC<-"NO"

M$Emb[which( M$V4 %in% filter(k, EMB>0)$PeakID )]<-"YES"
M$PreDFC[which( M$V4 %in% filter(k, DFC.F>0)$PeakID )]<-"YES"
M$PostDFC[which( M$V4 %in% filter(k, DFC.I>0 | DFC.A>0)$PeakID )]<-"YES"
M$PreCBC[which( M$V4 %in% filter(k, CBC.F>0)$PeakID )]<-"YES"
M$PostCBC[which( M$V4 %in% filter(k, CBC.I>0 | CBC.A>0)$PeakID )]<-"YES"

### Overlap enh and HiC data

bedTools.2in(bed1=bedTools.2sort(bed1=M[,-c(5,6)]), bed2=bedTools.2sort(bed1=dataLbed), opt.string='-wa -wb') -> M.L
bedTools.2in(bed1=bedTools.2sort(bed1=M[,-c(5,6)]), bed2=bedTools.2sort(bed1=dataRbed), opt.string='-wa -wb') -> M.R

cbind(as.character(dataLbed$l0c), dataLbed$start, dataRbed$end) -> composite
bedTools.2merge( bed1= bedTools.2sort(bed1= rbind(dataLbed, dataRbed)  )) -> all

### Overlap with TSS

tss <- read.table('Gencode/Genes.TSS.hg19.bed')
tssProm <- tss 
tssProm$V2<-tss$V2-1000 
tssProm$V3<-tss$V3+1000 

### Require Tss to overlap a peak
bedTools.2in(bed1=bedTools.2sort(bed1=tssProm), bed2=bedTools.2sort(bed1=M), opt.string='-wa -wb') -> openProm
openProm[c(1:4,8,11:15)] -> openProm

bedTools.2in(bed1=bedTools.2sort(bed1=dataLbed), bed2=bedTools.2sort(bed1=openProm), opt.string='-wa -wb') -> PromL
bedTools.2in(bed1=bedTools.2sort(bed1=dataRbed), bed2=bedTools.2sort(bed1=openProm), opt.string='-wa -wb') -> PromR

write.table(PromL, file='PromL.open.bed', quote=F, col.names=F, row.names=F)
write.table(PromR, file='PromR.opne.bed', quote=F, col.names=F, row.names=F)

unique(PromL[,c(4:9,13:19)]) -> PromLred
unique(PromR[,c(4:9,13:19)]) -> PromRred
rbind(PromLred, PromRred) -> PromTOTAL

# ddply(PromLred,'V4',summarize,X=paste(V9,collapse=",") ) ->PromLredSum
# ddply(PromRred,'V4',summarize,X=paste(V9,collapse=",") ) ->PromRredSum
# PromLredSum$dir<-"LEFT"
# PromRredSum$dir<-"RIGHT"
# 
# rbind(PromLredSum, PromRredSum) -> PromTOTAL

### This part requieres that the end loop with the promoter overlap at least one peak
#merge(PromLredSum, M.L[,c(4,8)], by.x="V4", by.y="V8") -> peak.TSS.L
#merge(PromRredSum, M.R[,c(4,8)], by.x="V4", by.y="V8") -> peak.TSS.R
#colnames(peak.TSS.L)<-c("ID","GENE","DIR","PEAK")
#colnames(peak.TSS.R)<-c("ID","GENE","DIR","PEAK")

# merge(M.L[,c(4,8,9)], PromRredSum, by.x="V8", by.y="V4") %>% dplyr::select(V4) -> all.peaks.interactions.L
# #unique(c(as.character(all.peaks.interactions.L$V4), as.character(peak.TSS.L$PEAK))) -> peaks.any.interaction.L
# merge(M.R[,c(4,8,9)], PromLredSum, by.x="V8", by.y="V4") %>% dplyr::select(V4) -> all.peaks.interactions.R
# #unique(c(as.character(all.peaks.interactions.R$V4), as.character(peak.TSS.R$PEAK))) -> peaks.any.interaction.R
# unique(c(as.character(all.peaks.interactions.R$V4), as.character(all.peaks.interactions.L$V4))) -> peaks.with.interactions
# enh$V28[!enh$V28 %in% peaks.with.interactions] -> peaks.No.interact
# 

### Overlap peaks and loops

merge(M.L[,c(1:9,13,14)], PromRred, by.x="V13", by.y="V4") -> all.peaks.interactions.L.ALL
merge(M.R[,c(1:9,13,14)], PromLred, by.x="V13", by.y="V4") -> all.peaks.interactions.R.ALL

#merge( M.L, PromLredSum, by.x="V8", by.y="V4", all.x = TRUE ) -> ML_TSS
#merge( M.R, PromRredSum, by.x="V8", by.y="V4", all.x = TRUE ) -> MR_TSS

rbind(all.peaks.interactions.L.ALL, all.peaks.interactions.R.ALL) %>% arrange(V4)  -> loops
colnames(loops)[1]<-"ID"

loops$AdultOnly<-"NO"
loops$AdultOnly[which(loops$V5.y=="Yes" & loops$V6.y=="No" & loops$V7.y=="No" & loops$V8.y=="No")]<-"YES"
loops$FetalOnly<-"NO"
loops$FetalOnly[which(loops$V5.y=="No" & (loops$V6.y=="Yes" || loops$V7.y=="Yes" || loops$V8.y=="Yes"))]<-"YES"

### Annotated SNPs falling in enhancers

unique(sapply(strsplit(as.character(tss$V4), split=".", fixed=T),"[", 1)) -> ensGenes 
tss->converter
converter$Ens<-sapply(strsplit(as.character(tss$V4), split=".", fixed=T),"[", 1)
converter$Sym<-sapply(strsplit(as.character(tss$V4), split="|", fixed=T),"[", 2)

conv_ens_92<-read.table("../CLUMPS/All_Genes_Codes_Ensembl92.txt", header=T,sep="\t")

disease="BD"

load("/Users/Gabriel/Desktop/YALE/MASTER/HC/GERSHWIND_LOOPS/oTSS.rda")

ClosestTSS<-function(peak){
  filter(M, V4==peak) -> row
  filter(oTSS, as.character(V1) == row$V1) -> possibleTSS
  a<-which(row[c(7:11)] == "YES")
  which(rowSums(as.data.frame(possibleTSS[,c(6:10)][,a]) == "YES")>0) -> m
  bedTools.2closest(bed1= bedTools.2sort(bed1=row[,c(1:3)]),bed2=possibleTSS[m,],opt.string='-wb -d') ->cros
  return(cros)
}

for (disease in c("IQ_2018", "ADHD", "ASD", "CLOZUK", "HbA1c", "NEUROTICISM","HEIGHT","AD","BD","MDD","PD")) { 
  
  print(paste0("Annotating clups for ",disease, "........."))
  clumps<-read.table(paste0("../CLUMPS/Clumps/",disease,".clumps.bed"), header=TRUE)
  bedTools.2merge(bed1=bedTools.2sort(bed1=clumps)) -> clumps
  bedTools.2in(bed1=bedTools.2sort(bed1=clumps), bed2=bedTools.2sort(bed1=M), opt.string='-wb') %>% dplyr::select (V7,V8,V9) -> peaks_snps
  unique(peaks_snps$V7) -> peaks2eval
  
  write.table(bedTools.2sort(bed1=clumps), file=paste0(disease,'.bed'), quote=F, col.names=F, row.names=F)
  
  ### Interacting peaks
  filter(loops,V4 %in% peaks2eval) -> l
  ## Filter discordant peaks for ChipSeq
  #sapply(c(1:nrow(l)), function(x) { sum(which(l[x,c(6:10)] == "YES") %in% which(l[x,c(19:23)] == "YES"))}) -> concordant
  sapply(c(1:nrow(l)), function(x) { sum(which(l[x,c(6,7,9)] == "YES") %in% which(l[x,c(19,20,22)] == "YES"))}) -> concordantFet
  sapply(c(1:nrow(l)), function(x) { sum(which(l[x,c(8,10)] == "YES") %in% which(l[x,c(21,23)] == "YES"))}) -> concordantAd
  ## Filter discordant HiC
  #sapply(c(1:nrow(l)), function(x) { sum(l[x,c(6,7,9,19,20,22)] == "YES")}) -> preLoop
  #sapply(c(1:nrow(l)), function(x) { sum(l[x,c(8,10,21,23)] == "YES")}) -> postLoop
  mask<-c(which(l$FetalOnly =="YES" & concordantFet==0),which(l$AdultOnly=="YES" & concordantAd==0), which(concordantFet+concordantAd==0))
  if(length(mask)>0){
    genes<-unique(l[-mask,]$V13)
    peakInteract<-unique(l[-mask,]$V4)
  } else { 
    genes<-unique(l$V13)
    peakInteract<-unique(l$V4)
  }
  
  HiC_Ensmbl.HiC<-unique(sapply(strsplit(as.character(genes), split="|", fixed=T),"[", 1))
  HiC_Gsymbol.HiC<-unique(sapply(strsplit(as.character(genes), split="|", fixed=T),"[", 2))
  HiC.Ensmbl<-data.frame(gene=HiC_Ensmbl.HiC, case="HiC_Giusti")
  HiC.Gsymbol<-data.frame(gene=HiC_Gsymbol.HiC, case="HiC_Giusti")
  
  ### Non-interacting enhancer peaks
  peaks2eval[!peaks2eval %in% peakInteract] -> peaks2evalInNonInteract
  peaks2evalInNonInteract<-peaks2evalInNonInteract[peaks2evalInNonInteract %in% enh$V28]
  NiInt.Ensmbl<-c()
  NiInt.Gsymbol<-c()
  if(length(peaks2evalInNonInteract) > 0) {
    cross.all.table<-c()
    for (i in peaks2evalInNonInteract){
      #print(i)
      cross.all.table<-rbind(cross.all.table, ClosestTSS(i))
    }
    
    HiC_Ensmbl.NoInt<-unique(sapply(strsplit(as.character(cross.all.table$V7), split="|", fixed=T),"[", 1))
    HiC_Gsymbol.NoInt<-unique(sapply(strsplit(as.character(cross.all.table$V7), split="|", fixed=T),"[", 2))
    NiInt.Ensmbl<-data.frame(gene=HiC_Ensmbl.NoInt, case="NoInt_Giusti")
    NiInt.Gsymbol<-data.frame(gene=HiC_Gsymbol.NoInt, case="NoInt_Giusti")
  }  
  ### Peaks in promotor
  peaks2eval[peaks2eval %in% prom] -> peaks2evalInProm
  filter(Mhg38, V28 %in% peaks2evalInProm) %>% dplyr::select(allOverlapAnyENSEMBL,allOverlapAnyGSymbol) -> peaks2evalInPromGenes
  HiC_Ensmbl.Prom<-unique(unlist(strsplit(as.character(peaks2evalInPromGenes$allOverlapAnyENSEMBL), ',')))
  HiC_Gsymbol.Prom<-unique(unlist(strsplit(as.character(peaks2evalInPromGenes$allOverlapAnyGSymbol), ',')))
  Prom.Ensmbl<-data.frame(gene=HiC_Ensmbl.Prom, case="Promoter")
  Prom.Gsymbol<-data.frame(gene=HiC_Gsymbol.Prom, case="Promoter")
  
  ### Exonic variants
  s<-read.table(paste0("../CLUMPS/",disease,".toClump.output.clumped.rsid"), header=FALSE)$V1
  #filter(snp, Biotype=="protein_coding",Variant.consequence %in% c("missense_variant", "stop_gained", "splice_acceptor_variant","splice_region_variant","splice_donor_variant","3_prime_UTR_variant","5_prime_UTR_variant")) -> temp
  filter(snp, id %in% s, Variant.consequence %in% cat) -> temp
  
  Exo.Ensmbl<-c()
  Exo.Gsymbol<-c()
  if(dim(temp)[1]>0){
    unique(temp$ens) -> ge
    unique(filter(converter, Ens %in% ge)$V4) -> exonicVariants
    print(paste0(disease, "...." ,mean(ge %in% converter$Ens), "...", ge[!ge %in% converter$Ens]))
    #print(paste0(disease, "...." ,mean(exonicVariantsSym %in% converter$Sym), "...", exonicVariantsSym[!exonicVariantsSym %in% converter$Sym]))
    HiC_Ensmbl.Exo<-unique(sapply(strsplit(as.character(exonicVariants), split="|", fixed=T),"[", 1))
    HiC_Gsymbol.Exo<-unique(sapply(strsplit(as.character(exonicVariants), split="|", fixed=T),"[", 2))
    Exo.Ensmbl<-data.frame(gene=HiC_Ensmbl.Exo, case="Genic")
    Exo.Gsymbol<-data.frame(gene=HiC_Gsymbol.Exo, case="Genic")
  }
  
  HiC_Ensmbl<-ddply( rbind(Exo.Ensmbl,  Prom.Ensmbl, NiInt.Ensmbl, HiC.Ensmbl)  ,'gene',summarize,X=paste(case,collapse=",") )
  HiC_Gsymbol<-ddply( rbind(Exo.Gsymbol,  Prom.Gsymbol, NiInt.Gsymbol, HiC.Gsymbol)  ,'gene',summarize,X=paste(case,collapse=",") )
  
  write.table(HiC_Gsymbol, file=paste0("GENES/",disease,"_NC2.Gsymbol.2018.direction.PromOpen.E78.ClosestOpen.Concord.CORRECT2.txt"), quote=F, col.names=F, row.names=F)
  write.table(HiC_Ensmbl, file=paste0("GENES/",disease,"_NC2.Ensmbl.2018.direction.PromOpen.E78.ClosestOpen.Concord.CORRECT2.txt"), quote=F, col.names=F, row.names=F)
}

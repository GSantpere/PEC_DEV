

### METHYLATION
options("width"=300)

################## Load normalized data for METYL-ARRAY an Gencode38 annotation

load("all-drach-fresco-norm-beta-05272015.rda")
load("gencodeAnnot38.rda") 

### Read original phenotype data (includes info for Spiers paper)
pheno<-read.table("allsamplePheno.sorted.all.txt", header=T)

### Post conception days have been recalculated and some number needs to be updated in the pheno file
newPCD<-read.table("New_PCD.txt", header=T)

for (id in newPCD$Braincode){
  newPCD[which(newPCD$Braincode==id),"PCD"]
  pheno[which(pheno$Brain_ID==id),"age_pcd"] <- newPCD[which(newPCD$Braincode==id),"PCD"]		
}

### Make masks for fetal and post-natal groups

which(pheno$Brain_Region_Group=="Pre-Cortex" & pheno$age_period < 8) ->NCX.pre
which(pheno$Brain_Region_Group=="NCX" & pheno$age_period >= 8) -> NCX.post


### Produce linear model for each metylathion sites, separatedly in prenatal and postnatal

log(pheno[NCX.pre,"age_pcd"], 2)->logPre
t(apply(all.norm.beta[, NCX.pre], 1, function(x){coef(summary(lm(x~logPre + pheno[NCX.pre,"sex"])))["logPre",c(1,4)]})) -> lm.pre.sex.l
p.adjust(lm.pre.sex.l[,2], method="bonferroni") -> lm.pre.sex.l.ad

log(pheno[NCX.post,"age_pcd"], 2)->logPost
t(apply(all.norm.beta[, NCX.post], 1, function(x){coef(summary(lm(x~logPost + pheno[NCX.post,"sex"])))["logPost",c(1,4)]})) -> lm.post.sex.l
p.adjust(lm.post.sex.l[,2], method="bonferroni") -> lm.post.sex.l.ad

data.frame( DM.PRE=lm.pre.sex.l.ad, DM.PRE.noBonf=lm.pre.sex.l[,2], SIGN.PRE=lm.pre.sex.l[,1], 
            DM.POST=lm.post.sex.l.ad, DM.POST.noBonf=lm.post.sex.l[,2], SIGN.POST=lm.post.sex.l[,1] ) -> result.eff

max(logPre)->max.logPre
min(logPre)->min.logPre
apply(result.eff, 1, function(x){(max.logPre-min.logPre)*x[3]})->result.eff$Delta.PRE

max(logPost)->max.logPost
min(logPost)->min.logPost
apply(result.eff, 1, function(x){(max.logPost-min.logPost)*x[6]})->result.eff$Delta.POST

### Filter significant probes

pre.sig.pos.eff <- result.eff[which(result.eff$DM.PRE<=0.05 & result.eff$SIGN.PRE>0 & abs(result.eff$Delta.PRE) >= 0.1),]
pre.sig.neg.eff <- result.eff[which(result.eff$DM.PRE<=0.05 & result.eff$SIGN.PRE<0 & abs(result.eff$Delta.PRE) >= 0.1),]
post.sig.pos.eff <- result.eff[which(result.eff$DM.POST<=0.05 & result.eff$SIGN.POST>0 & abs(result.eff$Delta.POST) >= 0.1),]
post.sig.neg.eff <- result.eff[which(result.eff$DM.POST<=0.05 & result.eff$SIGN.POST<0 & abs(result.eff$Delta.POST) >= 0.1),] 
pre.nosig.eff <- result.eff[which(result.eff$DM.PRE.noBonf>0.05 & abs(result.eff$Delta.PRE)<0.1),]
post.nosig.eff <- result.eff[which(result.eff$DM.POST.noBonf>0.05 & abs(result.eff$Delta.POST)<0.1),]

#### ---> Create categories

list(pre.sig.pos.eff, pre.sig.neg.eff, post.sig.pos.eff, post.sig.neg.eff, pre.nosig.eff, post.nosig.eff) ->subsets
names(subsets)<-c("pre.UP","pre.DOWN","post.UP","post.DOWN","pre.EQUAL","post.EQUAL")
matrix(nrow=6, ncol=6)->mat2
colnames(mat2)=names(subsets)
rownames(mat2)<-colnames(mat2)
sub.list<-list()
for (i in names(subsets)){
  for (j in names(subsets)){
    mat2[i,j]<-length(which(rownames(as.data.frame(subsets[i])) %in% rownames(as.data.frame(subsets[j]))))
    rownames(as.data.frame(subsets[i]))[which(rownames(as.data.frame(subsets[i])) %in% rownames(as.data.frame(subsets[j])))]->id
    if(length(id)>0){
      sub.list[[paste(i,j,sep="_")]]<-id
    }
  }
}	

save(sub.list, file="sub.list.rda")


### Label order

ids<-c("pre.EQUAL_post.EQUAL","pre.UP_post.EQUAL", "pre.DOWN_post.EQUAL", "pre.EQUAL_post.UP", "pre.UP_post.UP", "pre.DOWN_post.UP", 
       "pre.DOWN_post.DOWN","pre.EQUAL_post.DOWN","pre.UP_post.DOWN")

### ---> Methylation level per category (local: use cpg as the list of Cpg)
met.levels.pre <- list()
met.levels.post <- list()
for (i in ids){
  print(i)
  temp.pre <- all.norm.beta[which(rownames(all.norm.beta[, NCX.pre]) %in% sub.list[[i]]), NCX.pre]
  temp.post <- all.norm.beta[which(rownames(all.norm.beta[, NCX.post]) %in% sub.list[[i]]), NCX.post]
  rowMeans(temp.pre) -> met.levels.pre[[i]]
  rowMeans(temp.post) -> met.levels.post[[i]]	
}

#### Plot the methylation levels
source("/Users/gabriel/Desktop/Yale/MET/DNAmethylation-Gabriel-09022015/vioplot2.R")
library(vioplot)
par(mfrow=c(3,4))
c<-1
for (i in ids){
  if(c %% 2 == 0){
    vioplot2(met.levels.pre[[i]], met.levels.post[[i]], col=c("#585758", "#585758"), ylim=c(0,1), names=c("Pre-Natal", "Post-Natal"))
    mtext(i,3) 
  }
  else {
    vioplot2(met.levels.pre[[i]], met.levels.post[[i]], col=c("#207935","#207935"), ylim=c(0,1), names=c("Pre-Natal", "Post-Natal")) 
    mtext(i,3)
  }
  c<-c+1
}

#### --> Plot CpG in categories (all probes, to run locally)
i=ids[-1]
par(mfrow=c(4,4))
c=1
### Prenatal probes
for (i in ids[-1]){
  if(c%%2==0){col="#207935"} else {col="#585758"}
  sample(sub.list[[i]], 200)->sam
  plot(logPre, rep(1,length(logPre)), pch=19, col="white", xlab="Period", ylab="beta", main="Pre-natal", ylim=c(0,1), cex=0.8)
  for (probe in sam){
    abline(lm(all.norm.beta[probe, NCX.pre]~ logPre), col= col, lwd=0.1)
  }
  ### PostNatal Probes
  plot(logPost, rep(1,length(logPost)), pch=19, col="white", xlab="Period", main="Post-natal", ylim=c(0,1), cex=0.8, ylab="", yaxt="n")
  for (probe in sam){
    abline(lm(all.norm.beta[probe, NCX.post]~ logPost), col= col, lwd=0.1)
  }
  c=c+1
}


### Create the list of genes associated to each category	

genes_GSymbl<-list()
genic_ENSG<-list()
for (i in ids){
  genes_GSymbl[[i]] <- unique(gencodeAnnot38[which(gencodeAnnot38$PeakID %in% sub.list[[i]] & gencodeAnnot38$GSymbl != "0"), "GSymbl"])
  a<-unique(gencodeAnnot38[which(gencodeAnnot38$PeakID %in% sub.list[[i]] & gencodeAnnot38$genic_ENSG != "None"), "genic_ENSG"])
  genic_ENSG[[i]] <- sapply(strsplit(a, "\\."), `[[`, 1) 
}	

### Get only common probes to calculate enrichments on gencode annotations
cpg.clean <- list()
for (i in ids){	
  cpg.clean[[i]]<-sub.list[[i]][which(sub.list[[i]] %in% gencodeAnnot38$PeakID)]
}

annot<-c()
counts<-c()
for (i in ids){
  a<-c(i,round(100*prop.table(table(gencode[which(gencodeAnnot38$PeakID %in% cpg.clean[[i]]), "genic_location"])),2))
  b<-c(i, table(gencodeAnnot38[which(gencodeAnnot38$PeakID %in% cpg.clean[[i]]), "genic_location"])) 
  rbind(annot,a)->annot
  rbind(counts,b)->counts	
}		


#### Cell specific probes
### METHYLATION --> Kozlenkov et al. 2014
load("CellgencodeAnnot-12152015.NEW.rda")
Cellannot[which(Cellannot$neur.mean>Cellannot$glia.mean),] -> NeuronHyper
Cellannot[which(Cellannot$neur.mean<Cellannot$glia.mean),] -> GliaHyper
cell<-c()
for (i in ids){
  NeuHyp <- length(which(NeuronHyper$PeakID %in% sub.list[[i]]))
  GliHyp <- length(which(GliaHyper$PeakID %in% sub.list[[i]]))
  vec<-c(NeuHyp, GliHyp)
  rbind(cell,c(i,vec))->cell	
}


### Newborn vs centenarian (local) --> Heyn et al. 2012
cent<-read.table("CpG.cent.txt", header=T, sep="\t") ## Careful it is hg19
cell.cent<-c()
HPO<-cent[which(cent$Direction.in.Y103=="Hypomethylated"),"Target.ID"]
HPER<-cent[which(cent$Direction.in.Y103=="Hypermethylated"),"Target.ID"]
for (i in ids){
  Hyper <- length(which(HPER %in% sub.list[[i]]))
  Hypo <- length(which(HPO %in% sub.list[[i]]))
  vec<-c(Hyper, Hypo)
  rbind(cell.cent,c(i,vec))->cell.cent	
}	

### Cell type specific expression
ct <- read.table("CellTypeList.BONA.txt", header=T, sep="\t")
#glia <- c("astrocyte", "oligodendrocyte", "microglia")
#neuron <- c("neuron", "interneuron", "pyramidalNeuron", "S1_pyramidalNeuron", "CA1_pyramidalNeuron", "striatalNeuron")
glia <- c("astrocyte", "oligodendrocyte")
neuron <- c("neuron", "interneuron", "pyramidalNeuron", "S1_pyramidalNeuron")

save(ct, glia, neuron, file='../BRAINSPAN/GeneLists/CellTypeSpecificRob.rda')

glia.gene <- unique(ct[which(ct$celltype %in% glia), "GeneID_human"])
neuron.gene <- unique(ct[which(ct$celltype %in% neuron), "GeneID_human"])

### Create MIXED categories
genic_ENSG[['AdultHIPO']]<-unique(unlist(genic_ENSG[c('pre.DOWN_post.DOWN', 'pre.EQUAL_post.DOWN', 'pre.UP_post.DOWN')]))
genic_ENSG[['AdultHIPER']]<-unique(unlist(genic_ENSG[c('pre.DOWN_post.UP', 'pre.EQUAL_post.UP', 'pre.UP_post.UP')]))

cell<-c()
for (i in c(ids,'AdultHIPER','AdultHIPO')){
  Neu <- length(which(neuron.gene %in% genic_ENSG[[i]]))
  Gli <- length(which(glia.gene %in% genic_ENSG[[i]]))
  vec<-c(Neu, Gli)
  rbind(cell,c(i,vec))->cell	
}


### ---> ENHANCERS (Changed to new enhancers)

## Bash
for i in $(cat LIST.txt); do intersectBed -a ${i}.cpg.bed -b ../../MASTER/BEDs/DFC.ENH.bed | wc -l; done
for i in $(cat LIST.txt); do intersectBed -a ${i}.cpg.bed -b ../../MASTER/BEDs/EqualDFC.ENH.bed | wc -l; done
for i in $(cat LIST.txt); do intersectBed -a ${i}.cpg.bed -b ../../MASTER/BEDs/DFC_Fetal.ENH.bed | wc -l; done
for i in $(cat LIST.txt); do intersectBed -a ${i}.cpg.bed -b ../../MASTER/BEDs/DFC_Adult.ENH.bed | wc -l; done


### ---> Non-CpG sites
#### Separate non-CpG sites
row.names(all.norm.beta)[-grep("cg", row.names(all.norm.beta))] -> nonCpG
non.CpG<-c()
for (i in ids){
  nonCpG.v <- length(which(nonCpG %in% sub.list[[i]]))
  rbind(non.CpG,c(i, nonCpG.v))-> non.CpG	
}	


all.norm.beta[which(rownames(all.norm.beta) %in% nonCpG), c(NCX.pre, NCX.post)] ->nocpg
apply(nocpg, 2, function(x) {length(which(x>0.5))} ) -> met

nonCpG[which(nonCpG %in% sub.list[["pre.EQUAL_post.UP"]])] -> nonCpG.inCat
nonCpG[which(nonCpG %in% Cellannot$PeakID)] -> nonCpG.cs
nonCpG.cs[which(nonCpG.cs %in% sub.list[["pre.EQUAL_post.UP"]])] -> nonCpG.inCat.cs
nonCpG[which(nonCpG %in% NeuronHyper$PeakID)] -> nonCpG.GUM
nonCpG[which(nonCpG %in% GliaHyper$PeakID)] -> nonCpG.NUM
nonCpG[-which(nonCpG %in% Cellannot$PeakID)] -> nonCpG.nocs

## Get a random sample of 2079 CpG
all.norm.beta[which(rownames(all.norm.beta) %in% sample(rownames(all.norm.beta), 2079)), c(NCX.pre, NCX.post)] ->random
apply(random, 2, function(x) {length(which(x>0.5))} ) -> met.r
data.frame(met=met.r, period=pheno[c(NCX.pre, NCX.post),"age_period"], var=rep("Random", 365)) -> met.r.p					
data.frame(met=met, period=pheno[c(NCX.pre, NCX.post),"age_period"], var=rep("NonCpG", 365)) -> met.p
rbind(met.p, met.r.p)->data

boxplot(100*met.r.p$met/2079~met.r.p$period, col="grey", ylim=c(0,60))
par(new=TRUE)
boxplot(100*met.p$met/2079~met.p$period,ylim=c(0,60),col=c(rep("tomato",7), rep("cadetblue",7)), xlab="Periods", ylab="% met C")

cbind(t(nocpg), pheno[c(NCX.pre, NCX.post),"age_period"])->nocpg.p
colnames(nocpg.p)[2080]<-"period"
melt(as.data.frame(nocpg.p), id=c("period"))->a
ggplot(a, aes(x=as.factor(period), y=value)) + scale_fill_manual(values=c(rep("tomato",7), rep("cadetblue",7))) + geom_violin()

data.frame(cpg=all.norm.beta["ch.1.101940785F", NCX.pre], age=logPre)->df
qplot(data=df, age, cpg, stat="identity",ordered=TRUE) + geom_smooth(method = "loess", size = 1.53) + theme_bw() 

par(mfrow=c(1,2))
plot(logPre, rep(1,length(logPre)), pch=19, col="white", xlab="Period", ylab="beta", main="Pre-natal", ylim=c(0,1), cex=0.8)
for (i in nonCpG.cs){
  lines(lowess(logPre, all.norm.beta[i,NCX.pre]), col="#207935", lwd=0.1)
}
for (i in nonCpG.inCat.cs){
  lines(lowess(logPre, all.norm.beta[i,NCX.pre]), col="darkred", lwd=0.3)
}
plot(logPost, rep(1,length(logPost)), pch=19, col="white", xlab="Period", ylab="beta", main="Pre-natal", ylim=c(0,1), cex=0.8)
for (i in nonCpG.cs){
  lines(lowess(logPost, all.norm.beta[i,NCX.post]), col="#207935", lwd=0.1)
}
for (i in nonCpG.inCat.cs){
  lines(lowess(logPost, all.norm.beta[i,NCX.post]), col="darkred", lwd=0.3)
}



plot(logPre, rep(1,length(logPre)), pch=19, col="white", xlab="Period", ylab="beta", main="Pre-natal", ylim=c(0,1), cex=0.8)
for (i in nonCpG.nocs){
  lines(lowess(logPre, all.norm.beta[i,NCX.pre]), col="#207935", lwd=0.1)
}
plot(logPost, rep(1,length(logPost)), pch=19, col="white", xlab="Period", ylab="beta", main="Pre-natal", ylim=c(0,1), cex=0.8)
for (i in nonCpG.nocs){
  lines(lowess(logPost, all.norm.beta[i,NCX.post]), col="#207935", lwd=0.1)
}
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
#### Enrichments and FIGURE


library(plotrix)
library(ggplot2)
library(scales)		

ids<-c("pre.EQUAL_post.EQUAL","pre.UP_post.EQUAL", "pre.DOWN_post.EQUAL", "pre.EQUAL_post.UP", "pre.UP_post.UP", "pre.DOWN_post.UP", 
       "pre.DOWN_post.DOWN","pre.EQUAL_post.DOWN","pre.UP_post.DOWN")

counts<-read.table("COUNTS.2017.txt", header=T, sep="\t")
counts[which(counts$Category=="pre.EQUAL_post.EQUAL"),]->all
final.or<-c()
final.p<-c()
for (id in ids){
  for (element in colnames(counts)[-c(1:4)]){
    
    counts[which(counts$Category==id),]->temp
    dimnames=list(c("Test","Ref"), c("Si","No"))
    
    TestYes<- as.numeric(temp[element])
    RefYes <- as.numeric(all[element])
    TestNo <- as.numeric(temp["Total"]) - as.numeric(temp[element])
    RefNo  <- as.numeric(all$Total)-as.numeric(all[element])
    
    if( element=="Neuronal" | element=="Glial" ){
      TestNo <- as.numeric(temp["Total.gene"]) - as.numeric(temp[element])
      RefNo  <- as.numeric(all$Total.gene)-as.numeric(all[element])		
    }
    
    if( element=="TSS" | element=="Exon" | element=="Intron" | element=="Intergenic" | element=="Enhancers" | element=="EnhFetal" | element=="EnhAdult" ){
      TestNo <- as.numeric(temp["TotalAnot"]) - as.numeric(temp[element])
      RefNo  <- as.numeric(all$TotalAnot)-as.numeric(all[element])		
    }
    
    
    fisher <- matrix(c(TestYes, RefYes, TestNo, RefNo), nrow=2, dimnames=dimnames)
    a <- fisher.test(fisher)
    or=(TestYes * RefNo) / (RefYes * TestNo)
    if (TestYes==0 | RefYes==0 | TestNo==0 | RefNo==0){
      or=((TestYes+0.5) * (RefNo+0.5)) / ((RefYes+0.5) * (TestNo+0.5)) #Haldane
    }	
    final.or <-rbind(final.or, c(id, element, round(or,2)))
    final.p  <-rbind(final.p, c(id, element, a$p.value))
  }
}

final.p<-final.p[-which(final.p[,1]=="pre.EQUAL_post.EQUAL"),]
final.or<-final.or[-which(final.or[,1]=="pre.EQUAL_post.EQUAL"),]

merge(as.data.frame(final.or, stringsAsFactors=FALSE), as.data.frame(final.p,stringsAsFactors=FALSE), by=c("V1","V2"))->df
colnames(df)<-c("Cat","Elem","OR","P")
df<-data.frame(Cat=df$Cat,Elem=df$Elem,OR=as.numeric(df$OR),P=as.numeric(df$P))  

df$logOR<-round(log2(df$OR),2)

levels = c("TSS","Exon","Intron","Intergenic", "Enhancers","EnhFetal","EnhAdult", "NUM","GUM","Neuronal","Glial")
with(df,factor(Elem,levels = levels))->df$Elem

df$Y1 <- cut(df$logOR,breaks = c(-Inf,-3,-1,-0.5,0,0.5,1,3,Inf),right = FALSE)

bonf<-0.05/dim(df)[1]
rep(0.2, dim(df)[1])->vec
vec[which(df$P<=bonf)]<-1
df$border=as.factor(vec)

with(df, factor(Cat, levels=rev(levels(df$Cat)[c(7,2,5,8,3,1,4,6)])))->df$Cat

df$col<-rep(1,dim(df)[1])
df$col[which(df$logOR<0)]<-2
df$col=as.factor(df$col)
save(df, file='Integration.table.rda')

myColors=c("#C10834", "#334A96")
names(myColors) <- levels(df$col)
ggplot(df, aes(x=Elem, y=Cat, size=abs(logOR), fill=col, alpha=border),guide=FALSE)+
  geom_point(shape=21) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "white")) + scale_fill_manual(name = "col",values = myColors) + scale_size_continuous(range = c(2, 30)) + scale_alpha_discrete(range=c(0.5,1))


#### Simplified version


counts<-read.table("COUNTS.simple.2017.txt", header=T, sep="\t")
counts[which(counts$Category=="pre.EQUAL_post.EQUAL"),]->all
final.or<-c()
final.p<-c()
for (id in counts$Category){
  for (element in colnames(counts)[-c(1:4)]){
    
    counts[which(counts$Category==id),]->temp
    dimnames=list(c("Test","Ref"), c("Si","No"))
    
    TestYes<- as.numeric(temp[element])
    RefYes <- as.numeric(all[element])
    TestNo <- as.numeric(temp["Total"]) - as.numeric(temp[element])
    RefNo  <- as.numeric(all$Total)-as.numeric(all[element])
    
    if( element=="Neuronal" | element=="Glial" ){
      TestNo <- as.numeric(temp["Total.gene"]) - as.numeric(temp[element])
      RefNo  <- as.numeric(all$Total.gene)-as.numeric(all[element])		
    }
    
    if( element=="TSS" | element=="Exon" | element=="Intron" | element=="Intergenic" | element=="Enhancers" | element=="EnhFetal" | element=="EnhAdult" ){
      TestNo <- as.numeric(temp["TotalAnot"]) - as.numeric(temp[element])
      RefNo  <- as.numeric(all$TotalAnot)-as.numeric(all[element])		
    }
    
    
    fisher <- matrix(c(TestYes, RefYes, TestNo, RefNo), nrow=2, dimnames=dimnames)
    a <- fisher.test(fisher)
    or=(TestYes * RefNo) / (RefYes * TestNo)
    if (TestYes==0 | RefYes==0 | TestNo==0 | RefNo==0){
      or=((TestYes+0.5) * (RefNo+0.5)) / ((RefYes+0.5) * (TestNo+0.5)) #Haldane
    }	
    final.or <-rbind(final.or, c(id, element, round(or,2)))
    final.p  <-rbind(final.p, c(id, element, a$p.value))
  }
}

final.p<-final.p[-which(final.p[,1]=="pre.EQUAL_post.EQUAL"),]
final.or<-final.or[-which(final.or[,1]=="pre.EQUAL_post.EQUAL"),]

merge(as.data.frame(final.or, stringsAsFactors=FALSE), as.data.frame(final.p,stringsAsFactors=FALSE), by=c("V1","V2"))->df
colnames(df)<-c("Cat","Elem","OR","P")
df<-data.frame(Cat=df$Cat,Elem=df$Elem,OR=as.numeric(df$OR),P=as.numeric(df$P))  

df$logOR<-round(log2(df$OR),2)

levels = c("TSS","Enhancers","EnhFetal","EnhAdult", "GUM", "NUM","Neuronal","Glial")
with(df,factor(Elem,levels = levels))->df$Elem

df$Y1 <- cut(df$logOR,breaks = c(-Inf,-3,-1,-0.5,0,0.5,1,3,Inf),right = FALSE)

bonf<-0.05/dim(df)[1]
rep(0.2, dim(df)[1])->vec
vec[which(df$P<=bonf)]<-1
df$border=as.factor(vec)

with(df, factor(Cat, levels=rev(levels(df$Cat)[c(4,3,2,1)])))->df$Cat

df$col<-rep(1,dim(df)[1])
df$col[which(df$logOR<0)]<-2
df$col=as.factor(df$col)
save(df, file='Integration.table.simple.rda')

myColors=c("#C10834", "#334A96")
names(myColors) <- levels(df$col)
ggplot(df, aes(x=Elem, y=Cat, size=abs(logOR), fill=col, alpha=border),guide=FALSE)+
  geom_point(shape=21) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "white")) + scale_fill_manual(name = "col",values = myColors) + scale_size_continuous(range = c(2, 30)) + scale_alpha_discrete(range=c(0.5,1))


#### ENHANCERS


### Load enhacers from DeSeq2 analysis

#load('/Users/Gabriel/Desktop/YALE/MASTER/definitive.ENH.ANYPEAK.rda')
load('/Users/Gabriel/Desktop/YALE/MASTER/BEDs/LIST_BRAINSPAN/definitive.ENH.BRAINSPAN.rda')
definitive.ENH2-> definitive.ENH
#load('/Users/Gabriel/Desktop/YALE/MASTER/definitive.ANYPEAK.rda')  ## I do all peaks
definitive-> definitive.ENH
load('/Users/Gabriel/Desktop/YALE/MASTER/BEDs/LIST_BRAINSPAN/definitive.BRAINSPAN.rda')
definitive.2-> definitive.ENH

definitive.ENH[['novarDFC']]<-definitive.ENH[['DFC']][-which( rownames(definitive.ENH[['DFC']]) %in% c(rownames(definitive.ENH[['DFC_Fetal']]),rownames(definitive.ENH[['DFC_Adult']]))),]
#definitive[['novarDFC']]<-definitive[['DFC']][-which( rownames(definitive[['DFC']]) %in% c(rownames(definitive[['DFC_Fetal']]),rownames(definitive[['DFC_Adult']]))),]


genes.ALL<-list() 	
for(i in c('DFC_Fetal','DFC_Adult','EqualDFC','novarDFC','UnvarDFC','DFC')){
  genes.ALL[[i]]<-unique(unlist(strsplit(as.character(definitive.ENH[[i]][,"allOverlapAnyGSymbol"]),split=',')))
}

read.table("dexDESeq2_developWindow/W4-vs-W9.DFC.all.DEX.xls", header=T)->file

file[which(file$geneType=="protein_coding"),1]->PC
sapply(strsplit(as.character(PC), split="|", fixed=T),"[", 2) ->PC.names
file[which(file$geneType!="protein_coding"),1]->NC
sapply(strsplit(as.character(NC), split="|", fixed=T),"[", 2) ->NC.names

genes.ALL[['DFC']][which(genes.ALL[['DFC']] %in% NC.names)] -> DFC.enh.genes.NC.exp
genes.ALL[['DFC']][which(genes.ALL[['DFC']] %in% PC.names)] -> DFC.enh.genes.PC.exp

file$Gene<-sapply(strsplit(as.character(file$Geneid), split="|", fixed=T),"[", 2) 
file<-file[which(file$baseMeanA.RPKM. >=1 | file$baseMeanB.RPKM. >=1),]

genes.ALL[['ALL']]<-sapply(strsplit(as.character(file[,1]), split="|", fixed=T),"[", 2) ## Genes expressed in DFC

#### Expression section

### Correlated genes and patterns for each subset of enhancer-gene pair

#### Match expression pattern and enhancer pattern
### This can be done with all genes annotated or only PC

mylist=genes.ALL ## I chose this one

M<-matrix(ncol=3, nrow=5)
rownames(M) <- c('DFC_Fetal','DFC_Adult','EqualDFC','novarDFC','UnvarDFC')
colnames(M) <- c("num_genes","Fetal","Adult")

for(i in c(1:5)){ 			
  for(j in c(1,2,3)){
    if(j==1){				
      temp<-file[which(file$Gene %in% mylist[[rownames(M)[i]]]),]  ## Change this lie 					
      M[i,1]<-dim(temp)[1]
    }	
    if(j==2){ 
      M[i,j]<-length(temp[which(temp$log2FoldChange<=-1 & temp$padj<0.01),'Gene'])
      write.table( file=paste('GENE_LISTS/', rownames(M)[i],'_',colnames(M)[j], '.txt', sep=''), temp[which(temp$log2FoldChange<=-1 & temp$padj<0.01),'Geneid'], quote=F, row.names=F, col.names=F)
    }
    if(j==3){ 
      M[i,j]<-length(temp[which(temp$log2FoldChange>=1 & temp$padj<0.01),'Gene'])
      write.table( file=paste('GENE_LISTS/', rownames(M)[i],'_',colnames(M)[j], '.txt', sep=''), temp[which(temp$log2FoldChange>=1 & temp$padj<0.01),'Geneid'], quote=F, row.names=F, col.names=F)  
    }
  }
  i<-i+1
}

# listU<-list()
# temp<-file[which(file$Gene %in% mylist[['DFC_Fetal']]),]
# listU[['DFC_Fetal']]<-temp[which(temp$log2FoldChange<=-1 & temp$padj<0.01),'Gene']
# temp<-file[which(file$Gene %in% mylist[['DFC_Adult']]),] 
# listU[['DFC_Adult']]<-temp[which(temp$log2FoldChange>=1 & temp$padj<0.01),'Gene']
# listU[['ALL']]<-file[which(file$Gene %in% mylist[['DFC']]),'Gene']

ct <- read.table("~/Desktop/Yale/INTEGRATION/CellTypeList.BONA.txt", header=T, sep="\t")
#glia <- c("astrocyte", "oligodendrocyte", "microglia")
#neuron <- c("neuron", "interneuron", "pyramidalNeuron", "S1_pyramidalNeuron", "CA1_pyramidalNeuron", "striatalNeuron")
glia <- c("astrocyte", "oligodendrocyte")
neuron <- c("neuron", "interneuron", "pyramidalNeuron", "S1_pyramidalNeuron")

glia.gene <- unique(ct[which(ct$celltype %in% glia), "GeneName"])
neuron.gene <- unique(ct[which(ct$celltype %in% neuron), "GeneName"])

cell<-c()
for (i in names(mylist)){
  Neu <- length(which(neuron.gene %in% mylist[[i]]))
  Gli <- length(which(glia.gene %in% mylist[[i]]))
  vec<-c(Neu, Gli)
  rbind(cell,c(i,vec))->cell	
} 
colnames(cell)<-c("pat", "Neuron", "Glia")

cbind(M,cell[-c(6,7),])[,-4]->M2

### GUMs and NUMs
source('/Users/Gabriel/Desktop/bin/BEDtoolsR.R')
load('/Users/Gabriel/Desktop/YALE/INTEGRATION/CellgencodeAnnot-12152015.NEW.rda')
load('/Users/Gabriel/Desktop/YALE/INTEGRATION/gencodeAnnot38.rda')
load('/Users/Gabriel/Desktop/YALE/INTEGRATION/sub.list.rda')

Cellannot[which(Cellannot$neur.mean>Cellannot$glia.mean),] -> GUM
Cellannot[which(Cellannot$neur.mean<Cellannot$glia.mean),] -> NUM

bedTools.2sort(bed1=GUM[,c(5,6,7)])->gum
bedTools.2sort(bed1=NUM[,c(5,6,7)])->num

bedTools.2sort(bed1=gencodeAnnot38[,c(2,3,4,1)])->array

t<-c()
for(i in c('DFC_Fetal','DFC_Adult','EqualDFC','novarDFC','UnvarDFC')){
  bedTools.2sort(bed1=definitive.ENH[[i]])->temp
  bedTools.2in(bed1=gum, bed2=temp)->tempg
  bedTools.2in(bed1=num, bed2=temp)->tempn
  bedTools.2in(bed1=array, bed2=temp)->tempa
  
  bedTools.2in(bed1=array[which(array$V4 %in% unlist(sub.list[c('pre.DOWN_post.DOWN', 'pre.EQUAL_post.DOWN', 'pre.UP_post.DOWN')])),], bed2=temp)->tempPostHypo
  bedTools.2in(bed1=array[which(array$V4 %in% unlist(sub.list[c('pre.DOWN_post.UP', 'pre.EQUAL_post.UP', 'pre.UP_post.UP')])),], bed2=temp)->tempPostHyper
  
  rbind(t, c(i, dim(tempa)[1], dim(tempg)[1], dim(tempn)[1], dim(tempPostHypo)[1], dim(tempPostHyper)[1]))->t
}

colnames(t)<-c('cat','num_cpgs','GUM','NUM', 'PostHypo', 'PostHyper')
t<-as.data.frame(t)
rownames(t)<-t$cat
t[,-1]->t
t<-apply(t,2, as.numeric)
rownames(t)<-c('DFC_Fetal','DFC_Adult','EqualDFC','novarDFC','UnvarDFC')


M2<-apply(M2, 2, as.numeric)
rownames(M2)<-rownames(M)



##### Test the enrichment and plot the results

### 1
dimnames<-list(rownames(M2)[c(1,2)], colnames(M2)[-1])

Mor<-matrix(ncol=4, nrow=2, dimnames=dimnames)
Mp <-matrix(ncol=4, nrow=2 ,dimnames=dimnames)

for (en.p in rownames(M2)[c(1,2)]){	
  for (ex.p in colnames(M2)[-1]){
    
    TestYes<- M2[en.p,ex.p]
    RefYes <- M2["novarDFC", ex.p]
    TestNo <- M2[en.p, 1] - M2[en.p,ex.p]
    RefNo  <- M2["novarDFC", 1]-M2["novarDFC", ex.p]
    
    fisher <- matrix(c(TestYes, RefYes, TestNo, RefNo), nrow=2)
    p <- fisher.test(fisher)
    or=(TestYes * RefNo) / (RefYes * TestNo)
    if (TestYes==0 | RefYes==0 | TestNo==0 | RefNo==0){
      or=((TestYes+0.5) * (RefNo+0.5)) / ((RefYes+0.5) * (TestNo+0.5)) #Haldane
    }
    Mor[en.p,ex.p]<-or
    Mp[en.p,ex.p]<-p$p.value		
  }			
}

df<-cbind(melt(Mor), melt(Mp))[,-c(4,5)]
colnames(df)<-c("EnhPat", "ExpPat", "OR", "P")


### 2
dimnames<-list(rownames(t)[c(1,2)], colnames(t)[-1])

Mor<-matrix(ncol=4, nrow=2, dimnames=dimnames)
Mp <-matrix(ncol=4, nrow=2 ,dimnames=dimnames)

for (en.p in rownames(t)[c(1,2)]){	
  for (ex.p in colnames(t)[-1]){
    
    TestYes<- t[en.p,ex.p]
    RefYes <- t["novarDFC", ex.p]
    TestNo <- t[en.p, 1] - t[en.p,ex.p]
    RefNo  <- t["novarDFC", 1]-t["novarDFC", ex.p]
    
    fisher <- matrix(c(TestYes, RefYes, TestNo, RefNo), nrow=2)
    p <- fisher.test(fisher)
    or=(TestYes * RefNo) / (RefYes * TestNo)
    if (TestYes==0 | RefYes==0 | TestNo==0 | RefNo==0){
      or=((TestYes+0.5) * (RefNo+0.5)) / ((RefYes+0.5) * (TestNo+0.5)) #Haldane
    }
    Mor[en.p,ex.p]<-or
    Mp[en.p,ex.p]<-p$p.value		
  }			
}

df2<-cbind(melt(Mor), melt(Mp))[,-c(4,5)]
colnames(df2)<-c("EnhPat", "ExpPat", "OR", "P")

df<-rbind(df, df2)

bonf<-0.05/dim(df)[1]	
#df$adjustedP<-p.adjust(df$P, "BH")

df$col<-rep(1,dim(df)[1])
df$col[which(df$P<0.05)]<-2
df$direction<-rep(1,dim(df)[1])
df$direction[which(df$OR<1)]<-2
df$direction=as.factor(df$direction)
#df$OR[which(df$OR<1)]<-1/df$OR[which(df$OR<1)]
with(df,factor(EnhPat,levels = rev(c("DFC_Fetal","DFC_Adult"))))->df$EnhPat
with(df,factor(ExpPat,levels = c("Fetal","Adult","PostHypo","PostHyper","Neuron", "Glia", 'GUM', 'NUM')))->df$ExpPat

myColors=c("#C10834", "#334A96", "white")
names(myColors) <- levels(df$col)
ggplot(df, aes(x=ExpPat, y= EnhPat, size=abs(log2(OR)), fill=direction, alpha=col),guide=FALSE)+
  geom_point(shape=21) + theme_bw() + theme(axis.line = element_line(colour = "black")) + scale_fill_manual(name = "direction",values = myColors) + scale_size_continuous(range = c(2, 30)) + scale_alpha(range = c(0.5, 1))   


#############################	
#### GO categories ##########     
############################# 

library(org.Hs.eg.db)
library(topGO)

mylist= genes.ALL
number= 5018  # For DFC

compare.topGO.DFC<-list()
totals.DFC<-c()
for (i in c("DFC_Fetal","DFC_Adult")){	    
  j="ALL"           
  if (i == j){ next }
  myInterestingGenes<-mylist[[i]]
  gNames<-unique(mylist[[j]])
  
  genemap <- select(org.Hs.eg.db, myInterestingGenes, "ENTREZID", "SYMBOL")
  genemap <- genemap[!is.na(genemap$ENTREZID),]
  univmap <- select(org.Hs.eg.db, gNames, "ENTREZID", "SYMBOL") 
  univmap <- univmap[!is.na(univmap$ENTREZID),]
  
  geneList <- factor(as.integer(univmap$ENTREZID %in% genemap$ENTREZID))
  names(geneList)<-univmap$ENTREZID
  
  GOdata <- new("topGOdata", description="Brainspan", ontology="BP", allGene=geneList, nodeSize=10, 
                annotationFun=annFUN.org, mapping="org.Hs.eg.db")
  
  totals.DFC <-rbind(totals.DFC,c(i,length(which(feasible(GOdata))), length(which(attributes(GOdata)$allScores[feasible(GOdata)]==1))))
  
  resultFisher <- runTest(GOdata,algorithm="classic",statistic="fisher")
  pValue.classic <- score(resultFisher)
  n<-i
  compare.topGO.DFC[[n]] <- GenTable(GOdata, classicFisher = resultFisher, orderBy = "classic", ranksOf = "classicFisher", topNodes = number) 
  write.table(compare.topGO.DFC[[n]], file=paste(n,".topGO.DFC.AllAnnots.WholeList.txt",sep=""), quote=FALSE, row.names=FALSE, sep="\t")
  
}	

selected<-read.table(file="SelectedGoterms.txt", header=F, sep="\t", stringsAsFactors = FALSE)
selected[which(selected$V1=='GO:0007268'),'V2']<-'chemical synaptic transmission'
makeTestToAllTopGO(compare.topGO.DFC, subset="TRUE", terms.subset=as.character(selected$V1), totals=totals.DFC, region="DFC") -> t.or.DFC

#### Whole figure

part1<-df[,-4]
part1$OR<-log2(part1$OR)
colnames(part1)<-c('Cat','Term','logOR','Sig','Direction')

t.or.DFC[,c(1,4,5,12,13)] ->part2
colnames(part2)<-c('Cat','Term','logOR','Sig','Direction')

wholeset<-rbind( part1, part2 ) 
wholeset$Term <- with(wholeset, factor(Term, levels=c('Fetal', 'Adult', 'PostHypo','PostHyper','Neuron', 'Glia', 'NUM', 'GUM',as.character(selected$V2)))) 

save(wholeset, file='Enhancer.Table2plotMet.July.2017.rda')
save(wholeset, file='Enhancer.Table2plotMet.July.2017.ONLY_ENH.rda')

myColors=c("#C10834", "#334A96")
ggplot(wholeset, aes(x=Term, y=Cat, size=abs(logOR), fill=factor(Direction), alpha=Sig),guide=FALSE)+
  geom_point(shape=21) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"), axis.text.x = element_text(angle = 45, hjust = 1)) + scale_alpha_discrete(range = c(0.5, 1)) + scale_fill_manual(values = myColors) + scale_size_continuous(range = c(1, 20)) + scale_colour_manual(values = "black") 

#### Check expression

rpkm<-read.table("brainSpan.hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.RPKM.normalized.CQNCombat.txt", header=T)
rownames(rpkm)<-rpkm$Geneid
meta<-read.table("METADATA.rpkm.txt", header=T)
meta$ID<-paste(meta$Braincode, meta$Regioncode, sep=".")

grep('DFC',meta[which(meta$Period==6),'ID'], value=TRUE) ->w4
grep('DFC',meta[which(meta$Period==13),'ID'], value=TRUE) ->w9
rpkm[,which(colnames(rpkm) %in% c(w4,w9))]->sub

#### FUNCTIONS

makeTestToAllTopGO<-function(listU, subset, terms.subset, totals, region){
  
  if(subset=="TRUE"){
    compare.topGO.terms<-list()
    for (i in names(listU)){
      compare.topGO.terms[[i]]<-listU[[i]][grep(paste(terms.subset,collapse="|"), listU[[i]][,1]),]
    }
  }
  
  else{compare.topGO.terms<-listU}
  
  totals<-as.data.frame(totals)
  totals<-data.frame(V1=as.character(totals$V1), V2=as.numeric(as.character(totals$V2)), V3=as.numeric(as.character(totals$V3)))
  
  table.or<-list()
  t.or<-c()
  for (i in names(compare.topGO.terms)){
    
    temp<-totals[which(totals$V1==i),]
    counts<-compare.topGO.terms[[i]][,4]
    size<-compare.topGO.terms[[i]][,3]
    exp<-compare.topGO.terms[[i]][,5]
    P<-compare.topGO.terms[[i]][,6]
    
    TestYes<-counts
    TestNo<- as.numeric(temp$V3) - counts
    RefYes<- size - TestYes
    RefNo<- as.numeric(temp$V2) - TestYes - TestNo - RefYes
    
    NoHaldane<-cbind(TestYes, TestNo, RefYes, RefNo)
    
    which(TestYes==0 | RefYes==0 | TestNo==0 | RefNo==0)->zeros
    sum=rep(0, length(TestYes))
    sum[zeros]<-0.5
    TestYes<-TestYes+sum
    TestNo<-TestNo+sum
    RefYes<-RefYes+sum
    RefNo<-RefNo+sum
    
    or=(TestYes * RefNo) / (RefYes * TestNo) #Haldane
    
    table.or[[i]]<-cbind(compare.topGO.terms[[i]][,c(1:2)], round(log2(or),2), NoHaldane)
    t.or<-rbind(t.or, cbind(rep(i, length(or)), P, table.or[[i]]))
  }
  
  colnames(t.or)<-c("Cat","P","GOID","Term","OR", "TestYes", "TestNo", "RefYes", "RefNo")
  
  dimnames=list(c("Test","Ref"), c("Si","No"))
  a.g<-c()
  a.l<-c()
  for(c in c(1:length(t.or$Cat))){
    fisher <- matrix(c(t.or$TestYes[c], t.or$RefYes[c], t.or$TestNo[c], t.or$RefNo[c]), nrow=2, dimnames=dimnames)
    a.g <- c(a.g, fisher.test(fisher, alternative='greater')$p.value)
    a.l <- c(a.l, fisher.test(fisher, alternative='less')$p.value)
  }
  ### When doing our Fisher we obtain same Pvalues than GOstats
  t.or$Pgreater=a.g
  t.or$Pless=a.l
  t.or$P<-as.numeric(as.character(t.or$P))
  t.or$P[t.or$OR<0]<-t.or$Pless[t.or$OR<0] 
  
  t.or$col<-rep(1,dim(t.or)[1])
  as.numeric(as.character(t.or$P))->t.or$P
  t.or$col[which(t.or$P<0.05)]<-2
  t.or$col=as.factor(t.or$col)
  t.or$direction<-rep(1,dim(t.or)[1])
  t.or$direction[which(2^t.or$OR<1)]<-2
  
  return(t.or)   
}  

#### FINAL INTEGRATED ANALYSIS

oad("EnsembleGenesBelenMarch2016.rda")  ## These are DMR identified genes
load('all-drach-fresco-norm-beta-05272015.rda')
names(genic_ENSG)
[1] "All"         "PreHy"       "PreHp"       "PostHy"      "PostHp"     
[6] "PreHyPostHy" "PreHyPostHp" "PreHpPostHy" "PreHpPostHp" "PreHyPostEq"
[11] "PreHpPostEq" "PreEqPostHy" "PreEqPostHp"

names(genic_ENSG)<-c("All","PreHyper","PreHypo","PostHyper","PostHypo","PreHyper_PostHyper","PreHyper_PostHypo", "PreHypo_PostHyper",
                     "PreHypo_PostHypo", "PreHyper_PostEqual", "PreHypo_PostEqual", "PreEqual_PostHyper", "PreEqual_PostHypo")

## Correct ALL category to include all genes targeted in the array
load("gencodeAnnot38.rda")
gencodeAnnot38[which(gencodeAnnot38$PeakID %in% rownames(all.norm.beta)),"genic_ENSG"]->all.genes
sapply(strsplit(all.genes, "\\."), `[[`, 1) -> all.genes.clean
unique(all.genes.clean)[-1]-> all.genes.clean
genic_ENSG[["All"]] <- all.genes.clean  

## For the simplified version  
genic_ENSG[['pre.UP_post.EQUAL']]<-genic_ENSG[['PreHyper_PostEqual']]
genic_ENSG[['pre.DOWN_post.EQUAL']]<-genic_ENSG[['PreHypo_PostEqual']]
genic_ENSG[['PostUP']]<-unique(c(genic_ENSG[['PreHyper_PostHyper']],genic_ENSG[['PreHypo_PostHyper']],genic_ENSG[['PreEqual_PostHyper']]))
genic_ENSG[['PostDOWN']]<-unique(c(genic_ENSG[['PreHypo_PostHypo']],genic_ENSG[['PreHyper_PostHypo']],genic_ENSG[['PreEqual_PostHypo']]))

### TopGO againts all

number= 6075
compare.topGO<-list()
totals<-c()
for (i in c("PreHyper_PostHyper","PreHyper_PostHypo", "PreHypo_PostHyper",
            "PreHypo_PostHypo", "PreHyper_PostEqual", "PreHypo_PostEqual", "PreEqual_PostHyper", "PreEqual_PostHypo", 'PostUP', 'PostDOWN')){
  j="All"           
  if (i == j){ next }
  myInterestingGenes<-genic_ENSG[[i]]
  gNames<-unique(genic_ENSG[[j]])
  
  genemap <- AnnotationDbi::select(org.Hs.eg.db, myInterestingGenes, "ENTREZID", "ENSEMBL")
  genemap <- genemap[-which(is.na(genemap$ENTREZID)),]
  univmap <- AnnotationDbi::select(org.Hs.eg.db, gNames, "ENTREZID", "ENSEMBL") 
  univmap <- univmap[-which(is.na(univmap$ENTREZID)),]
  
  geneList <- factor(as.integer(univmap$ENTREZID %in% genemap$ENTREZID))
  names(geneList)<-univmap$ENTREZID
  
  GOdata <- new("topGOdata", description="MiRNA rebuttal", ontology="BP", allGene=geneList, nodeSize=10, 
                annotationFun=annFUN.org, mapping="org.Hs.eg.db")
  
  
  totals<-rbind(totals,c(i,length(which(feasible(GOdata))), length(which(attributes(GOdata)$allScores[feasible(GOdata)]==1))))
  
  resultFisher <- runTest(GOdata,algorithm="classic",statistic="fisher")
  pValue.classic <- score(resultFisher)
  n<-i
  #compare[[n]] <- GenTable(GOdata, classicFisher = resultFisher, orderBy = "classic", ranksOf = "classicFisher", topNodes = l)
  compare.topGO[[n]] <- GenTable(GOdata, classicFisher = resultFisher, orderBy = "classic", ranksOf = "classicFisher", topNodes = number)
  compare.topGO[[n]]$Adjust<-p.adjust(compare.topGO[[n]]$classicFisher, method="BH") 
  write.table(compare.topGO[[n]], file=paste('GOJuliol/',n,".topGO.all.txt",sep=""), quote=FALSE, row.names=FALSE, sep="\t")
  
}	

data.frame(V1=totals[,1], V2=as.numeric(totals[,2]), V3=as.numeric(totals[,3]))->totals

table.or<-list()
t.or<-c()
for (i in names(compare.topGO)){
  
  temp<-totals[which(totals$V1==i),]
  counts<-compare.topGO[[i]][,4]
  size<-compare.topGO[[i]][,3]
  exp<-compare.topGO[[i]][,5]
  P<-compare.topGO[[i]][,6]
  
  TestYes<-counts
  TestNo<- as.numeric(temp$V3) - counts
  RefYes<- size - TestYes
  RefNo<- as.numeric(temp$V2) - TestYes - TestNo - RefYes
  
  NoHaldane<-cbind(TestYes, TestNo, RefYes, RefNo)
  
  which(TestYes==0 | RefYes==0 | TestNo==0 | RefNo==0)->zeros
  sum=rep(0, length(TestYes))
  sum[zeros]<-0.5
  TestYes<-TestYes+sum
  TestNo<-TestNo+sum
  RefYes<-RefYes+sum
  RefNo<-RefNo+sum
  
  or=(TestYes * RefNo) / (RefYes * TestNo) #Haldane
  
  table.or[[i]]<-cbind(compare.topGO[[i]][,c(1:2)], round(log2(or),2), NoHaldane)
  t.or<-rbind(t.or, cbind(rep(i, length(or)), P, table.or[[i]]))
}

colnames(t.or)<-c("Cat","P","GOID","Term","OR", "TestYes", "TestNo", "RefYes", "RefNo")

dimnames=list(c("Test","Ref"), c("Si","No"))
a.g<-c()
a.l<-c()
for(c in c(1:length(t.or$Cat))){
  fisher <- matrix(c(t.or$TestYes[c], t.or$RefYes[c], t.or$TestNo[c], t.or$RefNo[c]), nrow=2, dimnames=dimnames)
  a.g <- c(a.g, fisher.test(fisher, alternative='greater')$p.value)
  a.l <- c(a.l, fisher.test(fisher, alternative='less')$p.value)
}
### When doing our Fisher we obtain same Pvalues than GOstats
t.or$Pgreater=a.g
t.or$Pless=a.l

t.or$col<-rep(1,dim(t.or)[1])
as.numeric(as.character(t.or$P))->t.or$P
t.or$col[which(t.or$P<0.01)]<-2
t.or$col=as.factor(t.or$col)
t.or$direction<-rep(1,dim(t.or)[1])
t.or$direction[which(2^t.or$OR<1)]<-2


### Compare list with the one from enhancers to obtain GO categories to compare and plot

for (i in c("PreHyper_PostHyper","PreHyper_PostHypo", "PreHypo_PostHyper",
            "PreHypo_PostHypo", "PreHyper_PostEqual", "PreHypo_PostEqual", "PreEqual_PostHyper", "PreEqual_PostHypo")){
  temp<-compare.topGO[[i]]
  temp$Adjust<-p.adjust(temp$classicFisher, method="BH")
  write.table(temp[1:20,], file=paste("GO/TopGO/",i,".topGO.top20.txt",sep=""), quote=FALSE, row.names=FALSE, sep="\t")
}

for (i in c("Fetal_and_FetalInf.CBC.topGO.all.DFC.fused.txt","PostNatal_and_Adult.CBC.topGO.all.DFC.fused.txt",
            "Fetal_and_FetalInf.DFC.topGO.all.DFC.fused.txt","PostNatal_and_Adult.DFC.topGO.all.DFC.fused.txt")){
  temp<-read.table(paste("/GO_fused/",i,sep=""), header=T, sep="\t")
  temp$Adjust<-p.adjust(temp$classicFisher, method="BH")
  write.table(temp[1:20,], file=paste("GO/TopGO/",i,".topGO.top20.txt",sep=""), quote=FALSE, row.names=FALSE, sep="\t")  
}

#### Selected categories from the comparison of the first 20 top categories

selected<-read.table(file="GO/SelectedGoterms.txt", header=F, sep="\t", stringsAsFactors = FALSE)
colnames(selected)<-c("GOID", "Term")
selected[which(selected$GOID=='GO:0007268'),'Term']<-'chemical synaptic transmission'
t.or[which(t.or$GOID %in% selected$GOID),]->t.or.sel

with(t.or.sel, factor(Cat, levels=rev(c("PreHyper_PostEqual", "PreHypo_PostEqual", "PostUP", "PostDOWN", "PreEqual_PostHyper", "PreHyper_PostHyper", "PreHypo_PostHyper", "PreHypo_PostHypo", "PreEqual_PostHypo", "PreHyper_PostHypo"))))  ->t.or.sel$Cat
factor(as.character(t.or.sel$Term), levels=selected$Term) ->t.or.sel$Term

save(t.or.sel, file='GO.categories.rda')

load('GO.categories.rda')
load('Integration.table.rda') ## Load df with values from Integration figure
load('Integration.table.simple.rda')

df$Cat2<-relabel(df$Cat, pre.UP_post.DOWN ='PreHyper_PostHypo', pre.EQUAL_post.DOWN='PreEqual_PostHypo', pre.DOWN_post.DOWN='PreHypo_PostHypo',
                 pre.DOWN_post.UP ='PreHypo_PostHyper', pre.UP_post.UP='PreHyper_PostHyper', pre.EQUAL_post.UP='PreEqual_PostHyper',
                 pre.DOWN_post.EQUAL='PreHypo_PostEqual', pre.UP_post.EQUAL='PreHyper_PostEqual')

## Simple version    				
df$Cat2<-relabel(df$Cat, pre.UP_post.EQUAL='PreHyper_PostEqual', pre.DOWN_post.EQUAL='PreHypo_PostEqual')

df$direction<-rep(1,dim(df)[1])
df$direction[which(df$OR<1)]<-2
bonf<-0.05/dim(df)[1]
df$col<-rep(1,dim(df)[1])
df$col[which(df$P<0.05)]<-2	

df[,c(9,2,5,8,10)] -> part1
colnames(part1)<-c('Cat','Term','logOR','Sig','Direction')
t.or.sel[,c(1,4,5,12,13)] ->part2
colnames(part2)<-c('Cat','Term','logOR','Sig','Direction')

wholeset<-rbind( part1, part2 )

## Complete figure for extended
wholeset<-rbind( part1, part2[-which(part2$Cat %in% c('PostUP', 'PostDOWN')),] )
with(wholeset, factor(Cat, levels=rev(c("PreHyper_PostEqual", "PreHypo_PostEqual", "PreEqual_PostHyper", "PreHyper_PostHyper", "PreHypo_PostHyper", "PreHypo_PostHypo", "PreEqual_PostHypo", "PreHyper_PostHypo"))))  -> wholeset$Cat

## Simple for main text
wholeset<-wholeset[which(wholeset$Cat %in% c("PreHyper_PostEqual", "PreHypo_PostEqual", "PostUP", "PostDOWN")),]
with(wholeset, factor(Cat, levels=rev(c("PreHyper_PostEqual", "PreHypo_PostEqual", "PostUP", "PostDOWN"))))  -> wholeset$Cat 
wholeset$Term <- with(wholeset, factor(Term, levels=c('TSS', 'Enhancers', 'EnhFetal', 'EnhAdult','Neuronal', 'Glial', 'NUM', 'GUM',as.character(selected$Term)))) 

myColors=c("#d10233", "#334A96")
ggplot(wholeset, aes(x=Term, y=Cat, size=abs(logOR), fill=factor(Direction), alpha=Sig),guide=FALSE)+
  geom_point(shape=21) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major.y = element_blank(), 
                                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"), axis.text.x = element_text(angle = 45, hjust = 1)) + scale_alpha_discrete(range = c(0.2, 1)) + scale_fill_manual(values = myColors) + scale_size_continuous(range = c(1, 20)) + scale_colour_manual(values = "black")  + theme(panel.grid.major.x = element_line(colour = "blue"))


### Compose the whole figure to share legend MET + ENH
wholeset->a

load<-load('Enhancer.Table2plotMet.July.2017.ONLY_ENH.rda')
wholeset.e<-get(load)
rm(load)

wholeset.e$Term<-relabel(wholeset.e$Term, Neuron='Neuronal', Glia='Glial', Fetal='EnhFetal', Adult='EnhAdult')

rbind(a, wholeset.e)->whole
with(whole, factor(Cat, levels=rev(c("PreHyper_PostEqual", "PreHypo_PostEqual", "PostUP", "PostDOWN", 'DFC_Fetal', 'DFC_Adult'))))  -> whole$Cat

myColors=c("#d10233", "#334A96")
ggplot(whole, aes(x=Term, y=Cat, size=abs(logOR), fill=factor(Direction), alpha=Sig),guide=FALSE)+
  geom_point(shape=21) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"), axis.text.x = element_text(angle = 45, hjust = 1)) + scale_alpha_discrete(range = c(0.2, 1)) + scale_fill_manual(values = myColors) + scale_size_continuous(range = c(1, 20)) + scale_colour_manual(values = "black")  

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
dir=args[1]
type=args[2]

# WGCNA module enrichment analysis ---------------------------------------------------------------

## function that calculates the significance of the overlap (Fisher's exact test)
testSigProportion2 <- function(set1, set2, Universe){
    
    set1 <- intersect(Universe, set1)
    set1.C <- setdiff(Universe, set1)
    set2 <- intersect(Universe, set2)
    
    set1.in <- intersect(set1, set2)
    set1.out <- setdiff(set1, set2)
    set1.C.in <- intersect(set1.C, set2)
    set1.C.out <- setdiff(set1.C, set2)
    
    a <- length(set1.in)
    b <- length(set1.out)
    c <- length(set1.C.in)
    d <- length(set1.C.out)
    count.dat <- rbind(c(a,b), c(c,d)) # matrix( cb(a, b, c, d), 2, 2)
    result <- fisher.test(count.dat, alternative="greater")
    print(count.dat)
    print(result)
    
    pval <- result$p.value
    OR <- result$estimate
    CI.m95 <- result$conf.int[1]
    CI.p95 <- result$conf.int[2]
    Common.genes.name <- paste(set1.in, collapse=",")
    
    list(p.value = pval, Odd.ratio = OR, CI.m95 = CI.m95 , CI.p95 = CI.p95, Common.genes = length(set1.in),
    set1.bg = length(set1), set2.bg = length(set2), overlap.bg = length(set1.in), bg=length(Universe), Common.genes.name=Common.genes.name)
}

testWGCNAenrichment2 <- function(allMods, Universe, target){
    
    ## Fot the test inside this function
    #allMods = allMods.sym
    #Universe <- allMods.genes.sym
    #target <- geneset.group
    
    M <- dim(allMods)[2]
    testInfo <- c("ModID", "OR", "CI.m95", "CI.p95", "Fisher.pvalue", "N.module",
    "N.target", "N.module.bg", "N.target.bg", "N.overlap.bg", "N.bg", "Bonferroni.pvalue", "FDR.pvalue", "Genes.overlap")
    testRes <- as.data.frame(matrix(nr=M, nc=length(testInfo)))
    colnames(testRes) <- testInfo
    
    for (kk in 1:M){
        modName <-  allMods[1,kk]
        set1 <- na.omit(unique(allMods[,kk]))[-c(1)]
        set2 <- target
        
        res <- testSigProportion2(set1, set2, Universe)
        testRes[kk, "ModID"] <- modName
        testRes[kk, "OR"] <- as.numeric(res$Odd.ratio)
        testRes[kk, "CI.m95"] <- as.numeric(res$CI.m95)
        testRes[kk, "CI.p95"] <- as.numeric(res$CI.p95)
        testRes[kk, "Fisher.pvalue"] <- res$p.value
        testRes[kk, "N.module"] <- length(set1)
        testRes[kk, "N.target"] <- length(set2)
        testRes[kk, "N.module.bg"] <- res$set1.bg
        testRes[kk, "N.target.bg"] <- res$set2.bg
        testRes[kk, "N.overlap.bg"] <- res$overlap.bg
        testRes[kk, "N.bg"] <- res$bg
        testRes[kk, "Genes.overlap"] <- res$Common.genes.name
    }
    testRes[, "Bonferroni.pvalue"] <- p.adjust(testRes[, "Fisher.pvalue"], "bonferroni")
    testRes[, "FDR.pvalue"] <- p.adjust(testRes[, "Fisher.pvalue"], "fdr")
    
    return(testRes)
}


## Load necessary libraries
library(gdata)
library(WriteXLS)
library(dplyr)

## Data and Result directories
data.dir <- "~/Desktop/YALE/BRAINSPAN/WGCNA_ENRICH/DATA/" #(Geneset file location - modified!)
result.dir <- "~/Desktop/YALE/BRAINSPAN/WGCNA_ENRICH/ANALYSIS/" #(output - modified!)

## BrainSpan WGCNA module --------------------------------------------------------------
df <- read.csv("~/Desktop/YALE/BRAINSPAN/WGCNA_ENRICH/brainSpan.genexp.brainWGCNA.geneModule.reassigned.sorted.csv", head=F, stringsAsFactors=F) # (WGCNA module - modified!)

## WGCNA genes per module: Ensemble ID|Gene symbol
allMods <- apply(df, 1, function(x) as.character(x))
allMods.genes <- unique( as.vector(t(allMods[-c(1),] ) ) )

## WGCNA genes per module: Ensemble ID
allMods.ens <- allMods
for (kk in 1:74){ allMods.ens[,kk] <- sapply(strsplit(allMods[,kk], "|", fixed=T), "[", 1) }
allMods.genes.ens <- na.omit(unique( as.vector(t(allMods.ens[-c(1),] ) ) ))

## WGCNA genes per module: Gene symbol
allMods.sym <- allMods
for (kk in 1:74){ allMods.sym[,kk] <- sapply(strsplit(allMods[,kk], "|", fixed=T), "[", 2) }
allMods.sym[1,] <- allMods[1,]
allMods.genes.sym <- na.omit(unique( as.vector(t(allMods.sym[-c(1),] ) ) ))


#######################################################################
##
## WGCNA enrichment analysis (Example below is for "GeschwindGenes")
##
## Genesets are under data.dir/Gabriel/GeschwindGenes/ENS for ensemble
## Genesets are under data.dir/Gabriel/GeschwindGenes/SYN for gene symbol
##
## For all the genesets - run with ensemble id
##
#######################################################################

## GWAS Gabriel (The rest of GWAS; not ASD)
dataset.path <- file.path(data.dir, dir, type) #(Geneset file location - modified!)
dataset.files <- list.files(dataset.path)
genesets <- list()

for (file in dataset.files){
   
    dataset.name <- sapply(strsplit(file, ".", fixed=T),"[[", 1)
    print(dataset.name)
    
    to.file <- file.path(dataset.path, file)
    genes.in.dataset <- read.table(to.file, as.is=T, head=F)$V1
    genes.in.dataset <- sapply(strsplit(genes.in.dataset, ".", fixed=T),"[[", 1)
    genesets[[dataset.name]] <- genes.in.dataset
}

groups <- names(genesets)
if(type=="ENS"){Universe <- allMods.genes.ens}
if(type=="GSY"){Universe <- allMods.genes.sym}
result <- list()

table<-c()
for(group in groups){
    
    geneset.group <- unique(unlist(genesets[[group]]))
    if(type=="GSY"){
	result[[group]] <- testWGCNAenrichment2(allMods = allMods.sym, Universe = Universe, target=geneset.group)
    }
    if(type=="ENS"){
        result[[group]] <- testWGCNAenrichment2(allMods = allMods.ens, Universe = Universe, target=geneset.group)
    }
	table<-cbind(table,result[[group]]$Fisher.pvalue)
}

filename <- paste("WGCNAenrichment", dir, "xls", sep=".") #(WGCNA result filename/location - modified!)
file.path <- file.path(result.dir, filename)
groups.label <- groups
WriteXLS(result, file.path, SheetNames=groups.label)
colnames(table)<-groups
write.table(file= paste0(result.dir,"/WGCNAenrichment_", dir, ".txt"), table, quote=F,sep="\t" )

## Get genes sets overlapping with modules
result.dir <- "~/Desktop/YALE/BRAINSPAN/WGCNA_ENRICH/GenesOverlap" #(Geneset overlap - modified!)
groups <- names(result)
ModID.of.interest <- "ME37"

for(group in groups){
    
    res <- result[[group]]
    
    res.sub <- subset(res, ModID == ModID.of.interest)
    ModGenes <- unlist(strsplit(res.sub$Genes.overlap, ",", fixed=T))
    
    cat(group, length(ModGenes), "\n")
    
    if (length(ModGenes) >=1 ){
        #ModGenes<-paste0(ModGenes,"$")
        ModGenes.matches <- unique (grep(paste(ModGenes,collapse="|"), allMods.genes, value=TRUE))
        
        outfilename <- paste(group, ModID.of.interest, "GenesOverlap.txt", sep=".")
        outpath <- file.path(result.dir, outfilename)
        write.table(ModGenes.matches, outpath, quote=F, row.names=F, col.names="GenesOverlap")
    }else{
        outfilename <- paste(group, ModID.of.interest, "GenesOverlap.txt", sep=".")
        outpath <- file.path(result.dir, outfilename)
        write.table("None", outpath, quote=F, row.names=F, col.names="GenesOverlap")
        
    }
    
}







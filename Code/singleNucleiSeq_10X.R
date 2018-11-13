##AIM: conduct single nuclei data processing, clustering, and differential expression analysis


###############----run cell ranger
###---build references for nucleus model
awk 'BEGIN{FS="\t"; OFS="\t"} $3 == "transcript"{ $3="exon"; print}' gencode.v21.annotation.gtf  > gencode.v21.annotation.edited.premrna.gtf
cellranger mkref --genome=hg38_nucleusModel --fasta=hg38.fa --genes=gencode.v21.annotation.edited.premrna.gtf
	   
###gene count
cellranger count --id=? --transcriptome=? --fastqs=? --sample=? --expect-cells=10000


######################
####seurat analysis

library(Seurat)
library(plyr)
library(dplyr)
library(Matrix)
library(Biobase)

para = 1.0
 
#####-----loading raw data
library("cellrangerRkit");
####kl
getGexpr <- function(samp_code){
	flag = samp_code;
	setFold = "./raw"
	inFile1 = paste(setFold, flag, "filtered_gene_bc_matrices/./matrix.mtx",sep="/");
	inFile2 = paste(setFold, flag, "filtered_gene_bc_matrices/./genes.tsv",sep="/");
	inFile3 = paste(setFold, flag, "filtered_gene_bc_matrices/./barcodes.tsv",sep="/");
	gbm <- NULL;
	gene <- NULL;
	gexpr.each <- NULL;
	gbm = load_cellranger_matrix_from_files(inFile1, inFile2, inFile3);
	gexpr.each =  exprs(gbm)
	gene = fData(gbm);
	rownames(gexpr.each) = toupper(gene$symbol);
	colnames(gexpr.each) = paste(flag,colnames(gexpr.each),sep="_");
	return(gexpr.each);
}
gexpr.kl <- NULL;
gexpr.op <- NULL;
gexpr.qr <- NULL;
gexpr.kl = getGexpr("kl");
gexpr.op = getGexpr("op");
gexpr.qr = getGexpr("qr");
gexpr.10x = cbind(gexpr.kl, gexpr.op, gexpr.qr);
sampList = colnames(gexpr.10x);

####-----save raw data
gexpr.10x.humanAdult = gexpr.10x;

###----
outLoc = "./mainCluster";

###----figure
figFile = paste0(outLoc, "/", "humanNuclei.10x.seuratReanalyze.HSB106DFC.pdf");
pdf(figFile,7,7);
par(omi=c(0.1,0.1,0.1,0.1));



##############------loading Seurat analysis pipeline
gexpr <- NULL
gexpr <- CreateSeuratObject(raw.data = gexpr.10x, min.cells = 50,  project = "humanNuclei")


# Calculate percentage of mitochondrial genes for filtering and add into object@meta.data,
mito.genes <- grep(pattern = "MT-", x = rownames(x = gexpr@data), value = TRUE)
percent.mito <- apply(gexpr@data[mito.genes, ], 2, sum) / apply(gexpr@data, 2, sum);
gexpr <- AddMetaData(object = gexpr, metadata = percent.mito, col.name = "percent.mito")
p0 <- VlnPlot(object = gexpr, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
print(p0);
GenePlot(object = gexpr, gene1 = "nUMI", gene2 = "nGene")

# We filter out cells that have unique gene countsover 2,500 and under 500, and > gexpr5% mitochondrial percentage
gexpr <- FilterCells(object = gexpr, subset.names = c("nGene", "nUMI", "percent.mito"),
					low.thresholds = c(300,300, -Inf), high.thresholds = c(7000,20000, 0.05))

# Normalize the data
gexpr <- NormalizeData(gexpr,normalization.method = "LogNormalize", scale.factor = 10000)

# Choose gene outliers on mean-variability plot
gexpr <- FindVariableGenes(object = gexpr, mean.function = ExpMean, dispersion.function = LogVMR, 
						x.low.cutoff = 0.1, x.high.cutoff = 8, y.cutoff = 1, do.contour = FALSE)
length(x = gexpr@var.genes)

##highly variable genes
# Perform negative-binomial regression on the variable genes, this sets their value in gexpr@scale.data, which is used for PCA/clustering
# We only do this on the variable genes to save time, but you can do this genome-wide
# We treat mitochondrial percentage, batch, and nUMI as confounding variables,
# We treat mitochondrial percentage, and nUMI as confounding variables,
#gexpr <- ScaleData(object = gexpr, vars.to.regress = c("percent.mito", "nUMI"), genes.use = gexpr@var.genes, model.use = "negbinom")
gexpr <- ScaleData(object = gexpr, vars.to.regress = c("percent.mito"),display.progress = FALSE, genes.use = gexpr@var.genes)


###----run PCA analysis
gexpr <- RunPCA(object = gexpr, pc.genes = gexpr@var.genes, pcs.compute = 100, do.print = F);

#PrintPCA(object = gexpr, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
#VizPCA(object = gexpr, pcs.use = 1:2)
PCAPlot(object = gexpr, dim.1 = 1, dim.2 = 2, group.by="orig.ident")

###projection
gexpr <- ProjectPCA(object = gexpr, do.print = FALSE)

###heatmap
PCHeatmap(object = gexpr, pc.use = c(1:25), cells.use = 500, do.balanced = TRUE, 
		label.columns = FALSE, use.full = FALSE)

###-----PCs significance
gexpr <- JackStraw(object = gexpr, num.replicate = 100)
p0 <- PCElbowPlot(object = gexpr, num.pc=100);
print(p0);

###----tsne analysis
gexpr <- RunTSNE(object = gexpr, reduction.use = "pca", dims.use = 1:25, 
				do.fast = TRUE,do.label=T, max_iter = 2000)
				
###----SNN-clip find clusters
gexpr <- FindClusters(object = gexpr, reduction.type = "pca", dims.use = 1:25, 
					resolution = para, print.output = FALSE)
#PrintFindClustersParams(object = gexpr)				
	
TSNEPlot(object = gexpr,pch=".", group.by = "orig.ident");
TSNEPlot(object = gexpr,pch=".", do.label=T,label.size = 8);
	
#####-----save data
save(gexpr, file=paste0(outLoc, "/", "humanNuclei.seuratReanalyze.HSB106DFC.RData"));

#####-----plot markers
##gene info
geneList = rownames(gexpr@data);
rec = unlist(strsplit(geneList, split="|",fixed=T))
dim(rec) = c(2, length(rec)/2);geneSymbol = rec[2,];

celltypeMarker = c("GAD1","GAD2", "SATB2","SLC17A7", "GFAP","AQP4","OLIG1","OLIG2","MBP"
,"CLDN5","FLT1","LY6C1","PTRF","B2M","IGFBP7"
,"APBB11P","CD74","GPR34","GPR183","P2RY12","HEXB","LRF5"
,"MKI67","TOP2A","HBA2");

for(k in 1:length(celltypeMarker)){
	xk = which(as.character(geneSymbol) == celltypeMarker[k]);
	if(length(xk) > 0){
		xk = xk[1];
		geneName = as.character(geneList[xk]); 
		FeaturePlot(object = gexpr, features.plot = geneName,pch.use=".", pt.size=0.1, no.legend = T,cols.use = c("grey", "red"))
	}
}	
dev.off()


##############------------------------
#####differential expression analysis
library(Seurat)
library(plyr)
library(dplyr)
library(Matrix)
library(Biobase)
source("spec.R")


###----
outLoc = "./HSB106DFC/mainCluster";


#####-----loading data
load(paste0(outLoc, "/", "humanNuclei.seuratReanalyze.HSB106DFC.RData"));
gexpr = gexpr;
gexprMat = gexpr@data
meta = gexpr@meta.data;
clsList = paste("cls", as.character(meta[,ncol(meta)]),sep="_");

####------calculate spec score
genexp = gexprMat;
header = as.character(clsList);
rec = getAllSpec(genexp, header)

###---save data
res = rec[[1]];
outFile = paste0(outLoc, "/", "humanNuclei.dexMaincluster.HSB106DFC_1.0.specScore.xls");
res2 = cbind(rownames(res), res);
colnames(res2) = c("Geneid", colnames(res));
write.table(res2, file=outFile, quote=F, col.names=T, row.names=F, sep="\t",append=F);











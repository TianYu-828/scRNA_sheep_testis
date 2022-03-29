###Scripts of data progressing for Single-Cell RNA Sequencing of the Mongonia Sheep Testis

########  Single-Library Analysis with CellRanger in server
########  CellRanger version=3.1.0

nohup /home/.../10xgenomics/cellranger-3.1.0/cellranger count --id=SampleName \
--localcores 40 \
--transcriptome=/lab412E/disk_A/sheep_NCBI_references/Oar_rambouillet \
--fastqs=/home/.../Singlecell_Rawdata \
--sample=SampleName &

#############################################################################################
         #########///////=========Seurat version=4.0.2========\\\\\\\#############
#############################################################################################
### Detail information see online vignettes of Seurat
### https://satijalab.org/seurat/vignettes.html

library(dplyr)
library(Seurat)
library(patchwork)

###Load Data
Seurat_Object <- Read10X(data.dir = "C:\\...\\filtered_feature_bc_matrix\\")
Seurat_Object <- CreateSeuratObject(counts = Seurat_Object, project = "Seurat_Object", min.cells = 3, min.features = 200)

###Value setting of filteration
Seurat_Object <- subset(Seurat_Object, subset =  nFeature_RNA > 200 & nCount_RNA >3500 & nCount_RNA < 100000)

###Perform integration
Object.anchors <- FindIntegrationAnchors(object.list = list(Object1,Object2), dims = 1:15)
Object.integrated <- IntegrateData(anchorset = Object.anchors, dims = 1:15)
DefaultAssay(object = Object.integrated) <- "integrated"
Object.integrated <- ScaleData(Object.integrated, verbose = FALSE)
Object.integrated <- RunPCA(Object.integrated, npcs = 30, verbose = FALSE)
Object.integrated <- RunUMAP(Object.integrated, reduction = "pca", dims = 1:15)

#### UMAP and Clustering
Object.integrated <- FindNeighbors(Object.integrated, reduction = "pca", dims = 1:20)
Object.integrated <- FindClusters(Object.integrated, resolution = 0.5)
p1 <- DimPlot(Object.integrated, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(Object.integrated, reduction = "umap", group.by = "seurat_clusters")
plot_grid(p1, p2)

###find markers
DefaultAssay(Object.integrated) <- "RNA"
Object.integrated.All.Markers <- FindAllMarkers(Object.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

Object.integrated<-BuildClusterTree(Object.integrated)
Tool(object = Object.integrated, slot = 'BuildClusterTree')
PlotClusterTree(Object.integrated)

###subset single cell type
cell_type <- subset(Object.integrated, idents = c("cell_type"))
saveRDS(cell_type, file = "C:\\...\\cell_type.rds")

########single cell type with seurat
Seurat_Germs <- readRDS(file = "C:\\...\\Germ_cells.rds")
DefaultAssay(Seurat_Germs) <- "integrated"
Seurat_Germs <- ScaleData(Seurat_Germs, verbose = FALSE)
Seurat_Germs <- RunPCA(Seurat_Germs, npcs = 30, verbose = FALSE)
Seurat_Germs <- FindNeighbors(Seurat_Germs, reduction = "pca", dims = 1:15)
Seurat_Germs <- FindClusters(Seurat_Germs, resolution = 0.5)
Seurat_Germs <- RunTSNE(object = Seurat_Germs, dims.use = 1:15, do.fast = TRUE)
Seurat_Germs <- RunUMAP(Seurat_Germs, reduction = "pca", dims = 1:15)
p3 <- DimPlot(Seurat_Germs, reduction = "umap", group.by = "orig.ident")
p4<- DimPlot(Seurat_Germs, reduction = "umap", label = TRUE,pt.size = 1.4)
plot_grid(p3)
plot_grid(p4)

###Species homeotic gene transformation
library(biomaRt)
count_raw <- as.matrix(Species@assays$RNA@counts)
usGenes  <- rownames(count_raw)
ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host="http://apr2018.archive.ensembl.org")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
Species = useMart("ensembl", dataset = "Species_gene_ensembl")
myout = getLDS(attributes = c("Species_symbol"), filters = "Species_symbol", values = usGenes , mart = Species, attributesL = c("hgnc_symbol","ensembl_gene_id"), martL = human, uniqueRows=T)

###Multi-species correlation analysis
library(ggcorrplot)
corr <-read.csv(file = "C:\\...\\all_Cor.csv")
head(corr)
cor_data <- cor(corr,method="pearson")
head(cor_data)
ggcorrplot(cor_data)

###Annotate gene information and calculate the ratio of X chromosome to autosomal genes
library(AnnotationHub )
count_raw <- as.matrix(Germpart@assays$RNA@counts)
geneid <- rownames(count_raw)
genes_chr <- mapIds(org.Oa.eg.db,keys = geneid,column = "CHR",keytype = "SYMBOL",multiVals = "first")

chr9_avg <- colSums(Germpart@assays$RNA@counts[chr9,])
chrX_avg <- colSums(Germpart@assays$RNA@counts[chrX,])
chr_avg <- colSums(Germpart@assays$RNA@counts[chr,])
Germpart$percent_chr9<- chr9_avg/chr_avg
Germpart$percent_chrX<- chrX_avg/chr_avg
VlnPlot(Germpart, features = c("percent_chr9","percent_chrX"),pt.size = 0)

###Analysis of differentially expressed genes in cell clusters
Differentially.expressed.genes <- FindMarkers(Germpart, ident.1 = "feature_1", ident.2 = "feature_2", verbose = FALSE,logfc.threshold = 0.25)

#############################################################################################
         #########///////=========Monocle version=2.18.0========\\\\\\\#############
#############################################################################################
### Detail information see online vignettes of Monocle
### http://cole-trapnell-lab.github.io/monocle-release/docs/
library(Seurat)
library(monocle)

Subset_Germs <- readRDS(file = "C:\\...\\Seurat.rds")

###Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(Subset_Germs@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = Subset_Germs@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
monoGERM <- newCellDataSet(data,
                           phenoData = pd,
                           featureData = fd,
                           lowerDetectionLimit = 0.5,
                           expressionFamily = negbinomial.size())
Germpart <- FindVariableFeatures(Subset_Germs, selection.method = "mean.var.plot", nfeatures = 2000)
VariableFeatures(Germpart)
seurat_var_genes = VariableFeatures(Germpart)
head(seurat_var_genes)
monoGERM <- estimateSizeFactors(monoGERM)
monoGERM <- estimateDispersions(monoGERM)
germ_seur_var = setOrderingFilter(monoGERM, seurat_var_genes)
plot_ordering_genes(germ_seur_var)
germ_cds <- reduceDimension(germ_seur_var, method = 'DDRTree')
germ_cds <- orderCells(germ_cds)
plot_cell_trajectory(germ_cds)
plot_cell_trajectory(germ_cds, color_by = "seurat_clusters")

###BEAM Function
BEAM_res <- BEAM(germ_cds, branch_point = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

library(colorRamps)
library(RColorBrewer)
cols<-colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(65)
plot_genes_branched_heatmap(germ_cds[row.names(subset(BEAM_res,qval < 1e-4)),],
                                         branch_point = 1,
                                         num_clusters = 4,
                                         cores = 1,
                                         hmcols = cols,
                                         use_gene_short_name = F,
                                         branch_colors = c("#bebebe","#009E73", "indianred2"),
                                         show_rownames = F,return_heatmap = T)


#############################################################################################
         #########///////=========SCENIC (version 1.2.4)========\\\\\\\#############
#############################################################################################
######The following steps are executed according to SCENIC pipeline with default parameters. 
#####Please refer to the url below: 
### https://github.com/aertslab/SCENIC
### https://rawcdn.githack.com/aertslab/SCENIC/6aed5ef0b0386a87982ba4cc7aa13db0444263a6/inst/doc/SCENIC_Running.html

library(Seurat)
library(monocle)
library(cowplot)
library(SCopeLoomR)
library(SCENIC)
library(GENIE3)
library(RcisTarget)
library(AUCell)
library(SingleCellExperiment)
library(zoo)
library(NMF)
library(doMC)
library(doRNG)
library(Rtsne)
library(R2HTML)
library(RcisTarget)
library(AUCell)
library(mixtools)
library(rbokeh)

Germ<-readRDS("C:\\...\\Germ.rds")
data <- as(as.matrix(Germ@assays$RNA@counts), 'sparseMatrix')
expr_Mat <- as.matrix(data)

cellInfo <- expr_Mat@meta.data
cellInfo$nGene <- colSums(expr_Mat>0)
cellInfo <- data.frame(cellInfo)

dir.create("int")
saveRDS(cellInfo, file="int/cellInfo.Rds")

colVars <- list(celltype=setNames(c("hotpink",
                                 "skyblue",
                                 "SeaGreen1"), 
                               c("cell_type_1", 
                                 "cell_type_2",
                                 "cell_type_3")))
								 
saveRDS(colVars, file="int/colVars.Rds")
plot.new(); legend(0,1, fill=colVars$State, legend=names(colVars$State))

org="hgnc" # or hgnc, or dmel
dbDir="./crisTarget_databases"  # RcisTarget databases location
myDatasetTitle="SCENIC of human germ cell" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]

scenicOptions <- initializeScenic(org="hgnc",dbDir="/lab412E/disk_A/crisTarget_databases",
                                  datasetTitle="SCENIC of human germ cell",
                                  nCores=30) 
								  
hgnc_dbs <- list('500bp' = 'hg19-500bp-upstream-7species.mc9nr.feather',
                 '10kb' = 'hg19-tss-centered-10kb-7species.mc9nr.feather')
db_mcVersion <- 'v9'
db_path <- './crisTarget_databases'
scenicOptions@settings$dbs <- hgnc_dbs
scenicOptions@settings$dbDir <- db_path
scenicOptions@settings$db_mcVersion <- db_mcVersion

scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"

saveRDS(scenicOptions, file="int/scenicOptions.Rds")

genesKept <- geneFiltering(expr_Mat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(expr_Mat),
                           minSamples=ncol(expr_Mat)*.01)
						   
exprMat_filtered <- expr_Mat[genesKept, ]	

runCorrelation <- function(exprMat_filtered,scenicOptions)
{ 
  corrMat <- cor(t(exprMat_filtered), method="spearman")
  saveRDS(corrMat, file=getIntName(scenicOptions, "corrMat"))
}

runCorrelation(exprMat_filtered, scenicOptions)		

exprMat_filtered <- log2(exprMat_filtered+1) 

# Run GENIE3
runGenie3(exprMat_filtered, scenicOptions)		

exprMat <- log2(expr_Mat+1)
dim(exprMat)
head(expr_Mat)[1:5,1:5]
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 20
scenicOptions@settings$seed <- 123	

scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] 
runSCENIC_1_coexNetwork2modules(scenicOptions)   
runSCENIC_2_createRegulons(scenicOptions) 
runSCENIC_3_scoreCells(scenicOptions, exprMat)

scenicOptions@settings$seed <- 123
# Run t-SNE with different settings:
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=c(5,15,50), perpl=c(5,15,50))
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=c(5,15,50), perpl=c(5,15,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
# Plot as pdf (individual files in int/): 
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=TRUE)

#and to view/compare them…
par(mfcol=c(3,3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=TRUE, varName="celltype", cex=.5)

# Using only "high-confidence" regulons (normally similar)
par(mfcol=c(3,3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_oHC_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=TRUE, varName="celltype", cex=.5)

###The chosen t-SNE can then be saved as default to use for plots (can also be “binary”, see below):
scenicOptions@settings$defaultTsne$aucType <- "AUC"
scenicOptions@settings$defaultTsne$dims <- 15
scenicOptions@settings$defaultTsne$perpl <- 50
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

####################################---regulon on/off---#################################
####Binarize the network activity (regulon on/off)

logMat <- exprMat # Better if it is logged/normalized
aucellApp <- plotTsne_AUCellApp(scenicOptions, logMat)
savedSelections <- shiny::runApp(aucellApp)

# Save the modified thresholds:
newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

runSCENIC_4_aucell_binarize(scenicOptions)

#############################################################################################
         #########///////=========CellChat (version 1.1.0)========\\\\\\\#############
#############################################################################################
###Sheep single-cell data should first be converted to human homologous genes

data.input<-GetAssayData(Seurat_Object,assay = "RNA",slot = "data")
meta<-Seurat_Object@meta.data
cellchat <- createCellChat(object = data.input,meta=meta,group.by = "celltype")
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels")
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.human
colnames(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

p = netVisual_bubble(cellchat, sources.use =c(7,8,9,10), targets.use = c(1,2,3,4,5,6), remove.isolate = FALSE,font.size =20,line.size = 0.8,font.size.title = 20)
p

###Species comparison
Species_1<-readRDS("C:\\...\\Species_1.rds")
Species_2<-readRDS("C:\\...\\Species_2.rds")
cco.list <- list(pbmc=Species_1_cellchat, til=Species_2_cellchat)
cellchat <- mergeCellChat(cco.list, add.names = names(cco.list), cell.prefix = TRUE)

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count") 
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight") 
p <- gg1+gg2
p

par(mfrow = c(1,2)) 
netVisual_diffInteraction(cellchat, weight.scale = T) 
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE) 
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE) 
p <- gg1+gg2
p

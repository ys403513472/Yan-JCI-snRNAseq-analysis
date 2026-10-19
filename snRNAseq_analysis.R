library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(scDblFinder)
library(SingleCellExperiment)

#load the dataset#
Control.counts <- Read10X((data.dir = "data/CFG4504/filtered_feature_bc_matrix/"))
Control <- CreateSeuratObject(counts = Control.counts)
Control@meta.data$Sample <- "Control"
Control@meta.data$treat <- "Control"

Cocaine.counts <- Read10X((data.dir = "data/CFG4506/filtered_feature_bc_matrix/"))
Cocaine <- CreateSeuratObject(counts = Cocaine.counts)
Cocaine@meta.data$Sample <- "Cocaine"
Cocaine@meta.data$treat <- "Cocaine"

#############################
##### doublets removal ######
#############################
### control ####
sce <- as.SingleCellExperiment(Control)
sce <- scDblFinder(sce)
Control$scDblFinder_class <- colData(sce)$scDblFinder.class
Control$scDblFinder_score <- colData(sce)$scDblFinder.score
Control_singlets <- subset(Control, subset = scDblFinder_class == "singlet")
table(Control_singlets@meta.data$orig.ident) 
table(Control@meta.data$orig.ident)

### Cocaine ####
sce <- as.SingleCellExperiment(Cocaine)
sce <- scDblFinder(sce)
Cocaine$scDblFinder_class <- colData(sce)$scDblFinder.class
Cocaine$scDblFinder_score <- colData(sce)$scDblFinder.score
Cocaine_singlets <- subset(Cocaine, subset = scDblFinder_class == "singlet")
table(Cocaine_singlets@meta.data$orig.ident) 
table(Cocaine@meta.data$orig.ident) 

############ Merge ###############
Merge <- merge(Control_singlets, y = Cocaine_singlets, add.cell.ids = c("ctrl", "coca"))
###### QC #######
Idents(Merge) <- "Sample"
Merge[["percent.mt"]] <- PercentageFeatureSet(Merge, pattern = "^mt-")
# Pre-calculate quantile thresholds
nFeat_high <- quantile(Merge$nFeature_RNA, 0.98)
nFeat_low  <- quantile(Merge$nFeature_RNA, 0.02)
nCount_high <- quantile(Merge$nCount_RNA, 0.98)
nCount_low  <- quantile(Merge$nCount_RNA, 0.02)

Merge.filter <- subset(
  Merge, 
  subset = nFeature_RNA >= 500 & 
    nFeature_RNA < nFeat_high & 
    nFeature_RNA > nFeat_low & 
    nCount_RNA < nCount_high & 
    nCount_RNA > nCount_low & 
    percent.mt < 15
)

#### normalization and run standard analysis workflow #####
Merge.filter <- NormalizeData(Merge.filter)
Merge.filter <- FindVariableFeatures(Merge.filter, selection.method = "vst", nfeatures = 2000)
Merge.filter <- ScaleData(Merge.filter)
Merge.filter <- RunPCA(Merge.filter)
Merge.filter <- FindNeighbors(Merge.filter, dims = 1:30, reduction = "pca")
Merge.filter <- FindClusters(Merge.filter, resolution = 0.2)
Merge.filter <- RunUMAP(Merge.filter, dims = 1:30, reduction = "pca", reduction.name = "umap")

##### cell type annotation ######
merge.label <- RenameIdents(Merge.filter, `15` ="EC",`2` ="Astro",
                              `9` ="Microglia",`12` ="OPC",`5` ="OL",
                              `14` ="exN",`16` ="MGE_Npy",`13` ="MGE_Maf",
                              `0` ="D2MSN_I",`10` ="D2MSN_II",`1` ="D1MSN_I",
                              `4` ="D1MSN_II",`3` ="eSPN",`7` ="D3MSN",`11` ="NPC",`8` ="inN",`6` ="immatureN")
merge.label$cell.type <- Idents(merge.label)
markers.to.plot <- c("Gad2","Drd1","Reln","Htr4","Tshz1","Drd2","Oprd1","Htr2a","Drd3","Dlx1",
                     "Trpc4","Maf","Npy","Slc17a7","Prr5l","Cspg4","Cx3cr1","Aldh1l1","Kdr")
levels <- c("D1MSN_I","D1MSN_II","eSPN","D2MSN_I","D2MSN_II","D3MSN","NPC","inN","MGE_Maf","MGE_Npy","immatureN","exN","OL","OPC","Microglia","Astro","EC")
merge.label@active.ident <- factor(merge.label@active.ident, levels = levels)

##### UMAP generation ######
Idents(merge.label) <- "cell.type"
DimPlot(merge.label, reduction = "umap", label = T, split.by = "Sample")

##### DEG analysis #######
# Ensure layers are joined if using Seurat v5
merge.label <- JoinLayers(merge.label)
Idents(merge.label) <- "treat"
merge.label$cluster.stim <- paste(Idents(merge.label), merge.label$cell.type, sep = "_")
Idents(merge.label) <- "cluster.stim"

cell <- read.csv("data/cell.type.csv", header = 1)
cell <- as.data.frame(cell)
cell <- split(cell, seq(nrow(cell)))

degs <- FindMarkers(merge.label, ident.1 = paste0("Cocaine_", cell[[1]]),
                    ident.2 = paste0("Control_", cell[[1]]), test.use = "MAST", min.pct = 0.2, logfc.threshold = 0.25)
degs$cluster <- as.character(cell[[1]])
degs$gene <- row.names(degs)
df <- degs

for (i in 2:18)
{
  degs <- FindMarkers(merge.label, ident.1 = paste0("Cocaine_", cell[[i]]),
                      ident.2 = paste0("Control_", cell[[i]]), test.use = "MAST", min.pct = 0.2, logfc.threshold = 0.25)
  degs$cluster <- as.character(cell[[i]])
  degs$gene <- row.names(degs)
  temp <- degs
  df <- rbind(df, temp)
}

# filter deg by p.adj
deg.sig <- df[df$p_val_adj < 0.05 & abs(df$avg_log2FC) > 0.38, ]

#### dotplot for Drd1 and Drd2 expression in astrocyte ###
Idents(merge.label) <- "cell.type"
astro <- subset(merge.label, idents = c("Astro"))

genes_to_plot <- c("Drd1","Drd2")
Idents(astro) <- "treat"
p <- DotPlot(astro, features = genes_to_plot,cols = c("grey70", "red")) + theme(
  axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
  axis.text.y = element_text(size = 16),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16))
print(p)



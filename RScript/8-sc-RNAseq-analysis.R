###################################################
## 
## Single cell validation
##
## - Analysis of data from Sade-Feldman M et al. 
## - Manual annotion of cells
## - Increasing resolution of B cell lineage
## - Mapping IMK gene signature
##
## Manuscript --> Figure 3, Figure S3, Table S7, Table S8
##
###################################################

set.seed(123)


# scRNA-seq analysis ------------------------------------------------------


suppressPackageStartupMessages(library(Seurat))
library(ggplot2)
library(dplyr)
library(ggpubr)
library(RColorBrewer)

org.options <- options

# add some convenience functions
setPlotSize <- function(height, width) {
  options(repr.plot.height=height, repr.plot.width = width)
}

resetPlotSize <- function() {
  options(org.options)
}

my_theme_bw <- theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()
        , axis.title.x=element_text(size=16, face = "plain")
        , axis.title.y=element_text(size=16, face = "plain")
        , axis.text.y=element_text(size=16)
        , axis.text.x=element_text(size=16)
        , plot.title=element_text(size=16, face = "bold", hjust = 0.5)
        , legend.text=element_text(size=16)
        , legend.title=element_text(size=16, face = "bold"))



GSE120575 <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE120575&format=file&file=GSE120575%5FSade%5FFeldman%5Fmelanoma%5Fsingle%5Fcells%5FTPM%5FGEO%2Etxt%2Egz"
GSE120575_id <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE120575&format=file&file=GSE120575%5Fpatient%5FID%5Fsingle%5Fcells%2Etxt%2Egz"

download.file(GSE120575, destfile = "../data/GSE120575_TPM.txt")
download.file(GSE120575_id, destfile = "../data/GSE120575_id.txt")

ID <- read.csv("../data/GSE120575_TPM.txt", header = F, sep ="\t", stringsAsFactors = FALSE, nrows = 2)
ID$V1 <- NULL
time <- as.character(ID[2, ])
ID <- as.character(ID[1, ])
exprs.data <- read.csv("../data/GSE120575_TPM.txt", header = F, sep ="\t", stringsAsFactors = FALSE, skip = 2)
rownames(exprs.data) <- exprs.data$V1
exprs.data$V1 <- NULL
colnames(exprs.data) <- ID

# load the sample metadata
sample.data <- read.csv("../data/GSE120575_id.txt", sep = "\t", stringsAsFactors = F, skip = 19, header = T)
sample.data <- sample.data[1:16291, ]
rownames(sample.data) <- sample.data[, "title"]
set.seed(123)
# keep pre-treatment patients
idt <- grep("Pre", sample.data$characteristics..patinet.ID..Pre.baseline..Post..on.treatment.)
pre.patient <- rownames(sample.data)[idt]

exprs.data <- exprs.data[, pre.patient]

# remove unused columns

# sanity check to make sure that all samples have corresponding annotations
if (!all(colnames(exprs.data) %in% rownames(sample.data))) {
  stop("Missing sample annotations")
}

# order the sample metadata like the columns in the epxression matrix
sample.data <- sample.data[colnames(exprs.data), ]


# Create seurat object 

sade <- CreateSeuratObject(counts = exprs.data, 
                           min.cells = 5, 
                           min.features = 200, 
                           project = "Sade-Feldman",
                           names.delim = " ",
                           meta.data = sample.data)

# compute MT
sade[["percent.mt"]] <- PercentageFeatureSet(sade, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(sade, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(sade, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sade, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# normalize
sade <- NormalizeData(object = sade, normalization.method = "LogNormalize", scale.factor = 10000)

# find the variable genes
set.seed(123)
sade <- FindVariableFeatures(sade, selection.method = "vst", nfeatures = 2000)


# correct for confounders
all.genes <- rownames(sade)
sade <- ScaleData(sade, features = all.genes)

# perform the PCA
set.seed(123)
sade <- RunPCA(sade, features = VariableFeatures(object = sade))

# Create cluster 

sade <- JackStraw(sade, num.replicate = 100)
sade <- ScoreJackStraw(sade, dims = 1:20)
JackStrawPlot(sade, dims = 1:20)
ElbowPlot(sade)


sade <- FindNeighbors(object = sade, dims=1:20, reduction = "pca", force.recalc = F)
sade <- FindClusters(object = sade, resolution = 0.8, print.output = 0, random.seed = 123)
sade <- RunUMAP(sade, dims = 1:23)
DimPlot(sade, reduction = "umap")

# Manually check annotation and compare with PanglaoDB 

sade$Lineage <- NA
sade$Lineage[sade$RNA_snn_res.0.8 == as.character(0)] = "T cell"
sade$Lineage[sade$RNA_snn_res.0.8 == as.character(1)] = "T cell"
sade$Lineage[sade$RNA_snn_res.0.8 == as.character(2)] = "T cell"
sade$Lineage[sade$RNA_snn_res.0.8 == as.character(3)] = "NK cells"
sade$Lineage[sade$RNA_snn_res.0.8 == as.character(4)] = "B cell"
sade$Lineage[sade$RNA_snn_res.0.8 == as.character(5)] = "T cell"
sade$Lineage[sade$RNA_snn_res.0.8 == as.character(6)] = "Macrophages"
sade$Lineage[sade$RNA_snn_res.0.8 == as.character(7)] = "T cell"
sade$Lineage[sade$RNA_snn_res.0.8 == as.character(8)] =  "T cell"
sade$Lineage[sade$RNA_snn_res.0.8 == as.character(9)] =  "T cell"
sade$Lineage[sade$RNA_snn_res.0.8 == as.character(10)] = "T cell"
sade$Lineage[sade$RNA_snn_res.0.8 == as.character(11)] = "Plasma cells"
sade$Lineage[sade$RNA_snn_res.0.8 == as.character(12)] = "B cells"
sade$Lineage[sade$RNA_snn_res.0.8 == as.character(13)] = "DC"

sade$CellType <- NA
sade$CellType[sade$RNA_snn_res.0.8 == as.character(0)] = "CD4 cytotoxic"
sade$CellType[sade$RNA_snn_res.0.8 == as.character(1)] =  "CD8 T cytotoxic 1"
sade$CellType[sade$RNA_snn_res.0.8 == as.character(2)] =  "CD8 T exhausted"
sade$CellType[sade$RNA_snn_res.0.8 == as.character(3)] =  "NK cells"
sade$CellType[sade$RNA_snn_res.0.8 == as.character(4)] =  "B cells"
sade$CellType[sade$RNA_snn_res.0.8 == as.character(5)] =  "CD4 self renewing"
sade$CellType[sade$RNA_snn_res.0.8 == as.character(6)] =  "Macrophages"
sade$CellType[sade$RNA_snn_res.0.8 == as.character(7)] =  "CD8 T memory"
sade$CellType[sade$RNA_snn_res.0.8 == as.character(8)] =  "T gamma delta"
sade$CellType[sade$RNA_snn_res.0.8 == as.character(9)] =  "CD8 T cytotoxic 2"
sade$CellType[sade$RNA_snn_res.0.8 == as.character(10)] = "T reg"
sade$CellType[sade$RNA_snn_res.0.8 == as.character(11)] = "Plasma cells"
sade$CellType[sade$RNA_snn_res.0.8 == as.character(12)] = "B cells"
sade$CellType[sade$RNA_snn_res.0.8 == as.character(13)] = "DC plasmacytoid"



# B cell lineage ----------------------------------------------------------


sade@active.ident <- factor(sade$CellType)
b.cells <- subset(sade, CellType %in% c("B cells", "Plasma cells"))

# perform the PCA
set.seed(123)
b.cells <- RunPCA(b.cells, features = VariableFeatures(object = b.cells))

b.cells <- FindNeighbors(object = b.cells, dims=1:14, reduction = "pca", force.recalc = T)
b.cells <- FindClusters(object = b.cells, resolution = 0.5, print.output = 0, random.seed = 123)
b.cells <- RunUMAP(b.cells, dims = 1:14)

dummy <- prueba@meta.data[b.cells$title, ]
b.cells@meta.data[dummy$title, "dummy"] <- dummy$Cell_Type

b.cells$cluster_type <- as.numeric(b.cells$RNA_snn_res.0.5)-1
b.cells$cluster_type <- ifelse(b.cells$RNA_snn_res.0.5 %in% c(2, 5), "Plasma cells", b.cells$cluster_type)
b.cells$cluster_type <- ifelse(b.cells$RNA_snn_res.0.5 %in% c(1), "Plasmablast", b.cells$cluster_type)
b.cells$cluster_type <- ifelse(b.cells$RNA_snn_res.0.5 %in% c(0, 7), "Naive B cells", b.cells$cluster_type)
b.cells$cluster_type <- ifelse(b.cells$RNA_snn_res.0.5 %in% c(4, 6), "Naive B cells IGK-high", b.cells$cluster_type)
b.cells$cluster_type <- ifelse(b.cells$RNA_snn_res.0.5 %in% c(3), "Naive B cells IGL-high", b.cells$cluster_type)

b.cells@active.ident <- factor(b.cells$cluster_type)

sade$Cell_Type <- as.character(sade$CellType)
sade$Cell_Type[WhichCells(b.cells, idents = c("Plasma cells"))] <- "Plasma cells"
sade$Cell_Type[WhichCells(b.cells, idents = c("Plasmablast"))] <- "Plasmablast"
sade$Cell_Type[WhichCells(b.cells, idents = c("Naive B cells"))] <- "Naive B cells"
sade$Cell_Type[WhichCells(b.cells, idents = c("Naive B cells IGK-high"))] <- "Naive B cells IGK-high"
sade$Cell_Type[WhichCells(b.cells, idents = c("Naive B cells IGL-high"))] <- "Naive B cells IGL-high"

#saveRDS(sade, "../data/sade.RDS")
sade <- readRDS("../data/sade.RDS")

sade$response <- sade$characteristics..response
sade$patient <- gsub("Pre_", "", sade$characteristics..patinet.ID..Pre.baseline..Post..on.treatment.)
table(sade$patient, sade$response)
info <- sade@meta.data

sade$Response <- ifelse(sade$response=="Responder", "Good", "Bad")
my.col <- c('Good'="#EFC000FF", 'Bad'="#0073C2FF")
Figure.3A2 <- DimPlot(sade, group.by = "Response", cols = my.col)+my_theme_bw

response <- info[info$response=="Responder", c("patient", "Cell_Type")]
response <- reshape2::melt(table(response$patient, response$Cell_Type))
response$response <- "Responder"
non.response <- info[info$response=="Non-responder", c("patient", "Cell_Type")]
non.response <- reshape2::melt(table(non.response$patient, non.response$Cell_Type))
non.response$response <- "Non-responder"
response <- rbind(response, non.response)


my.list <- list()
p.res <- list()
my_comparisons <- list( c("Bad", "Good") )
for(i in unique(sade$Cell_Type)){
  p <- ggviolin(response[response[,2]==i, ], x = "response", y = "value", facet.by = "Var2",
                color = "response", palette = "jco", add = "jitter", group.by = "Var2")+
    stat_compare_means(comparisons = my_comparisons, label = "p.signif", label.y = max(response[response[,2]==i, 3])*1.2)+ #
    stat_compare_means(label.y = max(response[response[,2]==i, 3])*1.2)
  p.res[[i]] <- wilcox.test(value~response, response[response[,2]==i, ])$p.value
  my.list[[i]] <- p
}

TableS7 <- data.frame(Stromal.cell.type.scRNAseq = names(p.res), Wilcox.test = unlist(p.res))

Figure.3B <- ggpubr::ggarrange(plotlist=my.list[c("Naive B cells", "Naive B cells IGL-high", "Naive B cells IGK-high")], common.legend=T, nrow = 1)


length(unique(sade$CellType))
col <- c(brewer.pal(n = 9, name = "Set1"), brewer.pal(n = 4, name = "Set2"))
names(col) <- unique(sade$CellType)
col["Plasma cells"] <- "antiquewhite2"
col["Macrophages"] <- "deepskyblue"
Figure.3A1 <- DimPlot(sade, reduction = "umap", group.by = "CellType", cols = col)+my_theme_bw



bcell <- subset(sade, Cell_Type%in%c("Naive B cells", "Naive B cells IGK-high", "Naive B cells IGL-high", "Plasma cells", "Plasmablast"))
col <- c(brewer.pal(n = 5, name = "Paired"))
names(col) <- unique(bcell$Cell_Type)
col[1] <- "goldenrod3"
col[2] <- "lightsteelblue1"
col[3] <- "antiquewhite2"
col[4] <- "salmon1"
col[5] <- "darkolivegreen1"
Figure.3A1bcell <- DimPlot(bcell, reduction = "umap", group.by = "Cell_Type", cols = col, size = 10)+my_theme_bw+xlab("")+ylab("")+ggtitle("")+theme(
  axis.ticks.x=element_blank()
  , axis.text.x=element_blank()
  , axis.text.y=element_blank()
  , axis.ticks.y=element_blank()
)



markers <- FindAllMarkers(sade)
markers.top <- unique(as.data.frame(markers %>% filter(p_val_adj<0.05, avg_log2FC>0) %>% dplyr::group_by(cluster) %>% top_n(10, avg_log2FC) %>% select(gene))[, "gene"])

my_theme <-   theme(
  axis.title.x=element_text(size=20)
  , axis.title.y=element_text(size=20)
  , axis.text.x.top =element_text(size=16)
  , axis.text.y.left =element_text(size=16)
  , legend.text=element_text(size=16)
  , plot.title=element_text(size=20, face = "bold")
  , plot.subtitle=element_text(size=0, face = "plain")
  , legend.title=element_text(size=16, face = "bold")
)

Figure.S3 <- DoHeatmap(sade, features = markers.top)+my_theme

# Mapping IMK -------------------------------------------------------------

sade@active.ident <- factor(sade$response)
mks.pd1 <- FindMarkers(sade[, sade$characteristics..therapy=="anti-PD1"], ident.1 = "Responder", ident.2 = "Non-responder", test.use	= "bimod")

genes_IMK <- read.table("../data/Result_RNASeq_DESeq2_melanoma_cutaneous.txt", sep = "\t", stringsAsFactors = F, header = T)
genes_IMK$SYMBOL[genes_IMK$SYMBOL==""] <- genes_IMK$ENSEMBL[genes_IMK$SYMBOL==""] 
rownames(genes_IMK) <- genes_IMK$SYMBOL

genes_IMK$SYMBOL[!genes_IMK$SYMBOL%in%rownames(sade)]

intersect(rownames(mks.pd1[mks.pd1$p_val_adj<0.05,]), genes_IMK$SYMBOL)
TableS8 <- mks.pd1[intersect(rownames(mks.pd1[mks.pd1$p_val_adj<0.05,]), genes_IMK$SYMBOL), ]
TableS8$Gene <- rownames(TableS8)

sade.pd1 <- subset(sade, characteristics..therapy%in%c("anti-PD1"))

p1 <- VlnPlot(sade.pd1, features = c("IGHGP"), group.by = "Cell_Type", pt.size = 1, ncol = 1)+theme(
  axis.title.x=element_text(size=20)
  , axis.title.y=element_text(size=20)
  , axis.text.y=element_text(size=15)
  , axis.text.x=element_text(size=15)
  , plot.title =element_text(face=c("bold.italic"))
  , legend.position = "none"
)+xlab("")

p2 <- VlnPlot(sade.pd1, features = c("IGHG2"), group.by = "Cell_Type", pt.size = 1, ncol = 1)+theme(
  axis.title.x=element_text(size=20)
  , axis.title.y=element_text(size=20)
  , axis.text.y=element_text(size=15)
  , axis.text.x=element_text(size=15)
  , plot.title =element_text(face=c("bold.italic"))
  , legend.position = "none"
)+xlab("")+ylab("")

p3 <- VlnPlot(sade.pd1, features = c("POU2AF1"), group.by = "Cell_Type", pt.size = 1, ncol = 1)+theme(
  axis.title.x=element_text(size=20)
  , axis.title.y=element_text(size=20)
  , axis.text.y=element_text(size=15)
  , axis.text.x=element_text(size=15)
  , plot.title =element_text(face=c("bold.italic"))
  , legend.position = "none"
)+xlab("")+ylab("")


Figure.3E <- ggpubr::ggarrange(plotlist=list(p1, p2, p3), nrow = 1)


p4 <- VlnPlot(sade.pd1, features = c("CXCR5"), group.by = "Cell_Type", pt.size = 1, ncol = 1)+theme(
  axis.title.x=element_text(size=20)
  , axis.title.y=element_text(size=20)
  , axis.text.y=element_text(size=15)
  , axis.text.x=element_text(size=15)
  , plot.title=element_text(size=20, face = "bold", hjust = 0.5)
  , plot.subtitle=element_text(size=18, face = "italic", hjust = 0.5)
  , legend.position = "none"
)+xlab("")+ylab("")+ggtitle("Single-cell", subtitle = "CXCR5")


cut <- read.delim("../data/Result_RNASeq_DESeq2_melanoma_cutaneous_all.txt")
cxcr5 <- reshape::melt(cut[cut$GeneSymbol=="CXCR5", -1])
cxcr5[, "Bulk"] <- cxcr5$value
cxcr5$Response <- gsub("[.].*", "", cxcr5$variable)

p <- ggviolin(cxcr5, x = "Response", y = "Bulk",
              color = "Response", palette = "jco", add = "jitter", caption = "DE adj-pvalue = 0.0969", font.caption=15, title = "Bulk", subtitle = "CXCR5")+theme(
                axis.title.x=element_text(size=20)
                , axis.title.y=element_text(size=20)
                , axis.text.y=element_text(size=16)
                , axis.text.x=element_text(size=16)
                , legend.text=element_text(size=16)
                , plot.title=element_text(size=20, face = "bold", hjust = 0.5)
                , plot.subtitle=element_text(size=18, face = "italic", hjust = 0.5)
                , legend.title=element_text(size=16, face = "bold")
              )+ylab("Expression Level")+xlab("")


bcell$Response <-bcell$response
cxcr5.sc <- data.frame("Single cell"=bcell@assays$RNA@data["CXCR5", ], Response=bcell$Response)


p2 <- ggviolin(cxcr5.sc, x = "Response", y = "Single.cell", title = "Single-cell",
               color = "Response", palette = "jco", add = "jitter", caption = paste0("DE adj-pvalue = ", signif(3.039091e-39, 4)), font.caption=15, subtitle = "CXCR5")+theme(
                 axis.title.x=element_text(size=20)
                 , axis.title.y=element_text(size=20)
                 , axis.text.y=element_text(size=16)
                 , axis.text.x=element_text(size=16)
                 , legend.text=element_text(size=16)
                 , plot.title=element_text(size=20, face = "bold", hjust = 0.5)
                 , plot.subtitle=element_text(size=18, face = "italic", hjust = 0.5)
                 , legend.title=element_text(size=16, face = "bold")
               )+ylab("")+xlab("")

Figure.3F <-ggpubr::ggarrange(plotlist=list(p, p2, p4), common.legend=T, legend = "none", nrow = 1)


# Results -----------------------------------------------------------------

write.table(TableS7, "../results/Tables/Table-S7.1.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(TableS8, "../results/Tables/Table-S8.txt", col.names = T, row.names = F, sep = "\t", quote = F)

ggsave("../results/Figures/Figure-3A1.pdf", Figure.3A1, width = 8, height = 7, device = "pdf")
ggsave("../results/Figures/Figure-3A1bcell.pdf", Figure.3A1bcell, width = 8, height = 7, device = "pdf")
ggsave("../results/Figures/Figure-3A2.pdf", Figure.3A2, width = 8, height = 7, device = "pdf")
ggsave("../results/Figures/Figure-3B.pdf", Figure.3B, width = 8, height = 5, device = "pdf")
ggsave("../results/Figures/Figure-3E.pdf", Figure.3E, width = 15, height = 7, device = "pdf")
ggsave("../results/Figures/Figure-3F.pdf", Figure.3F, width = 15, height = 7, device = "pdf")
pdf("../results/Figures/Figure-S3.pdf", width = 20, height = 20)
Figure.S3
dev.off()

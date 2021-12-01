###################################################
## 
## DE Melanoma cutaneous cohort
##
## - Differential expression analysis for melanoma cutaneous cohort
## 
## Manuscript --> Table S3, Figure 1C, Figure 1D
##
###################################################

## load libraries
library("DESeq2")
library("RColorBrewer")
library("pheatmap")
library(ggplot2)
library(dplyr)
library('biomaRt')
library(ggrepel)
library(ComplexHeatmap)
library(dplyr)
library(circlize)

cts <- read.delim("../data/raw_counts.txt")
rownames(cts) <- cts$ENSEMBL
cts$ENSEMBL <- NULL
colnames(cts) <- gsub("[.]", "-", colnames(cts))

sampleTable <- read.delim("../data/clinical.txt")
rownames(sampleTable) <- paste0(sampleTable$Response, "-", rownames(sampleTable))
sampleTable$condition <- factor(sampleTable$Response)
sampleTable$batch <- factor(sampleTable$Batch)
sampleTable$Melanoma.type <- sampleTable$Diagnosis

sampleTable <- sampleTable[sampleTable$Diagnosis%in%c("cutaneous"), ]
cts <- cts[, rownames(sampleTable)]

dds_2 <- DESeqDataSetFromMatrix(countData = cts,colData = sampleTable,design= ~ batch +  condition ) ##  correcting for batch effect

# Eliminacion heteroescedacidad (para poder aplicar regresiones lineales)
dds <- estimateSizeFactors(dds_2)
idx <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= 3 ##filtering the genes that have low counts 
dds <- dds[idx,]
dds <- DESeq(dds)

res <- results(dds, contrast=c("condition","Good","Bad")) ##compare gene expression between good and bad responders

padj <- 0.05
lFC <- 1.5
baseMean <- 10

cat("Pval -> ", padj, " -- lFC -> ", lFC, " -- Base mean -> ", baseMean)
table(res$padj<padj & abs(res$log2FoldChange)>lFC  & res$baseMean>baseMean)  ##filter only significant genes

## Order by adjusted p-value
res <- res[order(res$padj), ]

## Merge with normalized count data to keep counts per sample for some downstream analyses
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
resdata <- resdata %>% filter(baseMean>10) %>% filter(abs(log2FoldChange)>1.5) %>% filter(padj<0.05)


## Add gene names 
ensembl  <- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL",  dataset="hsapiens_gene_ensembl")
genes <- as.character((resdata$Row.names))

genemap <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"),
                 filters = "ensembl_gene_id",
                 values = genes,
                 mart = ensembl )

names(genemap)<-c("ENSEMBL","SYMBOL", "ENTREZ")
resdata$Row.names[!resdata$Row.names%in%genemap$ENSEMBL]
genemap <- rbind(genemap, data.frame(ENSEMBL="ENSG00000280237", SYMBOL="MIR4697HG", ENTREZ=NA))

resdata$ENSEMBL<-resdata$Row.names
resdata<-as.data.frame(resdata)

resdata_1<-merge(resdata,genemap,by="ENSEMBL") 
resf <- resdata_1[order(resdata_1$padj), ]

resf$Row.names <- NULL

#resf$SYMBOL[resf$SYMBOL==""] <- resf$ENSEMBL[resf$SYMBOL==""]
resf$SYMBOL[resf$SYMBOL==""] <- c("AL445647.1", "AC051619.7", 
                                  "AL136298.1", "AC233755.1", 
                                  "AC007938.3", "AC087762.2", 
                                  "AC084082.1", "U62631.1"
)

rownames(resf) <- resf$ENSEMBL

norm.counts <- as.data.frame(t(as.data.frame(counts(dds, normalized=TRUE))[resf$ENSEMBL, ]))
colnames(norm.counts) <- resf[colnames(norm.counts), "SYMBOL"]
norm.counts$response <- NULL
norm.counts[rownames(sampleTable), "response"] <- sampleTable$condition


res.dummy <- as.data.frame(res)

sum(is.na(res.dummy$log2FoldChange))
sum(is.na(res.dummy$padj))
res.dummy <- res.dummy[!is.na(res.dummy$padj), ]
res.dummy <- res.dummy[res.dummy$baseMean>baseMean, ]

res.dummy$SYMBOL <- rownames(res.dummy)
res.dummy[resf$ENSEMBL, "SYMBOL"] <- resf$SYMBOL


# Figure 1C ---------------------------------------------------------------

tab = data.frame(logFC = res.dummy$log2FoldChange, negLogPval = -log10(res.dummy$padj), gene=res.dummy$SYMBOL)

lFC <- 1.5
padj <- 0.05

lfc = lFC
fdr = -log10(padj)

tab$color <- NULL
tab$color <- ifelse(tab$logFC > lfc, "red", "gray48")
tab$color <- ifelse(tab$logFC < -lfc, "steelblue", tab$color)
tab$color <- ifelse(10^(-tab$negLogPval) > padj, "gray48", tab$color)

tab$genelabels <- as.character(tab$gene)
tab$genelabels <- ifelse(tab$color=="gray48", "", tab$genelabels)

tab$genelabels <- ""

cols <- c("Upregulated" = "red", "Downregulated" = "steelblue", "Nonsignificant" = "gray48")
vol <- ggplot(tab, aes(x = logFC, y = negLogPval, labels=gene))
#pdf(paste(type, "plots", "volcano_plot.pdf", sep = "/"))
Figure.1C <- vol +   
  scale_colour_manual(values = cols) +
  geom_point(size = 1, alpha = 1, na.rm = T, fill = tab$color, color = tab$color) +
  theme_bw(base_size = 14) + 
  theme(legend.position = "right"
        , axis.title.x=element_text(size=20)
        , axis.title.y=element_text(size=20)
        , axis.text.x=element_text(size=16)
        , axis.text.y=element_text(size=16)
  ) + 
  xlab(expression(log[2]("FC"))) + 
  ylab(expression(-log[10]("padj"))) +
  geom_hline(yintercept = fdr, colour="#990000", linetype="dashed") + geom_vline(xintercept = lfc, colour="#990000", linetype="dashed") + geom_vline(xintercept = -lfc, colour="#990000", linetype="dashed") +
  scale_y_continuous(trans = "log1p")+
  scale_x_continuous(breaks = sort(c(seq(-15, 15, length.out=5), -1.5, 1.5)))+
  geom_text_repel(aes(x = logFC, y = negLogPval, label = genelabels), label.size = 0.01, segment.size = 0.01, size = 4, box.padding = unit(0.05, "lines"))


# Figure 1D ---------------------------------------------------------------


# get the top genes
sigGenes <- resf[which(resf$padj<padj & abs(resf$log2FoldChange)>lFC  & resf$baseMean>baseMean), "ENSEMBL"]
sigGenes_SYMBOL <- resf[which(resf$padj<padj & abs(resf$log2FoldChange)>lFC  & resf$baseMean>baseMean), "SYMBOL"]

# filter the data for the top 200 by padj in the LRT test
plotDat <- vst(dds_2)[sigGenes,] %>% 
  assay()
z.mat <- t(scale(t(plotDat), center=TRUE, scale=TRUE))

# colour palette
myPalette <- c("royalblue3",  "ivory", "red3")
myRamp = colorRamp2(c(-2, 0, 2), myPalette)

group <- resf[resf$ENSEMBL %in% sigGenes, "log2FoldChange"]
group <- ifelse(group>0, "Upregulated", "Downregulated")
names(group) <- sigGenes

# col = list(`Response` = c("Bad" = "gray88", "Good" = "khaki3"))

ha1 = HeatmapAnnotation(`Response` = colData(dds_2)[,c("condition")], 
                        col = list(`Response` = c("Bad" = "#0073C2FF", "Good" = "#EFC000FF"))
                        , annotation_name_side = "left"
                        , annotation_name_gp = gpar(fontsize = 24)
                        , annotation_legend_param = list(labels_gp = gpar(fontsize = 24), title_gp=gpar(fontsize = 24, fontface = "bold"))
)

rownames(z.mat) <- sigGenes_SYMBOL

Figure.1D <- Heatmap(z.mat, name = "Color Key",
        col = myRamp,            
        show_row_name = T,
        cluster_columns = T,
        column_names_gp = gpar(fontsize = 18),
        split=group, row_title_gp = gpar(fontsize = 24),
        top_annotation = ha1, row_names_gp = gpar(fontsize = 10, fontface = "italic"), heatmap_legend_param = list(labels_gp = gpar(fontsize = 24), title_gp=gpar(fontsize = 24, fontface = "bold")))

Table.S3 <- res.dummy[, c("SYMBOL", "log2FoldChange", "padj")]

# Saving file for future analysis
#write.table(resf, "../data/Result_RNASeq_DESeq2_melanoma_cutaneous.txt", col.names = T, row.names = F, sep = "\t", quote = F)


# Results -----------------------------------------------------------------

write.table(Table.S3, "../results/Tables/Table-S3.2.txt", col.names = T, row.names = F, sep = "\t", quote = F)
ggsave("../results/Figures/Figure-1C.pdf", Figure.1C, width = 7, height = 6, device = "pdf")
pdf("../results/Figures/Figure-1D.pdf", width = 18, height = 18)
Figure.1D
dev.off()


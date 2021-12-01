###################################################
## 
## TCR, BCR and HLA profiling
##
## - TCR and BCR analysis from MiXCR
## - HLA analysis from seq2HLA
## - Merging HLA, TCR and BCR
##
## Manuscript --> Table 1, Table S10, Figure 5, Figure S5, Figure S6
##
###################################################



library(gdata)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(gridExtra)
library(grid)
library(lattice)
library(immunarch)
library(circlize)
library(stringr)
library(maditr)
library(RColorBrewer)


# TCR BCR analysis --------------------------------------------------------


df_mixcr <- read.delim("../data/mixcr.txt")
df_mixcr <- df_mixcr[df_mixcr$Melanoma.type=="cutaneous", ]

table(df_mixcr$Sample, df_mixcr$Type)

df_mixcr <- rbind(df_mixcr, data.frame(
  Sample=c("IMK21", "IMK21")
  , Response=c("Bad","Bad")
  , Melanoma.type=c("cutaneous", "cutaneous")
  , cloneCount=c(0, 0)
  , cloneFraction=c(0, 0)
  , aaSeqCDR3=c("AAA", "AAC")
  , targetSequences=c("AAA", "AAC")
  , bestVHit=c("", "")
  , bestDHit=c("", "")
  , bestJHit=c("", "")
  , bestCHit=c("", "")
  , bestVGene=c("", "")
  , bestDGene=c("", "")
  , bestJGene=c("", "")
  , bestCGene=c("", "")
  , bestVFamily=c("", "")
  , bestDFamily=c("", "")
  , bestJFamily=c("", "")
  , bestCFamily=c("", "")
  , ID=c("AAA", "AAC")
  , Type=c("BCR", "TCR")
))

table(df_mixcr$Sample, df_mixcr$Type)


df_mixcr$Chain <- substring(df_mixcr$bestVFamily, 1, 3)
dummy <- df_mixcr[df_mixcr$Response=="Good" & df_mixcr$Chain!="", ]
my.col <- RColorBrewer::brewer.pal(n=7, name = "Set2")
names(my.col) <- unique(dummy$Chain)
dummy$Sample <- paste("Good", dummy$Sample, sep = "-")

my_theme <-   theme(
  axis.title.x=element_text(size=20)
  , axis.title.y=element_text(size=20)
  , axis.text.y=element_text(size=16)
  , axis.text.x=element_text(size=16)
  , legend.text=element_text(size=16)
  , plot.title=element_text(size=20, face = "bold")
  , plot.subtitle=element_text(size=0, face = "plain")
  , legend.title=element_text(size=16, face = "bold"), 
)


Figure.S6A <- ggplot(dummy)+geom_bar(aes(x=Sample, fill=Chain), position="fill")+scale_fill_manual(values = my.col)+theme_bw()+coord_flip()+
  ylab("Frequency")+xlab("")+my_theme

dummy2 <- as.data.frame(table(dummy$Sample, dummy$Chain)/rowSums(table(dummy$Sample, dummy$Chain)))

dummy2 <- dummy2[dummy2$Var2%in%c("IGK", "IGH", "IGL"), ]

dummy2$Chain <- dummy2$Var2

my_theme <-   theme(
  axis.title.x=element_text(size=20)
  , axis.title.y=element_text(size=20)
  , axis.text.y=element_text(size=16)
  , axis.text.x=element_text(size=16)
  , legend.text=element_text(size=16)
  , plot.title=element_text(size=20, face = "bold")
  , plot.subtitle=element_text(size=0, face = "plain")
  , legend.title=element_text(size=16, face = "bold") 
  , legend.position = "none", 
)


my_comparisons <- list( c("IGK", "IGH"), c("IGH", "IGL"), c("IGK", "IGL") )
Figure.S6B <- ggviolin(dummy2, x = "Chain", y = "Freq",
         color = "Chain", palette = "jco", add = c("boxplot", "jitter"), caption= "Kruskal-Wallis, p = 4.302e-05", font.caption=14)+ 
  stat_compare_means(comparisons = my_comparisons, size=4.5)+ # Add pairwise comparisons p-value
  ylab("Frequency")+my_theme



# Wilcox Analysis 

dummy <- aggregate(df_mixcr$cloneCount, by=list(Sample=df_mixcr$Sample), FUN=sum)
rownames(dummy) <- dummy$Sample
dummy$Response <- df_mixcr$Response[match(dummy$Sample, df_mixcr$Sample)]
colnames(dummy) <- c("Sample", "cloneCount", "Response")



imk.immdata <- list()
imk.immdata[['data']] <- list()
for(i in unique(df_mixcr$Sample[df_mixcr$Melanoma.type=="cutaneous"])){
  dummy <- df_mixcr[df_mixcr$Sample%in%i, ]
  imk.immdata[['data']][[i]] <- data.frame(Clones=dummy$cloneCount, 
                                           Proportion=dummy$cloneCount/sum(dummy$cloneCount), 
                                           CDR3.nt=dummy$targetSequences, 
                                           CDR3.aa=dummy$aaSeqCDR3, 
                                           V.name=dummy$bestVGene, 
                                           D.name=dummy$bestDGene, 
                                           J.name=dummy$bestJGene, 
                                           V.end=NA,
                                           D.start=NA,
                                           D.end=NA,
                                           J.start=NA,
                                           VJ.ins=NA,
                                           VD.ins=NA,
                                           DJ.ins=NA,
                                           Sequence=NA
  )
}

imk.immdata.tcr <- list()
imk.immdata.tcr[['data']] <- list()
for(i in unique(df_mixcr$Sample[df_mixcr$Melanoma.type=="cutaneous"])){
  dummy <- df_mixcr[df_mixcr$Sample%in%i & df_mixcr$Type=="TCR", ]
  imk.immdata.tcr[['data']][[i]] <- data.frame(Clones=dummy$cloneCount, 
                                               Proportion=dummy$cloneCount/sum(dummy$cloneCount), 
                                               CDR3.nt=dummy$targetSequences, 
                                               CDR3.aa=dummy$aaSeqCDR3, 
                                               V.name=dummy$bestVGene, 
                                               D.name=dummy$bestDGene, 
                                               J.name=dummy$bestJGene, 
                                               V.end=NA,
                                               D.start=NA,
                                               D.end=NA,
                                               J.start=NA,
                                               VJ.ins=NA,
                                               VD.ins=NA,
                                               DJ.ins=NA,
                                               Sequence=NA
  )
}

imk.immdata.bcr <- list()
imk.immdata.bcr[['data']] <- list()
for(i in unique(df_mixcr$Sample[df_mixcr$Melanoma.type=="cutaneous"])){
  dummy <- df_mixcr[df_mixcr$Sample%in%i & df_mixcr$Type=="BCR", ]
  imk.immdata.bcr[['data']][[i]] <- data.frame(Clones=dummy$cloneCount, 
                                               Proportion=dummy$cloneCount/sum(dummy$cloneCount), 
                                               CDR3.nt=dummy$targetSequences, 
                                               CDR3.aa=dummy$aaSeqCDR3, 
                                               V.name=dummy$bestVGene, 
                                               D.name=dummy$bestDGene, 
                                               J.name=dummy$bestJGene, 
                                               V.end=NA,
                                               D.start=NA,
                                               D.end=NA,
                                               J.start=NA,
                                               VJ.ins=NA,
                                               VD.ins=NA,
                                               DJ.ins=NA,
                                               Sequence=NA
  )
}



imk.immdata[["meta"]] <- unique(df_mixcr[, c("Sample", "Response")])
imk.immdata.tcr[["meta"]] <- unique(df_mixcr[, c("Sample", "Response")])
imk.immdata.bcr[["meta"]] <- unique(df_mixcr[, c("Sample", "Response")])

exp_vol <- repExplore(imk.immdata$data, .method = "volume")
exp_vol.tcr <- repExplore(imk.immdata.tcr$data, .method = "volume")
exp_vol.bcr <- repExplore(imk.immdata.bcr$data, .method = "volume")
imm_top <- repClonality(imk.immdata$data, .method = "top", .head = c(5, 20, 100, 500, 5000))

rownames(imk.immdata$meta) <- imk.immdata$meta$Sample
rownames(imm_top) <- paste(imk.immdata$meta[rownames(imm_top), "Response"], rownames(imm_top), sep = "-")


Figure.5B <- vis(imm_top)+xlab("")+theme(
  axis.title.x=element_text(size=20)
  , axis.title.y=element_text(size=20)
  , axis.text.y=element_text(size=16)
  , axis.text.x=element_text(size=16)
  , legend.text=element_text(size=16)
  , plot.title=element_text(size=20, face = "bold")
  , plot.subtitle=element_text(size=16, face = "plain")
  , legend.title=element_text(size=16, face = "bold")
)


my_theme <-   theme(
  axis.title.x=element_text(size=20)
  , axis.title.y=element_text(size=20)
  , axis.text.y=element_text(size=16)
  , axis.text.x=element_text(size=16, angle = 0, vjust = 0.5, hjust = 0.5)
  , legend.text=element_text(size=16)
  , plot.title=element_text(size=20, face = "bold")
  , plot.subtitle=element_text(size=0, face = "plain")
  , legend.title=element_text(size=16, face = "bold"), 
)

Figure.5A <- vis(exp_vol, .by = c("Response"), .meta = imk.immdata$meta) + 
  scale_colour_manual(values =c(Good="#EFC000FF", Bad="#0073C2FF"), aesthetics = c("colour", "fill"))+xlab("")+my_theme+
  vis(exp_vol.bcr, .by = c("Response"), .meta = imk.immdata.bcr$meta) + 
  ggtitle("Number of BCR clonotypes", subtitle = "") + scale_colour_manual(values =c(Good="#EFC000FF", Bad="#0073C2FF"), aesthetics = c("colour", "fill"))+xlab("")+my_theme+
  vis(exp_vol.tcr, .by = c("Response"), .meta = imk.immdata.tcr$meta) + 
  ggtitle("Number of TCR clonotypes", subtitle = "") + scale_colour_manual(values =c(Good="#EFC000FF", Bad="#0073C2FF"), aesthetics = c("colour", "fill"))+xlab("")+my_theme


# D50
div_d50 <- repDiversity(imk.immdata.bcr$data, "d50")
rownames(imk.immdata.bcr$meta) <- imk.immdata.bcr$meta$Sample
#rownames(div_d50) <- paste(imk.immdata.bcr$meta[rownames(div_d50), "Response"], rownames(div_d50), sep = "-")

vis(div_d50, .by = c("Response"), .meta = imk.immdata.bcr$meta)


# Differences BCR - Table 1

target.idt <- names(which(table(df_mixcr$targetSequences)>3))

df_mixcr[df_mixcr$targetSequences%in%target.idt, ]
df_mixcr$ID_clone_gene <- paste(df_mixcr$bestVGene, df_mixcr$bestDGene, df_mixcr$bestJGene, df_mixcr$bestCGene, sep = "-")
df_mixcr$ID_clone_family <- paste(df_mixcr$bestVFamily, df_mixcr$bestDFamily, df_mixcr$bestJFamily, df_mixcr$bestCFamily, sep = "-")

dummy <- aggregate(df_mixcr$cloneCount, by=list(Response=df_mixcr$Response, ID_clone_gene=df_mixcr$ID_clone_gene), FUN=sum)
dummy.good <- dummy[dummy$Response=="Good", ]
dummy.bad <- dummy[dummy$Response=="Bad", ]

top.good <- dummy.good[order(dummy.good$x, decreasing = T), "ID_clone_gene"][1:5]
dummy[dummy$ID_clone_gene%in%top.good, ]


# Circle Plot 

my.col <- unique(df_mixcr[, c("Sample", "Response")])
my.col$col <- ifelse(my.col$Response=="Good", "#eec12b", "#1174bf")
rownames(my.col) <- my.col$Sample


df_mixcr$bestVFamily <- gsub("\\-.*", "", df_mixcr$bestVGene )
df_mixcr$bestDFamily <- gsub("\\-.*", "", df_mixcr$bestDGene )
df_mixcr$bestJFamily <- gsub("\\-.*", "", df_mixcr$bestJGene )
df_mixcr$bestCFamily <- gsub("\\-.*", "", df_mixcr$bestCGene )

df_mixcr$Type <- NA
df_mixcr[grep("TR", df_mixcr$bestVFamily), "Type"] <- "TCR"
df_mixcr[grep("IG", df_mixcr$bestVFamily), "Type"] <- "BCR"

# bestVFamily 

dummy <- aggregate(df_mixcr$cloneCount, by=list(Sample=df_mixcr$Sample, Type=df_mixcr$Type), FUN=sum)
dummy$Response <- df_mixcr$Response[match(dummy$Sample, df_mixcr$Sample)]
colnames(dummy) <- c("Sample", "Type", "cloneCount", "Response")


rownames(my.col) <- paste(my.col$Response, my.col$Sample, sep = "-")

df_mixcr$Sample <- paste(df_mixcr$Response, df_mixcr$Sample, sep = "-")

dummy <- aggregate(df_mixcr$cloneCount, by=list(Sample=df_mixcr$Sample, bestVFamily=df_mixcr$bestVFamily), FUN=sum)
dummy <- dummy[dummy$x>10, ]
dummy <- dummy[dummy$Sample%in%names(table(dummy$Sample))[table(dummy$Sample)>1], ]

grid.col <- union(dummy$Sample, dummy$bestVFamily)
names(grid.col) <- grid.col

idt <- names(grid.col)[names(grid.col)%in%rownames(my.col)]
grid.col[idt] <- my.col[idt,"col"]
grid.col[-grep("IMK", names(grid.col))] <- "lightseagreen"

par(cex = 1, mar = c(0, 0, 0, 0))
chordDiagram(dummy, annotationTrack = "grid", preAllocateTracks = list(track.height = 0.1), grid.col = grid.col)


circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  if(abs(xplot[2] - xplot[1]) > 10) {
    circos.genomicText(data.frame(mean(xlim), mean(ylim)), value="hola", labels=("adios"), sector.name, facing = "inside",
                       niceFacing = TRUE, adj = c(0.5, 0))
    circos.axis("bottom", labels.cex = 0.4)
  }
}, bg.border = NA)


circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  if(abs(xplot[2] - xplot[1]) > 10) {
    circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
                niceFacing = TRUE, adj = c(0, 0.5), col = ifelse(sector.name %in% rownames(my.col), my.col[sector.name, "col"], "lightseagreen"))

    
  }
  
}, bg.border = NA)



# bestDFamily 

dummy <- aggregate(df_mixcr$cloneCount, by=list(Sample=df_mixcr$Sample, bestDFamily=df_mixcr$bestDFamily), FUN=sum)
dummy <- dummy[dummy$x>5, ]
dummy <- dummy[dummy$bestDFamily!="", ]


grid.col <- union(dummy$Sample, dummy$bestDFamily)
names(grid.col) <- grid.col

idt <- names(grid.col)[names(grid.col)%in%rownames(my.col)]
grid.col[idt] <- my.col[idt,"col"]
grid.col[-grep("IMK", names(grid.col))] <- "lightsalmon"

chordDiagram(dummy, annotationTrack = "grid", preAllocateTracks = list(track.height = 0.1), grid.col = grid.col)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  if(abs(xplot[2] - xplot[1]) < 10) {
    circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
                niceFacing = TRUE, adj = c(0, 0.5), col = ifelse(sector.name %in% rownames(my.col), my.col[sector.name, "col"], "lightsalmon"))
  } else {
    circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
                niceFacing = TRUE, adj = c(0, 0.5), col = ifelse(sector.name %in% rownames(my.col), my.col[sector.name, "col"], "lightsalmon"))
  }
  
}, bg.border = NA)

# bestJFamily 

dummy <- aggregate(df_mixcr$cloneCount, by=list(Sample=df_mixcr$Sample, bestJFamily=df_mixcr$bestJFamily), FUN=sum)
dummy <- dummy[dummy$x>10, ]
dummy <- dummy[dummy$Sample%in%names(table(dummy$Sample))[table(dummy$Sample)>1], ]


grid.col <- union(dummy$Sample, dummy$bestJFamily)
names(grid.col) <- grid.col

idt <- names(grid.col)[names(grid.col)%in%rownames(my.col)]
grid.col[idt] <- my.col[idt,"col"]
grid.col[-grep("IMK", names(grid.col))] <- "lightcoral"

chordDiagram(dummy, annotationTrack = "grid", preAllocateTracks = list(track.height = 0.1), grid.col = grid.col)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  if(abs(xplot[2] - xplot[1]) < 10) {
    circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
                niceFacing = TRUE, adj = c(0, 0.5), col = ifelse(sector.name %in% rownames(my.col), my.col[sector.name, "col"], "lightcoral"))
  } else {
    circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
                niceFacing = TRUE, adj = c(0, 0.5), col = ifelse(sector.name %in% rownames(my.col), my.col[sector.name, "col"], "lightcoral"))
  }
  
}, bg.border = NA)

# bestCFamily 
dummy <- aggregate(df_mixcr$cloneCount, by=list(Sample=df_mixcr$Sample, bestCFamily=df_mixcr$bestCFamily), FUN=sum)
dummy <- dummy[dummy$x>10, ]
dummy <- dummy[dummy$Sample%in%names(table(dummy$Sample))[table(dummy$Sample)>1], ]
dummy <- dummy[dummy$bestCFamily!="", ]

grid.col <- union(dummy$Sample, dummy$bestCFamily)
names(grid.col) <- grid.col

idt <- names(grid.col)[names(grid.col)%in%rownames(my.col)]
grid.col[idt] <- my.col[idt,"col"]
grid.col[-grep("IMK", names(grid.col))] <- "yellowgreen"

chordDiagram(dummy, annotationTrack = "grid", preAllocateTracks = list(track.height = 0.1), grid.col = grid.col)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  #if(abs(xplot[2] - xplot[1]) < 10) {
  circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
              niceFacing = TRUE, adj = c(0, 0.5), col = ifelse(sector.name %in% rownames(my.col), my.col[sector.name, "col"], "yellowgreen"))

  
}, bg.border = NA)


# HLA ---------------------------------------------------------------------



## RPKM
classI <- list.files(path = "../data/seq2HLA_output/",pattern = "-ClassI-class.expression", recursive = T, full.names = T)
nonclass <- list.files(path = "../data/seq2HLA_output/",pattern = "-ClassI-nonclass.expression", recursive = T, full.names = T)
classII <- list.files(path = "../data/seq2HLA_output/",pattern = "-ClassII.expression", recursive = T, full.names = T)

df.classI <- do.call('rbind', lapply(classI, function(file){
  df <- read.table(file, sep = ":", dec = ".", header = F)
  df[, 1] <- gsub(" ", "", df[, 1])
  df[, 2] <- gsub(" ", "", df[, 2])
  df[, 2] <- gsub("RPKM", "", df[, 2])
  colnames(df) <- c("Allele", "RPKM")
  df$Sample <- strsplit(file, "/")[[1]][[4]]
  df$Class <- rep("Class I", nrow(df))
  return(df)
}))

df.nonclass <- do.call('rbind', lapply(nonclass, function(file){
  df <- read.table(file, sep = ":", dec = ".", header = F)
  df[, 1] <- gsub(" ", "", df[, 1])
  df[, 2] <- gsub(" ", "", df[, 2])
  df[, 2] <- gsub("RPKM", "", df[, 2])
  colnames(df) <- c("Allele", "RPKM")
  df$Sample <- strsplit(file, "/")[[1]][[4]]
  df$Class <- rep("Non-Class I", nrow(df))
  return(df)
}))

df.classII <- do.call('rbind', lapply(classII, function(file){
  df <- read.table(file, sep = ":", dec = ".", header = F)
  df[, 1] <- gsub(" ", "", df[, 1])
  df[, 2] <- gsub(" ", "", df[, 2])
  df[, 2] <- gsub("RPKM", "", df[, 2])
  colnames(df) <- c("Allele", "RPKM")
  df$Sample <- strsplit(file, "/")[[1]][[4]]
  df$Class <- rep("Class II", nrow(df))
  return(df)
}))

df.all <- rbind(df.classI, df.nonclass)
df.all <- rbind(df.all, df.classII)

df.rpkm <- df.all

## HLA genotype

classI <- list.files(path = "../data/seq2HLA_output/", pattern = "-ClassI-class.HLAgenotype4digits", recursive = T, full.names = T)
nonclass <- list.files(path = "../data/seq2HLA_output/", pattern = "-ClassI-nonclass.HLAgenotype4digits", recursive = T, full.names = T)
classII <- list.files(path = "../data/seq2HLA_output/", pattern = "ClassII.HLAgenotype4digits", recursive = T, full.names = T)


df.classI <- do.call('rbind', lapply(classI, function(file){
  df <- read.table(file, sep = "\t", dec = ".", header = T, quote="\"")
  colnames(df) <- c('Locus',	'Allele-1',	'Confidence-1',	'Allele-2', 'Confidence-2')
  df$Sample <- strsplit(file, "/")[[1]][[4]]
  df$Class <- rep("Class I", nrow(df))
  return(df)
}))

df.nonclass <- do.call('rbind', lapply(nonclass, function(file){
  df <- read.table(file, sep = "\t", dec = ".", header = T, quote="\"")
  colnames(df) <- c('Locus',	'Allele-1',	'Confidence-1',	'Allele-2', 'Confidence-2')
  df$Sample <- strsplit(file, "/")[[1]][[4]]
  df$Class <- rep("Non-Class I", nrow(df))
  return(df)
}))

df.classII <- do.call('rbind', lapply(classII, function(file){
  df <- read.table(file, sep = "\t", dec = ".", header = T, quote="\"")
  colnames(df) <- c('Locus',	'Allele-1',	'Confidence-1',	'Allele-2', 'Confidence-2')
  df$Sample <- strsplit(file, "/")[[1]][[4]]
  df$Class <- rep("Class II", nrow(df))
  return(df)
}))

df.all <- rbind(df.classI, df.nonclass)
df.all <- rbind(df.all, df.classII)

df.confidence <- df.all


samples_data <- read.delim("../data/clinical.txt")
samples_data <- samples_data[samples_data$Diagnosis=="cutaneous", ]

## remove imk_worRNA form rpkm
df.rpkm$Sample <- gsub("_S[[:digit:]]_worRNA", "", df.rpkm$Sample)
df.rpkm$Sample <- gsub("IMK0", "IMK", df.rpkm$Sample)
df.rpkm <- df.rpkm[df.rpkm$Sample%in%rownames(samples_data), ]

# Batch inicial

df.rpkm$Response <- samples_data$Response[match(df.rpkm$Sample, rownames(samples_data))]

## remove imk_worRNA form confidence
df.confidence$Sample <- gsub("_S[[:digit:]]_worRNA", "", df.confidence$Sample)
df.confidence$Sample <- gsub("IMK0", "IMK", df.confidence$Sample)

df.confidence <- df.confidence[df.confidence$Sample%in%rownames(samples_data), ]
df.confidence$Response <- samples_data$Response[match(df.confidence$Sample, rownames(samples_data))]


df.rpkm.2 <- dcast(df.rpkm[, c(1, 2, 3)], Sample ~ Allele, value.var = "RPKM")
df.rpkm.2 <- as.data.frame(df.rpkm.2)
rownames(df.rpkm.2) <- df.rpkm.2$Sample
df.rpkm.2$Sample <- NULL
samples <- rownames(df.rpkm.2)
df.rpkm.2 <- apply(df.rpkm.2, 2, as.numeric)
rownames(df.rpkm.2) <- samples

df.rpkm$RPKM <- as.numeric(df.rpkm$RPKM)


# Analysis 

wilcox.test(RPKM ~ Response, df.rpkm)

my_comparisons <- list( c("Bad", "Good") )

a <- ggboxplot(df.rpkm, x = "Response", y = "RPKM",
               color = "Response", palette = "jco", add = "jitter", main="HLA")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", label.y = 100)+ #
  stat_compare_means(label.y = 140)

my_theme <-   theme(
  axis.title.x=element_text(size=20)
  , axis.title.y=element_text(size=20)
  , axis.text.y=element_text(size=16)
  , axis.text.x=element_text(size=16)
  , legend.text=element_text(size=16)
  , plot.title=element_text(size=20, face = "bold",  hjust = 0.5)
  , plot.subtitle=element_text(size=0, face = "plain")
  , legend.title=element_text(size=16, face = "bold")
  , legend.position = "none"
)


b <- ggboxplot(df.rpkm, x = "Response", y = "RPKM", facet.by = "Class",
               color = "Response", palette = "jco", add = "jitter", group.by = "facet", main = "HLA")+
  stat_compare_means(comparisons = my_comparisons, label.y = 100)+my_theme+xlab("") #


Figure.5C <- ggarrange(a, b, common.legend = T)


my.col <- unique(df.rpkm[, c("Sample", "Response")])
my.col$col <- ifelse(my.col$Response=="Good", "blue", "red")
rownames(my.col) <- my.col$Sample



# Chord figure
dummy <- aggregate(df.rpkm$RPKM, by=list(Sample=df.rpkm$Sample, Allele=df.rpkm$Allele, Response=df.rpkm$Response), FUN=sum)
dummy$Sample <- paste(dummy$Response, dummy$Sample, sep = "-")
rownames(my.col) <- paste(my.col$Response, my.col$Sample, sep = "-")

dummy <- dummy[!dummy$Allele%in%c("J", "K", "P", "V"), ]


group <- rep(NA, length(unique(c(dummy$Sample, dummy$Allele))))
names(group) <- unique(c(dummy$Sample, dummy$Allele))
group[rownames(samples_data)] <- samples_data$Response.to.IT
group[is.na(group)] <- "Allele" 


grid.col <- union(dummy$Sample, dummy$Allele)
names(grid.col) <- grid.col


my.col[my.col$Response=="Good", "col"] = "#eec12b"
my.col[my.col$Response=="Bad", "col"] = "#1174bf"

idt <- names(grid.col)[names(grid.col)%in%rownames(my.col)]
grid.col[idt] <- my.col[idt,"col"]
grid.col[-grep("IMK", names(grid.col))] <- "maroon"


chordDiagram(dummy, annotationTrack = "grid", preAllocateTracks = list(track.height = 0.1), grid.col = grid.col)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  if(abs(xplot[2] - xplot[1]) < 10) {
    circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
                niceFacing = TRUE, adj = c(0, 0.5), col = ifelse(sector.name %in% rownames(my.col), my.col[sector.name, "col"], "maroon"))
  } else {
    circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
                niceFacing = TRUE, adj = c(0, 0.5), col = ifelse(sector.name %in% rownames(my.col), my.col[sector.name, "col"], "maroon"))
  }
  
}, bg.border = NA)


# Merging HLA, BCR, TCR ----------------------------------------------------------------


hla <- c(brewer.pal(n = 12, name = "Set3"), brewer.pal(n = 2, name = "Set2"))

df.rpkm$ID <- df.rpkm$Sample

net_mixcr <- df_mixcr
net_mixcr$ID <- paste(net_mixcr$Response, net_mixcr$Sample, sep = "-")

col <- c(brewer.pal(n = 12, name = "Set3"), brewer.pal(n = 8, name = "Set2"), brewer.pal(n = 8, name = "Set1"))

p1 <- ggplot(data=df.rpkm[df.rpkm$Sample%in%c("IMK35", "IMK37", "IMK39"), ]) + aes(x=Allele, y=RPKM, color = Allele, fill=Allele) + geom_bar(stat = "identity")+
  coord_polar(theta = "y")+theme_bw()+
  scale_y_log10(limits = c(1,150)) + scale_color_manual(values = col) + scale_fill_manual(values = col) + 
  facet_grid(.~ID) + ylab("HLA")+ xlab("RPKM") + theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 0, colour = "black"))+theme(strip.text.x = element_text(size = 12))


p2 <- ggplot(data=df.rpkm[df.rpkm$Sample%in%c("IMK26", "IMK31", "IMK38"), ]) + aes(x=Allele, y=RPKM, color = Allele, fill=Allele) + geom_bar(stat = "identity")+
  coord_polar(theta = "y")+theme_bw()+
  scale_y_log10(limits = c(1,150)) + scale_color_manual(values = col) + scale_fill_manual(values = col) + 
  facet_grid(.~ID) + ylab("HLA")+ xlab("RPKM") + theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 0, colour = "black"))+theme(strip.text.x = element_text(size = 12))


p3 <- ggplot(data=net_mixcr[net_mixcr$Sample%in%c("IMK26", "IMK31", "IMK38") & net_mixcr$Type=="TCR" & net_mixcr$cloneCount>1, ]) + aes(x=bestVFamily, y=cloneCount, color = bestVFamily, fill=bestVFamily) + geom_bar(stat = "identity")+
  coord_polar(theta = "y")+theme_bw()+
  scale_y_log10(limits = c(1,1500))+ scale_color_manual(values = col) + scale_fill_manual(values = col) + 
  facet_grid(.~ID) + ylab("TCR")+ xlab("VFamily") + theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 0, colour = "black"))+theme(strip.text.x = element_text(size = 12))

p4 <- ggplot(data=net_mixcr[net_mixcr$Sample%in%c("IMK26", "IMK31", "IMK38") & net_mixcr$Type=="BCR" & net_mixcr$cloneCount>2, ]) + aes(x=bestVFamily, y=cloneCount, color = bestVFamily, fill=bestVFamily) + geom_bar(stat = "identity")+
  coord_polar(theta = "y")+theme_bw()+
  scale_y_log10(limits = c(1,1500))+ scale_color_manual(values = col) + scale_fill_manual(values = col) +
  facet_grid(.~ID) + ylab("BCR")+ xlab("VFamily") + theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 0, colour = "black"))+theme(strip.text.x = element_text(size = 12))

p5 <- ggplot(data=net_mixcr[net_mixcr$Sample%in%c("IMK35", "IMK37", "IMK39") & net_mixcr$Type=="TCR" & net_mixcr$cloneCount>0, ]) + aes(x=bestVFamily, y=cloneCount, color = bestVFamily, fill=bestVFamily) + geom_bar(stat = "identity")+
  coord_polar(theta = "y")+theme_bw()+
  scale_y_log10(limits = c(1,1500))+ scale_color_manual(values = col) + scale_fill_manual(values = col) + 
  facet_grid(.~ID) + ylab("TCR")+ xlab("VFamily") + theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 0, colour = "black"))+theme(strip.text.x = element_text(size = 12))

p6 <- ggplot(data=net_mixcr[net_mixcr$Sample%in%c("IMK35", "IMK37", "IMK39") & net_mixcr$Type=="BCR" & net_mixcr$cloneCount>0, ]) + aes(x=bestVFamily, y=cloneCount, color = bestVFamily, fill=bestVFamily) + geom_bar(stat = "identity")+
  coord_polar(theta = "y")+theme_bw()+
  scale_y_log10(limits = c(1,1500))+ scale_color_manual(values = col) + scale_fill_manual(values = col) +
  facet_grid(.~ID) + ylab("BCR")+ xlab("VFamily") + theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 0, colour = "black"))+theme(strip.text.x = element_text(size = 12))



Figure.5D1 <- grid.arrange(p2, p3, p4, nrow = 3)
Figure.5D1 <- grid.arrange(p1, p5, p6, nrow = 3)


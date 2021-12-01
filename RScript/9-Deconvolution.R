###################################################
## 
## Deconvolution of RNA-seq
##
## - MCPcounter
## - CIBERSORTx
## - Multiespectral
## - Correlation between Multiespectral-CIBERSORTx
##
## Manuscript --> Figure 3C-D, Figure 4C, Table S7, Table S9
##
###################################################



# MCPcounter --------------------------------------------------------------

library(reshape)
library(ggplot2)
library(dplyr)
library(tidyr)
library(MCPcounter)
library(RColorBrewer)

samples_data <- read.delim("../data/clinical.txt")
rownames(samples_data) <- paste0(samples_data$Response, ".", rownames(samples_data))

TPM <- read.table("../data/TPM_IMK_gene.txt", sep = "\t", stringsAsFactors = F, header = T)
TPM <- TPM[!duplicated(TPM$SYMBOL), ]
rownames(TPM) <- TPM$SYMBOL
TPM$SYMBOL <- NULL
IMK_ids <- rownames(samples_data)[samples_data$Diagnosis%in%"cutaneous"]
TPM <- TPM[, IMK_ids]



mcp.results_tpm <- MCPcounter.estimate(expression = TPM, featuresType = "HUGO_symbols")
genes=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"),sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE)

mcp.counter <- reshape2::melt(mcp.results_tpm)
colnames(mcp.counter) <- c("CellType", "Sample", "Score")
mcp.counter$Sample <- as.character(mcp.counter$Sample)
mcp.counter <- mcp.counter[order(mcp.counter$Sample), ]


col <- brewer.pal(n = 10, name = 'Set3')
names(col) <- unique(mcp.counter$CellType)

mcp.counter$Response <- gsub(".IMK[[:digit:]][[:digit:]]", "", mcp.counter$Sample)
mcp.counter$Sample <- gsub("[.]", "-", mcp.counter$Sample)
mcp.counter <- mcp.counter[order(mcp.counter$Response, decreasing = F), ]
mcp.counter$Sample <- factor(mcp.counter$Sample, levels = unique(mcp.counter$Sample))
ggplot(data = mcp.counter) + theme_bw() +
  geom_bar(stat="identity", width = 0.9, aes(y=Score, x=Sample, fill=CellType)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=16)
        , axis.title.x=element_text(size=20)
        , axis.title.y=element_text(size=20)
        , axis.text.y=element_text(size=16)
        , plot.title=element_text(size=20, face = "bold", hjust = 0.5)
        , legend.text=element_text(size=16)
        , legend.title=element_text(size=16, face = "bold")) + 
  scale_fill_manual("", values =  col[names(col) %in% unique(mcp.counter$CellType)]) +
  xlab("") + ylab("Estimates of MCP-counter")+ggtitle("MCP-counter")

mcp.res <- list()
for(i in levels(mcp.counter$CellType)){
  mcp.res[[i]] <- wilcox.test(Score~Response, mcp.counter[mcp.counter$CellType%in%i, ])$p.value
}
mcp.res


# CIBERSORTx --------------------------------------------------------------

cibersort <- read.delim("../data/CIBERSORTx.txt")

cut.idt <- c("Bad.IMK21" , "Bad.IMK22" , "Good.IMK23", "Good.IMK24", "Good.IMK26",
             "Good.IMK27", "Good.IMK30", "Good.IMK31", "Good.IMK32", "Bad.IMK34" ,
             "Bad.IMK35" , "Bad.IMK36" , "Bad.IMK37" , "Good.IMK38", "Bad.IMK39" ,
             "Good.IMK48")

cibersort <- reshape2::melt(cibersort)

cibersort$Response <- gsub(".IMK[[:digit:]][[:digit:]]", "", cibersort$Mixture)
cibersort[, 2] <- as.character(cibersort[, 2])
cibersort <- cibersort[cibersort[, 1]%in%c(cut.idt), ]

cibersort <- cibersort[!cibersort[,2]%in%c("P.value", "Correlation", "RMSE", "Absolute.score..sig.score."), ]

library(ggpubr)
my.list <- list()

my_comparisons <- list( c("Bad", "Good") )
for(i in unique(cibersort[, 2])){
  if(any(cibersort[cibersort[,2]==i, 3]!=0)){
    p <- ggviolin(cibersort[cibersort[,2]==i, ], x = "Response", y = "value", facet.by = "variable",
                  color = "Response", palette = "jco", add = "dotplot", group.by = "variable")+
      stat_compare_means(comparisons = my_comparisons, label = "p.signif", label.y = max(cibersort[cibersort[,2]==i, 3])*1.2)+ #
      stat_compare_means()
    my.list[[i]] <- p
  }
}

ggpubr::ggarrange(plotlist=my.list, common.legend=T)

cibersort$Mixture <- gsub("[.]", "-", cibersort$Mixture)
cibersort$variable <- gsub("[.]", " ", cibersort$variable)
length(unique(cibersort$variable))
col <- c(brewer.pal(n = 8, name = "Set1"), brewer.pal(n = 8, name = "Set3"))
names(col) <- unique(cibersort$variable)

ggplot(data = cibersort, aes(y=value, x=Mixture, fill=variable)) +
  geom_bar(stat="identity") + theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1, size=16)
                                                , axis.title.x=element_text(size=20)
                                                , axis.title.y=element_text(size=20)
                                                , axis.text.y=element_text(size=16)
                                                , legend.text=element_text(size=16)
                                                , plot.title=element_text(size=20, face = "bold", hjust = 0.5)
                                                , legend.title=element_text(size=16, face = "bold")) + xlab("") + ylab("Estimates of CIBERSORTx")+
  scale_fill_manual("", values =  col[names(col) %in% unique(cibersort$variable)])+ggtitle("CIBERSORTx")


p.res_mcpcounter <- list()
for(i in unique(mcp.counter[, 1])){
  if(any(mcp.counter[mcp.counter[,1]==i, 3]!=0)){
    p.res_mcpcounter[[i]] <- wilcox.test(Score~Response, mcp.counter[mcp.counter[,1]==i, ])$p.value
  }
}

p.res_cibersort <- list()
for(i in unique(cibersort[, 2])){
  if(any(cibersort[cibersort[,2]==i, 3]!=0)){
    p.res_cibersort[[i]] <- wilcox.test(value~Response, cibersort[cibersort[,2]==i, ])$p.value
  }
}

# Comparasion MCP-counter - CIBERSORTx -------------------------------------------------------------

cibersort.score <- aggregate(cibersort$value, by=list(Sample=cibersort$Mixture), FUN=sum)
mcp.counter.score <- aggregate(mcp.counter$Score, by=list(Sample=mcp.counter$Sample), FUN=sum)

colnames(cibersort.score)[2] <- c("CIBERSORTx")
colnames(mcp.counter.score)[2] <- c("MCPcounter")

rownames(cibersort.score) <- cibersort.score$Sample
rownames(mcp.counter.score) <- mcp.counter.score$Sample

mcp.counter.score$Samples <- mcp.counter.score$Sample
idt <- rownames(mcp.counter.score)

deconv <- as.data.frame(cbind(cibersort.score[idt, c("Sample","CIBERSORTx")], mcp.counter.score[idt,  c("Samples", "MCPcounter")]))

ggplot(deconv, aes(x=CIBERSORTx, y=MCPcounter)) + 
  geom_point(aes(col = Samples))+
  geom_smooth(method=lm)+xlab("CIBERSORTx - Score")+ylab("MCP counter - Score")+theme_bw()

cor(deconv$CIBERSORTx, deconv$MCPcounter)
cor.test(deconv$CIBERSORTx, deconv$MCPcounter)

df.cibersort <- reshape2::dcast(Mixture~variable, value.var = "value",data = cibersort)
rownames(df.cibersort) <- df.cibersort$Mixture
df.cibersort$Mixture <- NULL

df.mcp <- as.data.frame(t(mcp.results_tpm))
rownames(df.mcp) <- gsub("[.]", "-", rownames(df.mcp))


res <- matrix(NA, nrow = ncol(df.cibersort), ncol = ncol(df.mcp))
rownames(res) <- colnames(df.cibersort)
colnames(res) <- colnames(df.mcp)

p.res <- res

idt <- intersect(rownames(df.cibersort), rownames(df.mcp))

for(i in rownames(res)){
  for(j in colnames(res)){
    res[i, j] <- cor(df.cibersort[idt, i], df.mcp[idt, j])
    p.res[i, j] <- cor.test(df.cibersort[idt, i], df.mcp[idt, j])$p.value
  }
}

res <- na.omit(res)
p.res <- na.omit(p.res)

my.palette <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", 
                                 "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
                                 "#4393C3", "#2166AC", "#053061"))(200)
my.palette <- rev(my.palette)


corrplot::corrplot(na.omit(res), p.mat = na.omit(p.res), insig = "label_sig", pch.cex = 2, col = my.palette)

# Multiespectral ----------------------------------------------------------

multiespectral <- read.delim("/Volumes/TRANSCEND/IBIMA/Paper/data/multiespectral.txt")
samples_data <- read.delim("../data/clinical.txt")

rownames(multiespectral) <- multiespectral$Samples
multiespectral$Response <- samples_data[multiespectral$Samples, "Response"]
multiespectral$Melanoma.type <- samples_data[multiespectral$Samples, "Diagnosis"]
multiespectral$Batch_used <- samples_data[multiespectral$Samples, "Batch"]
multiespectral$Samples <- NULL

## Keep only Batch1-Batch4 and cutaneous
multiespectral <- multiespectral[multiespectral$Batch_used!="Batch5", ]
multiespectral$Melanoma.type <- NULL
multiespectral$Batch_used <- NULL

# Wilcox test
df.multiespectral <- multiespectral
multiespectral <- reshape2::melt(multiespectral)
res.test <- list()

for(var in unique(multiespectral$variable)){
  dummy <- multiespectral[multiespectral$variable%in%var & !is.na(multiespectral$value), ]
  res.test[[var]] <- wilcox.test(value~Response, data = dummy)$p.value
}

res.test[unlist(res.test)<0.06]

dummy <- multiespectral[multiespectral$variable%in%c("CD19.dens", "dens_CD8_Melanoma") & !is.na(multiespectral$value), ]
library(ggpubr)
my_comparisons <- list( c("Bad", "Good") )

dummy$variable <- factor(as.character(dummy$variable))
levels(dummy$variable) <- c("CD19", "CD8")

my_theme <-   theme(
  axis.title.x=element_text(size=20)
  , axis.title.y=element_text(size=20)
  , axis.text.y=element_text(size=16)
  , axis.text.x=element_text(size=16)
  , legend.text=element_text(size=16)
  , strip.text.x=element_text(size=16)
  , plot.title=element_text(size=20, face = "bold")
  , plot.subtitle=element_text(size=0, face = "plain")
  , legend.title=element_text(size=16, face = "bold"), 
)

ggboxplot(dummy, x = "Response", y = "value", fill = "Response", facet.by = "variable", palette = "jco", outlier.shape=NA)+
  stat_compare_means(comparisons = my_comparisons, label = "p.format", label.y = 400, size = 5)+xlab("")+ylab(expression("Cells"/ mm^2))+my_theme




# Correlation CIBERSORTx - Multiespectral ---------------------------------

rownames(df.cibersort) <- gsub("Bad-", "", gsub("Good-", "", rownames(df.cibersort)))

res <- matrix(NA, nrow = ncol(df.multiespectral)-1, ncol = ncol(df.cibersort))
rownames(res) <- colnames(df.multiespectral)[-ncol(df.multiespectral)]
colnames(res) <- colnames(df.cibersort)
p.res <- res
idt <- intersect(rownames(df.multiespectral), rownames(df.cibersort))

res.variables.spearman <- c()

for(i in rownames(res)){
  for(j in colnames(res)){
    non.na <- !is.na(df.multiespectral[idt, i])
    res[i, j] <- cor(df.multiespectral[idt, i][non.na], df.cibersort[idt, j][non.na], method = "spearman")
    p.res[i, j] <- cor.test(df.multiespectral[idt, i][non.na], df.cibersort[idt, j][non.na], method = "spearman")$p.value
    if(!is.na(p.res[i, j])) if(p.res[i, j]<0.05) res.variables.spearman <- c(res.variables.spearman, paste(i, j, sep = " - "))
  }
}

res <- res[, colSums(is.na(res))==0]
p.res <- p.res[, colSums(is.na(p.res))==0]

corrplot::corrplot(res, p.mat = p.res)






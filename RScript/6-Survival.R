###################################################
## 
## Survival analysis for gene expression
##
## - Logrank with KM based on quantile expression
## - Univariate cox proportional analysis
##
## Manuscript --> Table S6, Figure 2
##
###################################################

library(survival)
library(survminer)
library(ggplot2)

source("help_fn.R")

clinical <- read.delim("../data/clinical.txt")
rownames(clinical) <- paste0(clinical$Response, "-", rownames(clinical))
clinical <- clinical[clinical$Diagnosis%in%"cutaneous", ]
norm.counts <- read.delim("../data/Result_RNASeq_DESeq2_melanoma_cutaneous.txt")
rownames(norm.counts) <- norm.counts$SYMBOL
colnames(norm.counts) <- gsub("[.]", "-", colnames(norm.counts))
norm.counts <- t(norm.counts[, rownames(clinical)])

quantile(unlist(as.list(norm.counts)), c(.33,.66, .99))
# # #       33%         66%         99% 
# # #   5.649889   61.339891 4389.779549  
event_rna<-norm.counts
event_rna[norm.counts<5.649889]<-"low"
event_rna[norm.counts >= 5.649889   & norm.counts <  61.339891]<-"med"
event_rna[norm.counts>=61.339891]<-"high"


custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5)
      , legend.text  = element_text(size=16)
      , title = element_text(size=16, face = "bold.italic") 
      , axis.text.y = element_text(size=16)
      , axis.text.x = element_text(size=16)
    )
}

# Do KM analysis
df_surv <- cbind(clinical[rownames(event_rna), c("OS", "OS_event", "PFS", "PFS_event")], event_rna)
res_OS <- do_ggsurvival_gene(df_surv = df_surv, time="OS", event = "OS_event")
res_PFS <- do_ggsurvival_gene(df_surv = df_surv, time="PFS", event = "PFS_event")

gene_OS <- names(res_OS)[unlist(res_OS)<0.05]
gene_PFS <- names(res_PFS)[unlist(res_PFS)<0.05]
colnames(df_surv) <- gsub("[-]", ".", colnames(df_surv))
gene_OS <- gsub("[-]", ".", gene_OS)
gene_PFS <- gsub("[-]", ".", gene_PFS)

genes <- unique(c(gene_OS, gene_PFS))

# Do HR for significant gene in logrank test
res_HR_OS <- do_survival(datos = df_surv, vars = genes, time = "OS", event = "OS_event")
res_HR_PFS <- do_survival(datos = df_surv, vars = genes, time = "PFS", event = "PFS_event")

Table.S6 <- data.frame(
  Gene = res_HR_OS$variable, 
  Level.OS = res_HR_OS$level, 
  Logrank.OS = res_HR_OS$Logrank, 
  HR.OS = res_HR_OS$coef, 
  Logrank.PFS = res_HR_PFS$Logrank, 
  HR.PFS = res_HR_PFS$coef
)



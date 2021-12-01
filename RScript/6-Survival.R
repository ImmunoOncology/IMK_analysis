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

my.col  = c("Low"="firebrick2", "Med"="mediumseagreen", "High"="skyblue")

df_surv[, "B4GALNT2"] <- factor(df_surv[, "B4GALNT2"], labels = c("Low", "Med", "High"))
Figure.2.1 <- p <- ggsurvplot(
  survfit(Surv(PFS, PFS_event) ~ B4GALNT2, data = df_surv),                 
  conf.int = F,         
  size=1,#0.7,                    
  pval=TRUE, 
  palette = my.col[levels(df_surv[, "B4GALNT2"])[levels(df_surv[, "B4GALNT2"])%in%unique(df_surv[, "B4GALNT2"])]],
  pval.method=TRUE,   break.x.by = 500,    
  xlab="Time in days",
  ylab="Progression-free survival probability",
  ylim=c(0,1),
  surv.scale="percent",
  legend.labs = levels(df_surv[, "B4GALNT2"])[levels(df_surv[, "B4GALNT2"])%in%unique(df_surv[, "B4GALNT2"])],
  tables.col="strata",
  risk.table=F, 
  risk.table.col = "strata",
  risk.table.y.text = T,
  legend = c(.7,.8),
  tables.y.text = T, 
  legend.title="", ggtheme = custom_theme(), title = "B4GALNT2")

df_surv[, "KCNA1"] <- factor(df_surv[, "KCNA1"], labels = c("Low", "Med", "High"))
Figure.2.2 <- p <- ggsurvplot(
  survfit(Surv(PFS, PFS_event) ~ KCNA1, data = df_surv),                 
  conf.int = F,         
  size=1,#0.7,                    
  pval=TRUE, 
  palette = my.col[levels(df_surv[, "KCNA1"])[levels(df_surv[, "KCNA1"])%in%unique(df_surv[, "KCNA1"])]],
  pval.method=TRUE,   break.x.by = 500,    
  xlab="Time in days",
  ylab="Progression-free survival probability",
  ylim=c(0,1),
  surv.scale="percent",
  legend.labs = levels(df_surv[, "KCNA1"])[levels(df_surv[, "KCNA1"])%in%unique(df_surv[, "KCNA1"])],
  tables.col="strata",
  risk.table=F, 
  risk.table.col = "strata",
  risk.table.y.text = T,
  legend = c(.7,.8),
  tables.y.text = T, 
  legend.title="", ggtheme = custom_theme(), title = "KCNA1")

df_surv[, "LGR5"] <- factor(df_surv[, "LGR5"], labels = c("Low", "Med", "High"))
Figure.2.3 <- p <- ggsurvplot(
  survfit(Surv(PFS, PFS_event) ~ LGR5, data = df_surv),                 
  conf.int = F,         
  size=1,#0.7,                    
  pval=TRUE, 
  palette = my.col[levels(df_surv[, "LGR5"])[levels(df_surv[, "LGR5"])%in%unique(df_surv[, "LGR5"])]],
  pval.method=TRUE,   break.x.by = 500,    
  xlab="Time in days",
  ylab="Progression-free survival probability",
  ylim=c(0,1),
  surv.scale="percent",
  legend.labs = levels(df_surv[, "LGR5"])[levels(df_surv[, "LGR5"])%in%unique(df_surv[, "LGR5"])],
  tables.col="strata",
  risk.table=F, 
  risk.table.col = "strata",
  risk.table.y.text = T,
  legend = c(.7,.8),
  tables.y.text = T, 
  legend.title="", ggtheme = custom_theme(), title = "LGR5")

df_surv[, "IGLV3.21"] <- factor(df_surv[, "IGLV3.21"], labels = c("Low", "Med", "High"))
Figure.2.4 <- p <- ggsurvplot(
  survfit(Surv(PFS, PFS_event) ~ IGLV3.21, data = df_surv),                 
  conf.int = F,         
  size=1,#0.7,                    
  pval=TRUE, 
  palette = my.col[levels(df_surv[, "IGLV3.21"])[levels(df_surv[, "IGLV3.21"])%in%unique(df_surv[, "IGLV3.21"])]],
  pval.method=TRUE,   break.x.by = 500,    
  xlab="Time in days",
  ylab="Progression-free survival probability",
  ylim=c(0,1),
  surv.scale="percent",
  legend.labs = levels(df_surv[, "IGLV3.21"])[levels(df_surv[, "IGLV3.21"])%in%unique(df_surv[, "IGLV3.21"])],
  tables.col="strata",
  risk.table=F, 
  risk.table.col = "strata",
  risk.table.y.text = T,
  legend = c(.7,.8),
  tables.y.text = T, 
  legend.title="", ggtheme = custom_theme(), title = "IGLV3-21")

df_surv[, "IGLV6.57"] <- factor(df_surv[, "IGLV6.57"], labels = c("Low", "Med", "High"))
Figure.2.5 <- p <- ggsurvplot(
  survfit(Surv(OS, OS_event) ~ IGLV6.57, data = df_surv),                 
  conf.int = F,         
  size=1,#0.7,                    
  pval=TRUE, 
  palette = my.col[levels(df_surv[, "IGLV6.57"])[levels(df_surv[, "IGLV6.57"])%in%unique(df_surv[, "IGLV6.57"])]],
  pval.method=TRUE,   break.x.by = 500,    
  xlab="OS_event in days",
  ylab="Overall survival probability",
  ylim=c(0,1),
  surv.scale="percent",
  legend.labs = levels(df_surv[, "IGLV6.57"])[levels(df_surv[, "IGLV6.57"])%in%unique(df_surv[, "IGLV6.57"])],
  tables.col="strata",
  risk.table=F, 
  risk.table.col = "strata",
  risk.table.y.text = T,
  legend = c(.7,.8),
  tables.y.text = T, 
  legend.title="", ggtheme = custom_theme(), title = "IGLV6-57")


df_surv[, "IGKV4.1"] <- factor(df_surv[, "IGKV4.1"], labels = c("Low", "Med", "High"))
Figure.2.6 <- p <- ggsurvplot(
  survfit(Surv(OS, OS_event) ~ IGKV4.1, data = df_surv),                 
  conf.int = F,         
  size=1,#0.7,                    
  pval=TRUE, 
  palette = my.col[levels(df_surv[, "IGKV4.1"])[levels(df_surv[, "IGKV4.1"])%in%unique(df_surv[, "IGKV4.1"])]],
  pval.method=TRUE,   break.x.by = 500,    
  xlab="OS_event in days",
  ylab="Overall survival probability",
  ylim=c(0,1),
  surv.scale="percent",
  legend.labs = levels(df_surv[, "IGKV4.1"])[levels(df_surv[, "IGKV4.1"])%in%unique(df_surv[, "IGKV4.1"])],
  tables.col="strata",
  risk.table=F, 
  risk.table.col = "strata",
  risk.table.y.text = T,
  legend = c(.7,.8),
  tables.y.text = T, 
  legend.title="", ggtheme = custom_theme(), title = "IGKV4-1")


df_surv[, "IGHA1"] <- factor(df_surv[, "IGHA1"], labels = c("Low", "Med", "High"))
Figure.2.7 <- p <- ggsurvplot(
  survfit(Surv(OS, OS_event) ~ IGHA1, data = df_surv),                 
  conf.int = F,         
  size=1,#0.7,                    
  pval=TRUE, 
  palette = my.col[levels(df_surv[, "IGHA1"])[levels(df_surv[, "IGHA1"])%in%unique(df_surv[, "IGHA1"])]],
  pval.method=TRUE,   break.x.by = 500,    
  xlab="OS_event in days",
  ylab="Overall survival probability",
  ylim=c(0,1),
  surv.scale="percent",
  legend.labs = levels(df_surv[, "IGHA1"])[levels(df_surv[, "IGHA1"])%in%unique(df_surv[, "IGHA1"])],
  tables.col="strata",
  risk.table=F, 
  risk.table.col = "strata",
  risk.table.y.text = T,
  legend = c(.7,.8),
  tables.y.text = T, 
  legend.title="", ggtheme = custom_theme(), title = "IGHA1")


df_surv[, "FUT9"] <- factor(df_surv[, "FUT9"], labels = c("Low", "Med", "High"))
Figure.2.8 <- p <- ggsurvplot(
  survfit(Surv(OS, OS_event) ~ FUT9, data = df_surv),                 
  conf.int = F,         
  size=1,#0.7,                    
  pval=TRUE, 
  palette = my.col[levels(df_surv[, "FUT9"])[levels(df_surv[, "FUT9"])%in%unique(df_surv[, "FUT9"])]],
  pval.method=TRUE,   break.x.by = 500,    
  xlab="OS_event in days",
  ylab="Overall survival probability",
  ylim=c(0,1),
  surv.scale="percent",
  legend.labs = levels(df_surv[, "FUT9"])[levels(df_surv[, "FUT9"])%in%unique(df_surv[, "FUT9"])],
  tables.col="strata",
  risk.table=F, 
  risk.table.col = "strata",
  risk.table.y.text = T,
  legend = c(.7,.8),
  tables.y.text = T, 
  legend.title="", ggtheme = custom_theme(), title = "FUT9")


# Results -----------------------------------------------------------------

write.table(Table.S6, "../results/Tables/Table-S6.txt", col.names = T, row.names = F, sep = "\t", quote = F)
ggsave("../results/Figures/Figure-2.1.pdf", plot = print(Figure.2.1), width = 6.5, height = 7, onefile=F)
ggsave("../results/Figures/Figure-2.2.pdf", plot = print(Figure.2.2), width = 6.5, height = 7, onefile=F)
ggsave("../results/Figures/Figure-2.3.pdf", plot = print(Figure.2.3), width = 6.5, height = 7, onefile=F)
ggsave("../results/Figures/Figure-2.4.pdf", plot = print(Figure.2.4), width = 6.5, height = 7, onefile=F)
ggsave("../results/Figures/Figure-2.5.pdf", plot = print(Figure.2.5), width = 6.5, height = 7, onefile=F)
ggsave("../results/Figures/Figure-2.6.pdf", plot = print(Figure.2.6), width = 6.5, height = 7, onefile=F)
ggsave("../results/Figures/Figure-2.7.pdf", plot = print(Figure.2.7), width = 6.5, height = 7, onefile=F)
ggsave("../results/Figures/Figure-2.8.pdf", plot = print(Figure.2.8), width = 6.5, height = 7, onefile=F)



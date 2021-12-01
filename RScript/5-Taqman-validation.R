###################################################
## 
## Taqman - RNASeq validation
##
## - Pearson correlation between RNA-seq and Taqman over 
## 4 genes
##
## Manuscript --> Table S5, Figure S1
##
###################################################

taqman <- read.delim("../data/taqman.txt", dec = ",")
tpm <- read.delim("../data/TPM_IMK_gene.txt")

imk <- colnames(tpm)[-1]

names(imk) <- gsub("Good.IMK", "", gsub("Bad.IMK", "", imk))
taqman$sample%in%names(imk)

taqman$sample <- imk[as.character(taqman$sample)]
colnames(taqman)[2] <- "exp_2_deltaCt"

taqman$deltaCt <- log2(taqman$exp_2_deltaCt)

rownames(taqman) <- paste(taqman$sample, taqman$gene, sep = ":")

tpm <- tpm[tpm$SYMBOL%in%taqman$gene, ]
rownames(tpm) <- tpm$SYMBOL
tpm$SYMBOL <- NULL

taqman$tpm <- NA

res <- list()
res[["corr"]] <- list()
res[["corr.pvalue"]] <- list()


for(i in unique(taqman$gene)){
  dummy <-  taqman[taqman$gene==i, ]
  dummy2 <- unlist(tpm[i, dummy$sample])
  taqman[paste(names(dummy2), i, sep = ":"), "tpm"] <- log2(dummy2+0.0001)
  
  res[["corr"]][[i]] <- cor(taqman[paste(names(dummy2), i, sep = ":"), "tpm"], 
                            taqman[paste(names(dummy2), i, sep = ":"), "deltaCt"])
  res[["corr.pvalue"]][[i]] <- cor.test(taqman[paste(names(dummy2), i, sep = ":"), "tpm"], 
                                        taqman[paste(names(dummy2), i, sep = ":"), "deltaCt"])$p.value
  
}

res$corr
res$corr.pvalue

cor(taqman[taqman$gene!="CDR1", "tpm"], taqman[taqman$gene!="CDR1", "deltaCt"])
cor.test(taqman[taqman$gene!="CDR1", "tpm"], taqman[taqman$gene!="CDR1", "deltaCt"])


taqman$Gene <- taqman$gene

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


ggplot(data=taqman[taqman$Gene!="CDR1", ], aes(x=tpm, y=deltaCt,))+geom_point(aes( col=Gene), size=3)+ 
  geom_smooth(method='lm', alpha=0.2)+
  theme_bw()+
  xlab(expression(log[2]("TPM")))+ylab(expression(paste(Delta, "Ct", sep = "")))+my_theme





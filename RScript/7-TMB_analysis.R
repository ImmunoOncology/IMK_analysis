###################################################
## 
## TMB and somatic mutational signature quantificatio
##
## - TMB
## - Mutational Signature DBS
##
## Manuscript --> Figure S2
##
###################################################

library(ggplot2)
library(dplyr)
library(mutSigExtractor)
library(RColorBrewer)

# TMB ---------------------------------------------------------------------

clinical <- read.delim("../data/clinical.txt")
rownames(clinical) <- paste(clinical$Response, "-", rownames(clinical), sep = "")

tmb <- read.delim("../data/TMB.txt")
tmb$ID <- paste(tmb$Type, "-IMK", gsub("Sample ", "", tmb$Sample), sep = "")

tmb <- tmb[tmb$ID%in%rownames(clinical)[clinical$Diagnosis%in%"cutaneous"], ]


round(mean(tmb$TMB[tmb$Type=="Good"]), 2)
round(sd(tmb$TMB[tmb$Type=="Good"]), 2)

round(mean(tmb$TMB[tmb$Type=="Bad"]), 2)
round(sd(tmb$TMB[tmb$Type=="Bad"]), 2)


Figure.S2A <- ggplot(data=tmb)+geom_bar(aes(x=ID, y=TMB, fill=Type), stat = "identity")+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=16), legend.position = "none"
        , axis.title.x=element_text(size=20)
        , axis.title.y=element_text(size=20)
        , axis.text.y=element_text(size=16))+ scale_fill_manual(values=c("#1174bf", "#eec12b"))+
  annotate("text", x=9.7, y=41, label = expression(paste("TMB mean in good responders = 21.79 (", sigma, "=6.73)")), size=7)+
  annotate("text", x=10, y=36, label = expression(paste("TMB mean in bad responders = 24.82 (", sigma, "=9.5)")), size=7)+xlab("")



# MutationalSignatures ---------------------------------------------------------------------



DBS <- read.delim("../data/MutationalSignatures.DBS78.txt")
rownames(DBS) <- DBS$MutationType
DBS$MutationType <- NULL

DBS_SIGNATURE_PROFILES <- read.delim("../data/COSMIC_v3.2_DBS_GRCh38.txt", row.names = 1)

DBS.res <- list()
DBS.dummy <- rep(NA, nrow(DBS_SIGNATURE_PROFILES))
names(DBS.dummy) <- rownames(DBS_SIGNATURE_PROFILES)

all(rownames(DBS)%in%names(DBS.dummy))

for(i in colnames(DBS)){
  idt <- names(DBS.dummy)
  DBS.dummy <- DBS[idt, i]
  names(DBS.dummy) <- idt
  DBS.res[[i]] <- fitToSignatures(
    mut.context.counts=DBS.dummy, 
    signature.profiles=DBS_SIGNATURE_PROFILES
  )
}

DBS.res.df <- do.call('rbind', DBS.res)
DBS.res.df <- DBS.res.df[, colSums(DBS.res.df)>0]

DBS.res.norm <- as.data.frame(DBS.res.df)#as.data.frame(DBS.res.df/rowSums(DBS.res.df))
#DBS.res.norm$response <- gsub(".*_", "", rownames(DBS.res.norm))
DBS.res.norm$ID <- rownames(DBS.res.norm)
DBS.res.norm <- reshape2::melt(DBS.res.norm)
DBS.res.norm$ID <- unlist(lapply(lapply(strsplit(DBS.res.norm$ID, "_"), rev), paste, collapse = "_", sep = ""))

cols <- brewer.pal(n = 11, name = "Set3")
names(cols) <- levels(DBS.res.norm$variable)

ms <- rev(levels(DBS.res.norm$variable))
DBS.res.norm$variable <- factor(as.character(DBS.res.norm$variable), levels = ms)

DBS.res.norm$ID <- gsub("[_]", "-", DBS.res.norm$ID)

DBS.res.norm$DBS <- DBS.res.norm$variable

Figure.S2B <- ggplot(DBS.res.norm, aes(x=ID, y=value, fill=DBS)) + geom_bar(stat = "identity")+ 
  scale_fill_manual(values=cols)+
  ggplot2::coord_flip()+xlab("")+ylab("")+theme_bw()+theme(
    axis.title.x=element_text(size=20)
    , axis.title.y=element_text(size=20)
    , axis.text.y=element_text(size=16)
    , axis.text.x=element_text(size=16)
    , legend.text=element_text(size=16)
    , legend.title=element_text(size=16, face = "bold")
  )

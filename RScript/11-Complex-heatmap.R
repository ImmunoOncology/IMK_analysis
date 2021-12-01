
library(ComplexHeatmap)
library(circlize)
library(dplyr)

load("../data/complexHeatmapRData.RData")

mat2 = t(apply(mat, 1, function(x) {
  q10 = quantile(x, 0.1)
  q90 = quantile(x, 0.9)
  x[x < q10] = q10
  x[x > q90] = q90
  scale(x)
}))


go.bcell <- unique(unlist(go.bp[go.bp$V1=="GOBP_B_CELL_RECEPTOR_SIGNALING_PATHWAY", ]))
go.bcell <- go.bcell[go.bcell!=""]
go.bcell.l <- sigGenes_SYMBOL %in% go.bcell
go.bcell_gene <-  rownames(mat)[go.bcell.l]

go.ph <- unique(unlist(go.bp[go.bp$V1=="GOBP_PHAGOCYTOSIS_RECOGNITION", ]))
go.ph <- go.ph[go.ph!=""]
go.ph.l <- sigGenes_SYMBOL %in% go.ph
go.ph_gene <-  rownames(mat)[go.ph.l]


go.ab <- unique(unlist(go.mf[go.mf$V1=="GOMF_ANTIGEN_BINDING", ]))
go.ab <- go.ab[go.ab!=""]
go.ab.l <- sigGenes_SYMBOL %in% go.ab
go.ab_gene <-  rownames(mat)[go.ab.l]


go.plasma <- unique(unlist(go.cc[go.cc$V1=="GOCC_EXTERNAL_SIDE_OF_PLASMA_MEMBRANE", ]))
go.plasma <- go.plasma[go.plasma!=""]
go.plasma.l <- sigGenes_SYMBOL %in% go.plasma
go.plasma_gene <-  rownames(mat)[go.plasma.l]

survival_gene <- survival$Gene[survival$OS<0.05]
survival_gene.l <- sigGenes_SYMBOL %in% survival_gene
survival_gene <- rownames(mat)[survival_gene.l]

survival_gene.flag <- survival_gene.l+0
survival_gene.flag[survival_gene.l & sigGenes_SYMBOL%in%resf$SYMBOL[resf$log2FoldChange<0]] <- -1


survival_pfs_gene <- survival_pfs$Gene[survival_pfs$PFS<0.05]
survival_pfs_gene.l <- sigGenes_SYMBOL %in% survival_pfs_gene
survival_pfs_gene <- rownames(mat)[survival_pfs_gene.l]

survival_pfs_gene.flag <- survival_pfs_gene.l+0
survival_pfs_gene.flag[survival_pfs_gene.l & sigGenes_SYMBOL%in%resf$SYMBOL[resf$log2FoldChange<0]] <- -1


tmb$ID <- paste(tmb$Type, "-IMK", gsub("Sample ", "", tmb$Sample), sep = "")
rownames(tmb) <- tmb$ID

colnames(mat2) = colnames(mat)


SBS5 <- mutational_signatures[mutational_signatures$variable %in% "SBS5", ]
SBS5$ID <- gsub("[_]", "-", SBS5$ID)
rownames(SBS5) <- SBS5$ID
SBS5 <- SBS5[colnames(mat2), "value"]

DBS1 <- mutational_signatures[mutational_signatures$variable %in% "DBS1", ]
DBS1$ID <- gsub("[_]", "-", DBS1$ID)
rownames(DBS1) <- DBS1$ID
DBS1 <- DBS1[colnames(mat2), "value"]

library(RColorBrewer)

brewer.pal(n = 8, name = "BuGn")
brewer.pal(n = 8, name = "BuPu")
col_location <- brewer.pal(n = 3, name = "Set3")
names(col_location) <- unique(clinical[colnames(mat2), "Location"])


col_hp = colorRamp2(c(0, 2, 4), c("#E5F5F9", "#66C2A4", "#238B45"))
col_ms_1 = colorRamp2(c(0, max(c(DBS1))), c("#F7FCFD", "#6E016B"))
col_ms_2 = colorRamp2(c(0, max(c(SBS5))), c("#F7FCFD", "#6E016B"))

ha_bottom = HeatmapAnnotation(TMB = anno_barplot(tmb[colnames(mat2), "TMB"])
                              #, "SBS5" = SBS5
                              , "DBS1" = DBS1
                              , "Location" = clinical[colnames(mat2), "Location"]
                              , "BRAF" = clinical[colnames(mat2), "BRAF.mut"]
                              , "Num M1" = clinical[colnames(mat2), "num.M1"]
                              , "No. Previous lines" = clinical[colnames(mat2), "No..Previous.lines"]
                              , "Maximum toxicity grade" = clinical[colnames(mat2), "Maximum.toxicity.grade"]
                              , `Response` = clinical[colnames(mat2), "Responder.type"] 
                              , col = c(
                                "Response" = c("Bad" = "#0073C2FF", "Good" = "#EFC000FF")
                                , "BRAF" = c("No" = "#1B9E77", "Yes"="#D95F02")
                                , "Location" = col_location
                                , "Num M1" = col_hp
                                , "No. Previous lines" = col_hp
                                , "Maximum toxicity grade" = col_hp
                                , "DBS1" = col_ms_1
                                #, "SBS5" = col_ms_2
                              ), annotation_name_side = "left", annotation_legend_param = 
                                list(
                                  "Num M1" = list(title = "Clinical grade")
                                  , "DBS1" = list(title = "Mutational signature score")
                                  , "BRAF" = list(title = "BRAF mutated")
                                )
                              , annotation_label = list(#"SBS5" = "SBS5 (clock-like signature)",
                                "DBS1" = "DBS1 (Ultraviolet light exposure)")
                              , show_legend = c(F, F, F, F, F, F, F))
#, show_legend = c(T, T, T, T, F, F, T))


base_mean_good = rowMeans(mat[, grep("Good", colnames(mat))])
base_mean_bad = rowMeans(mat[, grep("Bad", colnames(mat))])


library(GetoptLong)
library(Seurat)


mks.sig <- mks[mks$p_val_adj<0.05, ]

mks.sig$gene[mks.sig$gene%in%rownames(mat2)]
table(mks.sig$cluster[mks.sig$gene%in%rownames(mat2)])

plasmacells <- mks.sig$gene[mks.sig$gene%in%rownames(mat2) & mks.sig$cluster=="Plasma cells"]
plasmacells.l <- rownames(mat) %in% plasmacells
plasmacells_gene <-  rownames(mat)[plasmacells.l]

plasmablast <- mks.sig$gene[mks.sig$gene%in%rownames(mat2) & mks.sig$cluster=="Plasmablast"]
plasmablast.l <- rownames(mat) %in% plasmablast
plasmablast_gene <-  rownames(mat)[plasmablast.l]

bigk <- mks.sig$gene[mks.sig$gene%in%rownames(mat2) & mks.sig$cluster=="Naive B cells IGK-high"]
bigk.l <- rownames(mat) %in% bigk
bigk_gene <-  rownames(mat)[bigk.l]

bigl <- mks.sig$gene[mks.sig$gene%in%rownames(mat2) & mks.sig$cluster=="Naive B cells IGL-high"]
bigl.l <- rownames(mat) %in% bigl
bigl_gene <-  rownames(mat)[bigl.l]

naive <- mks.sig$gene[mks.sig$gene%in%rownames(mat2) & mks.sig$cluster=="Naive B cells"]
naive.l <- rownames(mat) %in% naive
naive_gene <-  rownames(mat)[naive.l]

survival_gene.flag.label <- ifelse(survival_gene.flag==0, "No significant", survival_gene.flag)
survival_gene.flag.label <- ifelse(survival_gene.flag.label=="1", "Associated with better survival", survival_gene.flag.label)
survival_gene.flag.label <- ifelse(survival_gene.flag.label=="-1", "Associated with worser survival", survival_gene.flag.label)

survival_pfs_gene.flag.label <- ifelse(survival_pfs_gene.flag==0, "No significant", survival_pfs_gene.flag)
survival_pfs_gene.flag.label <- ifelse(survival_pfs_gene.flag.label=="1", "Associated with better survival", survival_pfs_gene.flag.label)
survival_pfs_gene.flag.label <- ifelse(survival_pfs_gene.flag.label=="-1", "Associated with worser survival", survival_pfs_gene.flag.label)


ht1 <- Heatmap(mat2, col = colorRamp2(c(-1.5, 0, 1.5), c("royalblue3", "white", "red3")), 
               name = "Gene scaled expression", #column_title = qq("relative expression for @{nrow(mat)} genes"),
               show_column_names = FALSE, width = unit(8, "cm"),
               heatmap_legend_param = list(title = "Gene scaled expression")
               , right_annotation = rowAnnotation("B cell receptor signaling pathway"=c(go.bcell.l + 0)
                                                  , "Phagocytosis, recognition" = c(go.ph.l + 0)
                                                  , "Antigen binding" = c(go.ab.l + 0)
                                                  , "PFS" = c(survival_pfs_gene.flag.label)
                                                  , "OS" = c(survival_gene.flag.label)
                                                  , "Plasma cells" = c(plasmacells.l + 0)
                                                  , "Plasmablast" = c(plasmablast.l + 0)
                                                  , "Naive B cells IGK-high" = c(bigk.l + 0)
                                                  , "Naive B cells IGL-high" = c(bigl.l + 0)
                                                  , "Naive B cells" = c(naive.l + 0)
                                                  , col = list(
                                                    "B cell receptor signaling pathway" =c("0" = "white", "1" = "red2")
                                                    , "Phagocytosis, recognition" =c("0" = "white", "1" = "red2")
                                                    , "Antigen binding" =c("0" = "white", "1" = "red2")
                                                    , "PFS" =c("No significant" = "white", "Associated with better survival" = "#1B9E77", "Associated with worser survival"="azure4")
                                                    , "OS" =c("No significant" = "white", "Associated with better survival" = "#1B9E77", "Associated with worser survival"="azure4")
                                                    , "Plasma cells" =c("0" = "white", "1" = "gold")
                                                    , "Plasmablast" =c("0" = "white", "1" = "gold")
                                                    , "Naive B cells IGK-high" =c("0" = "white", "1" = "gold")
                                                    , "Naive B cells IGL-high" =c("0" = "white", "1" = "gold")
                                                    , "Naive B cells" =c("0" = "white", "1" = "gold")
                                                  )
                                                  , annotation_legend_param = list(OS = list(title = "Survival (PFS and OS)"))
                                                  #, show_legend = c(F, F, F, F, T, F, F, F, F, F), annotation_name_side = "top")
                                                  , show_legend = c(F, F, F, F, F, F, F, F, F, F), annotation_name_side = "top")
               , top_annotation=ha_bottom
               #, bottom_annotation = ha_bottom
               , show_row_names = F
               , show_heatmap_legend = F
)

z.mat <- t(scale(hla, center=TRUE, scale=TRUE))
z.mat2 = t(apply(z.mat, 1, function(x) {
  q10 = quantile(x, 0.1)
  q90 = quantile(x, 0.9)
  x[x < q10] = q10
  x[x > q90] = q90
  scale(x)
}))

colnames(z.mat2) <- colnames(z.mat)

ht2 = Heatmap(z.mat2, name = "HLA scaled expression", col = colorRamp2(c(-1.5, 0, 1.5), c("yellow", "white", "green")), cluster_columns = F, show_row_names = F, show_heatmap_legend = F)


rownames(vdj_vol) <- paste(vdj_vol$Response, vdj_vol$Sample, sep = "-")
vdj_vol <- vdj_vol[colnames(z.mat2), ]
imm_top <- imm_top[colnames(z.mat2), ]
imm_top[1, ] <- 0

dummy <- imm_top

for(i in 2:ncol(imm_top)){
  if(i>2)
    dummy[, i] = dummy[, i] - rowSums(dummy[, seq(i-1, 1)])
  else
    dummy[, i] = dummy[, i] - (dummy[, (i-1)])
}

dummy <- dummy[, ncol(dummy):1]

net_mixcr$ID <- paste(net_mixcr$Condition, net_mixcr$Sample, sep = "-")
net_mixcr <- net_mixcr[net_mixcr$ID %in% colnames(z.mat), ]
net_mixcr$VDJ <- substr(net_mixcr$bestVGene, 1, 3)
net_mixcr <- aggregate(net_mixcr$cloneCount, by=list(ID=net_mixcr$ID, VDJ=net_mixcr$VDJ), FUN=sum)
net_mixcr <- reshape2::dcast(formula = ID~VDJ, net_mixcr, value.var = "x")
rownames(net_mixcr) <- net_mixcr$ID
net_mixcr$ID <- NULL
net_mixcr <- net_mixcr[colnames(z.mat2), ]
rownames(net_mixcr) <- colnames(z.mat2)
net_mixcr[is.na(net_mixcr)] <- 0

vdj = HeatmapAnnotation(VDJ = anno_barplot(vdj_vol[, c("TCR", "BCR")], gp=gpar(fill=c("green" ,"orange"))), annotation_name_side = "left")
vdj = HeatmapAnnotation(VDJ = anno_barplot(net_mixcr, gp=gpar(fill=2:7)), annotation_name_side = "left")
breaks = c(1, 5, 50, 100, 500, 1000, 10000)
vdj = HeatmapAnnotation(VDJ = anno_barplot(log10(net_mixcr+1), gp=gpar(fill=2:7), axis_param = list(at = log10(breaks), labels =breaks)), annotation_name_side = "left")

vdj.mat <- t(scale(log10(net_mixcr+1), center=TRUE, scale=TRUE))
vdj.mat2 = t(apply(vdj.mat, 1, function(x) {
  q10 = quantile(x, 0.1)
  q90 = quantile(x, 0.9)
  x[x < q10] = q10
  x[x > q90] = q90
  scale(x)
}))


brewer.pal(n = 8, name = "Accent")

vdj = Heatmap(vdj.mat2,  name = "VDJ scaled expression", col = colorRamp2(c(-1.5, 0, 1.5), c("purple", "white", "orange")), cluster_columns = F, show_row_names = F, show_heatmap_legend = F)
vdj2 = HeatmapAnnotation(VDJ = anno_barplot(dummy,
                                            gp = gpar(fill=c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0")))
                         , annotation_name_side = "left", show_legend=F)

lgd = Legend(labels = c(
  "Clonotype indexes [1:5)", 
  "Clonotype indexes [5:20)", 
  "Clonotype indexes [20:100)", 
  "Clonotype indexes [100:500)", 
  "Clonotype indexes [500:5000)"
), legend_gp = gpar(fill = c("#386CB0", "#FFFF99", "#FDC086", "#BEAED4", "#7FC97F"), col = c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0")), 
title = "VDJ clonal expansion", ncol = 1)

lgd2 = Legend(labels = c(
  "Genes associated with pathways", 
  "Genes significant for prognosis survival", 
  "Top markers in scRNA-seq"
), legend_gp = gpar(fill = c("red2", "#1B9E77", "gold"), col = c("red2", "#1B9E77", "gold")), 
title = "Gene annotation", ncol = 1)

ht1 %v% vdj2 %v% vdj %v% ht2
draw(lgd, x = unit(0.6, "npc"), y = unit(0.2, "npc"))
draw(lgd2, x = unit(0.8, "npc"), y = unit(0.2, "npc"))


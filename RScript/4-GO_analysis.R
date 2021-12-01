###################################################
## 
## Gene Ontology enrichment
##
## - BP enrichGO
## - CC enrichGO
## - MF enrichGO
##
## Manuscript --> Table S4
##
###################################################


library(clusterProfiler)
library(org.Hs.eg.db)

genemap <- read.table("../data/Result_RNASeq_DESeq2_melanoma_cutaneous.txt", stringsAsFactors = F, header = T)
rownames(genemap) <- genemap$ENSEMBLENSEMBL

up <- genemap[genemap$log2FoldChange>0, "ENTREZ"]
down <- genemap[genemap$log2FoldChange<0, "ENTREZ"]

result_go <- list()

for(i in c("BP", "CC", "MF")){
  
  ego.up <- enrichGO(gene       = up,
                     OrgDb         = org.Hs.eg.db,
                     ont           = i,
                     pAdjustMethod = "BH",
                     readable      = TRUE)
  
  ego.up.raw <- ego.up
  ego.up <- clusterProfiler::simplify(ego.up, by="p.adjust")
  
  ego.down <- enrichGO(gene       = down, 
                       OrgDb         = org.Hs.eg.db,
                       ont           = i,
                       pAdjustMethod = "BH",
                       readable      = TRUE)
  
  ego.down.raw <- ego.down
  ego.down <- clusterProfiler::simplify(ego.down, by="p.adjust")
  
  result_go[[i]] <- list(up=ego.up, up.raw=ego.up.raw, down=ego.down
                                   , down.raw = ego.down.raw)
}

result_go$BP$up@result$GO_type <- "BP"
result_go$CC$up@result$GO_type <- "CC"
result_go$MF$up@result$GO_type <- "MF"

go_up <- rbind(result_go$BP$up@result, result_go$CC$up@result, result_go$MF$up@result)
go_down <- rbind(result_go$BP$down@result, result_go$CC$down@result, result_go$MF$down@result)

Table.S4 <- go_up[, c("Description", "Count", "p.adjust", "GO_type")]

# Results -----------------------------------------------------------------

write.table(Table.S4, "../results/Tables/Table-S4.txt", col.names = T, row.names = F, sep = "\t", quote = F)


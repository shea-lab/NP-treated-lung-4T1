#' Run Cluster Profiler and generate plots for output 
#' 
#' @param sample seurat object with idents set 
#' @param organism organism for genes: human, mouse
#' 
#' @return results_list list of clusterProfiler GO, KEGG, DO, and Pathway clusters



prep_clusterProfiler <- function(sample, organism) {
  # load libraries 
  library(Seurat)
  library(clusterProfiler)
  library(dplyr)
  
  # #' @importFrom Seurat Read10X
  
  # 
  # Idents(subgroup) <- subgroup$phenotype
  
  if (organism == "human") {
    library(org.Hs.eg.db)
    org_db = "org.Hs.eg.db"
    org = 'hsp'
    org_long = 'human'
  } else if (organism == "mouse") {
    library(org.Mm.eg.db)
    org_db = "org.Mm.eg.db"
    org = "mmu"
    org_long = 'mouse'
  } else {
    return(print("organism not supported"))
  }

  
  markers <- FindAllMarkers(subgroup, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
  
  markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
  
  ##################################################################
  #Subsetting top 100 markers with adjusted p values lower than .05#
  ##################################################################
  top100 <- markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
  top100pval <- subset(top100, rowSums(top100[5] < 0.05) > 0)
  
  df <- top100pval[,7:6]
  dfsample <- split(df$gene,df$cluster)
  length(dfsample)
  
  dfsample <- split(df$gene,df$cluster)
  for (i in 1:length(dfsample)) {
    genes = dfsample[[i]]
    dfsample[[i]] = bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb= org_db)
  }
  
  # dfsample$`Mature` = bitr(dfsample$`Mature`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  # dfsample$`Inflammatory` = bitr(dfsample$`Inflammatory`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  # dfsample$`Unknown 1` = bitr(dfsample$`Unknown 1`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  # dfsample$`Transitional` = bitr(dfsample$`Transitional`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  # dfsample$`Unknown 2` = bitr(dfsample$`Unknown 2`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  # dfsample$`Immature` = bitr(dfsample$`Immature`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  # dfsample$`Nonclassical Monocytes` = bitr(dfsample$`Nonclassical Monocytes`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  # dfsample$`Proliferating` = bitr(dfsample$`Proliferating`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  # dfsample$`8` = bitr(dfsample$`8`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  # dfsample$`9` = bitr(dfsample$`9`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  # dfsample$`10` = bitr(dfsample$`10`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  # dfsample$`11` = bitr(dfsample$`11`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  # dfsample$`12` = bitr(dfsample$`12`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  # dfsample$`13` = bitr(dfsample$`13`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  # dfsample$`14` = bitr(dfsample$`14`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  # dfsample$`8` = bitr(dfsample$`8`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  
  genelist = dfsample
  
  for (i in 1:length(dfsample)) {
    name = names(dfsample)[[i]]
    genelist[[name]] <- dfsample[[i]]$ENTREZID
  }
  
  
  #do the same here, a line like below for each cluster
  # genelist <- list("Mature" = dfsample$`Mature`$ENTREZID, 
  #                  "Inflammatory" = dfsample$`Inflammatory`$ENTREZID,
  #                  "Unknown 1" = dfsample$`Unknown 1`$ENTREZID,
  #                  "Transitional" = dfsample$`Transitional`$ENTREZID,
  #                  "Unknown 2" = dfsample$`Unknown 2`$ENTREZID,
  #                  "Immature" = dfsample$`Immature`$ENTREZID,
  #                  "Nonclassical Monocytes" = dfsample$`Nonclassical Monocytes`$ENTREZID,
  #                  "Proliferating" = dfsample$`Proliferating`$ENTREZID
                   # "8" = dfsample$`8`$ENTREZID, 
                   # "9" = dfsample$`9`$ENTREZID,
                   # "10" = dfsample$`10`$ENTREZID,
                   # "11" = dfsample$`11`$ENTREZID,
                   # "12" = dfsample$`12`$ENTREZID,
                   # "13" = dfsample$`13`$ENTREZID,
                   # "14" = dfsample$`14`$ENTREZID
                   
  return(genelist)

  # GOclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = org_db)
  # dotplot(GOclusterplot)
  # 
  # results <- GOclusterplot@compareClusterResult
  # 
  # KEGGclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichKEGG", org = org)
  # dotplot(KEGGclusterplot)
  # 
  # Pathwayclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichPathway", organism = org_long)
  # dotplot(Pathwayclusterplot)
  # 
  # DOclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichDO")
  # dotplot(DOclusterplot)
  # 
  # results_table <- GOclusterplot$
  #   
  # results_list <- list(GOclusterplot, KEGGclusterplot, DOclusterplot)
  
  # return(GOclusterplot)
  # return(KEGGclusterplot)
  # return(Pathwayclusterplot)
  # return(DOclusterplot)
  # return(results_list)

}

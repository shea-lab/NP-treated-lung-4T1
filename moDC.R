# Setup -------------------------------------------------------------------
library(hypeR)
library(Seurat)
library(dplyr)
library(readr)
library(tidyverse)

# Load libraries
library(SingleCellExperiment)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(clusterProfiler)
library(RColorBrewer)
library(msigdbr)
library(readxl)
library(ComplexHeatmap)

library(DoMultiBarHeatmap)
library(dittoSeq)
library(GO.db)
library(patchwork)

convert_mouse_to_human <- function(gene_list) { 
  output = c()
  mouse_human_genes = read.csv("https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
  
  for(gene in gene_list) {
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name == "mouse, laboratory"))[['DB.Class.Key']]
    if( !identical(class_key, integer(0)) ) {
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
      for(human_gene in human_genes) {
        output = rbind(c(gene, human_gene), output)
      }
    }
  }
  return (output)
}


# Load Subgroup -----------------------------------------------------------


# not needed if you load stromal-free macrophage lung_all_macrophage_no_stromal.rds
Idents(data) <- "phenotype"
subgroup <- subset(data, subset = phenotype == "moDC")


# setup files 
subgroup.list <- c("microglia.rds",  "monocyte.rds",  
                   "clean_integrated_neutrophil.rds",   "integrated_stromal.rds",    "clean_integrated_vascular.rds" )
cell.list <- c("Microglia", "Monocyte", "Neutrophil", "Stromal", "Vascular")
subgroup.files <- c("Microglia_deg.csv",  "monocyte_deg.csv",   "neutrophil_deg.csv", "Stromal_deg.csv" ,   "vascular_deg.csv")

results.path = "~/Documents/Kate/NMMscRNAseq/Results"
data.path = "~/Documents/Kate/NMMscRNAseq/Data/"


subgroup.list <- c("Untreated", "")


# Generate DEGs -----------------------------------------------------------

for (i in 1:length(subgroup.list)) {
  # load Seurat object
  subgroup <- readRDS(paste0(data.path, subgroup.list[i]))
  cell = cell.list[i]
  
  # set comparisons 
  Idents(subgroup) <- "sample"
  
  # prepare for DEG comparison
  subgroup <- PrepSCTFindMarkers(subgroup)
  # DefaultAssay(object = subgroup) <- "SCT"
  # 
  # subgroup[["RNA"]] <- JoinLayers(subgroup[["RNA"]])
  # subgroup <- JoinLayers(subgroup)
  
  # Run DEG Comparison
  homeostatic.de <- FindMarkers(subgroup, ident.1 = "d2.pbs", ident.2 = "d2.np", min.pct=0.25, logfc.threshold=0.25) %>%    
    mutate(comp = "d2.pbs_v_d2.np")   %>% rownames_to_column(var = "gene")
  
  new <-FindMarkers(subgroup, ident.1 = "d7.pbs", ident.2 = "d7.np", min.pct=0.25, logfc.threshold=0.25) %>% 
    mutate(comp = "d7.pbs_v_d7.np") %>% rownames_to_column(var = "gene")
  homeostatic.de <- bind_rows(homeostatic.de, new)
  
  # new <-FindMarkers(subgroup, ident.1 = "d2.pbs", ident.2 = "d7.pbs", min.pct=0.25, logfc.threshold=0.25) %>% mutate(comp = "d2.pbs_v_d7.pbs") %>% rownames_to_column(var = "gene")
  # homeostatic.de <- bind_rows(homeostatic.de, new)
  # 
  # new <-FindMarkers(subgroup, ident.1 = "d2.np", ident.2 = "d7.np", min.pct=0.25, logfc.threshold=0.25) %>% mutate(comp = "d2.np_v_d7.np") %>% rownames_to_column(var = "gene")
  # homeostatic.de <- bind_rows(homeostatic.de, new)
  
  write.csv(homeostatic.de, file=paste0(results.path,cell,"_condition_deg.csv"))
  rm(new, subgroup)
}


data2$samplename <- data2$shortname

Idents(data2) <- "shortname"
data2 <- RenameIdents(data2, 'lung_np' = "NP+", 'lung_imm' = "NP-", 'lung_tb' = "Untreated")
data2@meta.data$shortname <- data2@active.ident 
# subgroup@meta.data$identity <- NULL

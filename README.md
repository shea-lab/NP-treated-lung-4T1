# NP-treated-lung-4T1
Repository for "Induction of antigen-presenting monocyte-derived dendritic cells by nanoparticles inhibits metastasis and relieves immunosuppression in the metastatic niche"

## Files 
CellChat_Comparison_Analysis.Rmd - compare NP treated and untreated groups with cellchat 

CellChat_Figures.Rmd - generate cellchat figures 
Cellchat_Object_Prep.Rmd - convert seurat object to cellchat objects
DEG Analysis and cluster Profiler.Rmd - 
Flex_Data_Preprocessing.Rmd - run QC, cleanup data, combine samples to aggregate object 
Generate_Figures.Rmd - generate plots 
Label_Clusters.Rmd - label clusters (incl lung and PT)
Lung_Label_Clusters.Rmd - label clusters in the lung, cell identification
Lung_Label_Clusters_NoStro.Rmd - label clusters in the lung, no stromal sample (CD45-)
moDC_Th17_Mechanism.Rmd - generate moDC and Th17 DEGs and plots 
Myeloid_Analysis.Rmd - analyze myeloid cells in lung 
Seurat_Processing.Rmd - create and pre-process Seurat object (filter and cluster, pre-cluster ID) 
v2_CellChat_Comparison_Analysis.Rmd - compare updated cellchat to old cellchat 
moDC.R - examine moDC population 

## Functions
find_doublets.R - run findDoublets on per sample basis 
remove_soup.R - run the SoupX package on per sample basis
run_ddqcR.R - run ddqcR clean up on per sample basis 
prep_clusterProfiler.R - take seurat object, find DEGs based on set ident, return object for clusterProfiler 
run_clusterProfiler.R - defunct 

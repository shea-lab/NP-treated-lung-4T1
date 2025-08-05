# NP-treated-lung-4T1
Repository for "Induction of antigen-presenting monocyte-derived dendritic cells by nanoparticles inhibits metastasis and relieves immunosuppression in the metastatic niche"

## Data 
Raw sequencing data is available at GEO under accession number GSE302362. Processed Seurat object is available [here](https://www.dropbox.com/scl/fi/jxbu7i7u78pkxthvavepl/labeled_lung_int_merge_sct_filt_tnbc.rds?rlkey=y2vn60d8g0ytfr7bwmbxw9abf&st=trqx3omm&dl=0). Lungs were isolated from female BALB/c 4T1-tumor bearing mice at D14 post-inoculation (inoculated with 50 uL of 2e6 4T1 cells/mL). Mice were treated with saline (control) or 1mg of NPs every three days starting at D1 post-inocluation (D1, 4, 7, 10, 13 NP treatments). Cells were then prepared for flow-assisted cell sorting by staining with AF488-CD45 and sorted into CD45+NP+, CD45+NP-, and CD45- populations on the Bigfoot Cell Sorter. 
<p> Sorted samples were fixed with the 10x Chromium Next GEM Single Cell Fixed RNA Sample Preparation Kit following manufacturer instructions. Fixed samples were submitted to the University of Michigan Advanced Genomics Core for library preparation with the 10x Single Cell Fixed RNA Hybridization & Library Kit and Single Cell Fixed RNA Mouse Transcriptome Probe Kit following manufacturer guidelines for the 10x Single Cell Gene Expression Flex platform. Samples were sequenced on the NovaSeq X 10B flow cells (300 cycle) at a depth of approximately 25,000 reads per cell. CellRanger v8.0.0 was used to align reads to mouse reference genome GRCm38 (mm10-2020-A). </p>


## Files 
In order of processing: 
+ Flex_Data_Preprocessing.Rmd - run QC, cleanup data, combine samples to aggregate object
+ Seurat_Processing.Rmd - create and pre-process Seurat object (filter and cluster, pre-cluster ID)
+ Lung_Label_Clusters.Rmd - label clusters in the lung, cell identification
+ Generate_Figures.Rmd - generate plots
+ Cellchat_Object_Prep.Rmd - convert seurat object to cellchat objects
+ CellChat_Comparison_Analysis.Rmd - compare NP treated and untreated groups with cellchat
+ v2_CellChat_Comparison_Analysis.Rmd - v2 script uses updated cellchat objects without dropped unknowns (unknowns were ID'd) 4/7/25
+ CellChat_Figures.Rmd - generate cellchat figures 
+ moDC_Th17_Mechanism.Rmd - generate moDC and Th17 DEGs and plots 
+ Myeloid_Analysis.Rmd - contains more granular markers for myeloid subtype identification
+ moDC.R - examine moDC DEGs, generate .csv  

## Functions
+ find_doublets.R - run findDoublets on per sample basis 
+ remove_soup.R - run the SoupX package on per sample basis
+ run_ddqcR.R - run ddqcR clean up on per sample basis 
+ prep_clusterProfiler.R - take seurat object, find DEGs based on set ident, return object for clusterProfiler 

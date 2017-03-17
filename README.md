# TCGA-OV-miRNA-DHS-motif
ECE 286 Final Project in collaboration with Kevin Blansit.


These scripts are dependent on the following installations:
  R (3.3.3)  
  
  R libraries:  
  'dplyr'  
  'data.table'  
  'tidyverse'  
  'dtplyr'  
  'survival'  
  'randomForestSRC'  
  'TCGAbiolinks'
  
  MEME Suite (4.11.2)  
  HOMER (v4.9) with hg19 dataset 
  
# Execution
Run the TCGA-OV-miRNA-DHS-motif.ipynb jupyter notebook until the step "Get gene names of important genes"  

Copy the following files into "kevin/data":  
DHS_subset.miRNA_quantifications.reads_per_million_miRNA_mapped.transposed.txt  
clinical.patient.drug.txt  
clinical.patient.txt  
miRNA_quantifications.reads_per_million_miRNA_mapped.transposed.txt  

At this point, enter the "kevin" directory, then run the R script for random forest (data_pipeline.R)  

Copy all files from "kevin/out" into the main directory, then continue running the notebook.

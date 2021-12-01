## IMK analysis

Code for the manuscript "A gene signature of high IGKC expressing intratumoral Plasma cells predict response to Immune Checkpoint Blockade". The scripts were separated in several sections:

  1. Clinico-pathological-analysis.R: descriptive, bivariate and survival analysis for the clinico pathological variables of the melanoma cohort and cutaneous melanoma cohort. 
  2. DESeq2_melanoma.R: differential expression analysis with DESeq2 for the melanoma cohort. 
  3. DESeq2_melanoma_cutaneous.R: differential expression analysis with DESeq2 for the cutaneous melanoma cohort. 
  4. GO_analysis.R: Gene Ontology enrichment analysis with the DE genes from the cutaneous melanoma cohort. 
  5. Taqman-validation.R: Validation of 4 DE genes by pearson correlation from TPM values of RNA-seq and Taqman.
  6. Survival.R: Survival analysis for the expression value of the 140 DE genes from the cutaneous melanoma cohort. 
  7. TMB_analysis: TMB and mutational signature for the cutaneous melanoma cohort. 
  8. sc-RNAseq-analysis.R: Single cell RNA-seq analysis of Sade-Feldman M et al. data. 
  9. Deconvolution.R: in silico cytometry method to obtained quantification of cell types in our RNA-seq.
  10. TCR-BCR-HLA.R: TCR, BCR and HLA profiling of the cutaneous melanoma cohort. 
  11. Machine-Learning.R: Machine Learning analysis to establish the predictive power of our gene signature.



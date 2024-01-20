# HB_singlecell

Below is a list of files in this package.

Data files: 
	1. Lung map (CountMatrix.Rdata): Single Nucleus Multiomic Profiling Reveals Age-Dynamic Regulation of Host Genes Associated with SARS-CoV-2 Infection (GSE161382, Wang et al. );
	2. IPF scRNA-seq (GSE135893_meta.csv): Single-cell RNA sequencing reveals profibrotic roles of distinct epithelial and mesenchymal lineages in pulmonary fibrosis (GSE135893, Habermann et al. );
	3. IPF bulk RNA-seq (Bulk_IPF_GSE150910.csv): Chronic hypersensitivity pneumonitis, an interstitial lung disease with distinct molecular signatures (GSE150910, Furusawa et al. );
	4. Prob_mean.csv: Probability for a cell being an alveolar macrophage cell obtained by HierXGB;
	5. ind_data_1.csv: pseudo bulk RNA data calculated based on IPF scRNA-seq;
	6. HypoMap (Neuro_normal_chow.RData): HypoMap—a unified single-cell gene expression atlas of the murine hypothalamus (University of Cambridge’s Apollo Repository (doi:10.17863/CAM.87955), Steuernagel et al.);

Code for IPF data analysis:
	1. info_parallel_modified.R: generate results for genes with informative prior;
	2. non_info_parallel.R: generate results for genes with  non-informative prior;
	3. d230203_combine_results_across_genes.Rmd: combine results from info_parallel_modified.R and non_info_parallel.R;

Code for simulation: 
	1. DE_hybrid_info_2.5.R: hybrid with informative prior mean of 2.5;
	2. DE_hybrid_info.R: hybrid with informative prior mean of 5;
	3. DE_hybrid_noninfo.R: hybrid with non-informative prior;
	4. DE_MAST.R: code for running MAST;
	5. DE_nebula.R: code for running nebula;
	6. DE_pseudo_mean_DESeq2.R: code for running DESeq2 on pseudo bulk data;
	7. DE_pseudo_mean_edgeR.R: code for running edgeR on pseudo bulk data;
	8. DE_pseudo_mean_limma.R: code for running limma on pseudo bulk data;
	9. main_fn.R: hybrid method code implemented in DE_hybrid_info_2.5.R, DE_hybrid_info.R, and DE_hybrid_noninfo.R.


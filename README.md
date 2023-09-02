# (Temporary) Dual-axis RNA subtyping of pancreatic ductal adenocarcinoma reflecting tumor and its microenvironment using SNUBH cohort with  unbiased real-world stage
A New Perspective on Pancreatic Cancer Subtyping and Biological Insights from RNA Sequencing

# Method
![Method2](https://github.com/JoonghoLee/PDAC_Subtyping/assets/35910715/b91963d1-e26c-403e-ae90-8aa76d0275dd)
Supplemental Figure 1. scRNA-seq. annotation results from each individual data. A: UMAP of annotation results. The dots represent individual single cells. The colors indicate, in order, the major cell-type, whether the cell is a tumor cell, and the patient samples. B: Cell-type proportions of PDAC samples. Endothelial cells, acinar cells, ductal cells, fibroblasts, stellate cells, myeloid cells, T cells, B cells, and mast cells were detected in both datasets. In scRNA-seq. data from GSE154778 and GSE156405, epithelial cells (ductal cells) were commonly detected across multiple patients, while cells such as B cells or acinar cells were not commonly detected in all patients. Specifically, B cells could only be detected in two patients in the GSE154778 dataset, and the deconvolution result using bulk RNA sequencing data in this study also confirmed B cells in very few patients.
![Method](https://github.com/JoonghoLee/PDAC_Subtyping/assets/35910715/839e4560-8627-46ea-b7b7-78d860cf9bdf)
Supplemental Figure 2. Schematic overview of the workflow to identify the Tumor/TME subtype.
A: Extraction process of pre-selected gene set for Tumor/TME subtyping. After performing cell-type annotation and tumor cell annotation, we extracted pre-selected genes for Tumor/TME subtyping. B: 2-stage ledein clustering process using a pre-selected gene set. We performed 2-stage ledein clustering using the pre-selected gene set and finally obtained a total of 110 genes for Tumor subtyping and 51 genes for TME subtyping.


## Subtyping procedure
In this work, two different subtypes are considered, i.e., tumor subtype and TME subtype. The former is to stratify PDAC according to characteristics of tumor cell itself, whereas the latter is to count immune reaction to the tumor cell invasion. The subtyping should be accompanied by their respective gene sets, based on which the transcriptomic profile of PDAC patients (and, hopefully, their response to anti-cancer therapies) can be well stratified. Our approach to find a good gene set (even if not the best) is briefly depict in Supplemental Figure 2. It can be divided in two steps:
1	Pre-selection of gene set (for tumor and TME subtyping)
2	Gene set refinement and subtyping via consensus clustering
In the first step, a large set of genes is selected first, separately for tumor subtype and TME subtype. Then, in the second step, the gene sets are refined by selecting those ones that stratify well the entire cohort into two groups of patients. For the latter step, three data sets were used to avoid any possible data dependencies. The processing details of the two steps are as follows.
## Pre-selection of gene set
For tumor subtype, we used public single-cell RNA-seq data to extract tumor-cells and perform DEG to finally collect a pool of genes to be used in the second step.
- Single-cell RNA-seq data for pre-selection of gene set for tumor subtype: For tumor subtyping, we need tumor markers of PDAC cancer cells. Two single-cell RNA sequencing (scRNA-seq) datasets, GSE154778 (Lin, Noel et al. 2020) and GSE156405 (Lee, Bernard et al. 2021), of PDAC patients were downloaded and used to select genes that are overexpressed in cancer cells. The genes were then used for subtyping of PDAC. GSE154778 consisted of 10 primary tumor samples and 6 metastatic tumor samples, of which we utilized 10 primary tumor samples for our analysis. GSE156405 included 5 primary tumor samples, 4 metastatic tumor samples, 8 peripheral blood mononuclear cell (PBMC) samples, and 2 organoid samples, and we used 5 primary tumor samples for our analysis.
- Preprocessing and filtering of single-cell RNA-seq data: SCANPY (Python package) (Wolf, Angerer et al. 2018) was used for data preprocessing; cells were selected based on three conditions: (1) total number of genes expressed is more than 400 and less than 4,000; (2) total number of unique molecular identifiers (UMI) is less than 15,000; and (3) mitochondrial gene percentage is less than 20%. Subsequently, counts per million (CPM) normalization and log-transformation were applied before further analysis.
- Cell-type identification: We utilized HiCAT (Lee, Kim et al. 2023), marker-based cell-type annotation tool, to annotate major cell type. The marker genes used in the cell-type annotation are shown in Supplemental Table 1. The gene set analysis score was calculated for major-types (T cell, B cell, Myeloid cell, Fibroblast, Mast cell, Ductal cell, Endothelial cell, Granulocyte, Acinar cell) and the highest scored cell-type was assigned to each cluster. For tumor cell annotation, python implementation of InferCNV (https://icbi-lab.github.io/infercnvpy) was used to infer copy number alteration (CNA). We used normal cells including immune cells and stromal cells as reference to estimate CNA of tumor cells. And then highest scored CNA cluster was assigned to tumor cell. In the case of GSE154779 dataset, 1410 tumor cells were identified, and in the case of GSE156405 dataset, 2560 tumor cells were identified. Annotation results and cell-type proportions are shown in Supplemental Figure 1.
- Pre-selection procedure for tumor subtype: To find PDAC tumor marker genes, we first selected primary cancer samples from the two datasets, then performed cell type annotation separately for the two by using HiCAT, for which we used cell markers for pancreas tissue used in MarkerCount(Kim, Lee et al. 2022)  (Supplemental Table 1). The cell marker database contains markers of immune cells and those cell types of pancreas tissue, including ductal cells, whereas it does not provide tumor cell markers to identify tumor cells. To find tumor cells from the single-cell RNA-seq datasets, we used the InferCNV (Patel, Tirosh et al. 2014, Tirosh, Izar et al. 2016) by setting those cells that are annotated as immune cells (T cells, B cells and myeloid cells) to normal reference. With the CNV scores obtained, we set the cells from the clusters having high CNV scores as tumor cells, which were mostly the ductal cells, according to HiCAT annotations. To find tumor markers, SCANPY was then used to perform differentially expressed genes (DEG) analysis by comparing tumor cells with others, excluding the ductal cells that were not chosen as tumor cells. We filtered out those DEGs with adjusted p-value greater than 0.05 and, finally, the tumor markers were identified by taking intersection of filtered DEGs from the two datasets. Through this process, we obtained 3,702 genes commonly overexpressed (adjusted p-value<0.05) in tumor cells from two datasets.
- Pre-selection of gene set for TME subtype: For TME subtype identification, we simply used marker genes of immune and stromal cells, including, T cell, B cell, myeloid cell, Fibroblast and smooth muscle cells. The marker genes of these cells were collected from R&D systems (https://www.rndsystems.com/resources/cell-markers), where various subset markers of immune and stromal cells are available. Specifically, we used 332 marker genes. (Supplemental Figure 2A).

## Gene set refinement and subtyping via consensus clustering
In the second step, we start from the pre-selected genes to refine for PDAC subtype identification using , and performed validation for them.
- Bulk RNA-seq data used: We used three datasets, including our collection (SNUBH), TCGA (Nawy 2018) and Bailey (Bailey, Chang et al. 2016) dataset. TCGA transcriptomic data (n=175) were downloaded from cBioPortal (https://www.cbioportal.org) along with their clinical annotation. Bailey dataset (n=91) were downloaded from its original study.
- Gene expression quantification of bulk RNA-seq data: To quantify gene and transcript expression, we used RSEM v1.2.25. We first generated reference transcriptome index files using “rsem-prepare-reference”. Then, “rsem-calculate-expression” was used to align paired-end reads of bulk RNA-seq data (FASTQ files) from 63 PDAC patients to the transcriptome to finally obtain gene and transcript expression profiles. TPM (transcripts per million)-normalized gene expression was used for further analysis with 1-augmented natural log transformation.
- Consensus clustering procedure: We performed 2-stage Leiden clustering using the pre-selected gene set. To select genes that are not dependent on specific dataset, we used SNUBH data as well as additional TCGA and Bailey data. Before clustering, we obtained the normalized variance of genes using the highly variable gene selection method (Satija, Farrell et al. 2015) in each data, and used the normalized variance to rank the genes. After that, we explored the number of genes and the number of dimensions to reduce through principal component analysis (PCA) for clustering. Then, we calculated the neighborhood graph (Becht, McInnes et al. 2018) using the dimension reduced features, and performed leiden clustering using the neighborhood graph. In the first clustering with the pre-selected gene set, we selected the clustering parameters where the gene expression of the clusters was similar in the three datasets. In the case of Tumor subtyping, we checked the inclusion of well-known gene such as GATA6. Finally, we extracted the common important genes from the three datasets in the second clustering, and obtained 110 genes for Tumor subtyping and 51 genes for TME subtyping. By using this gene set, we confirmed that the second clustering results showed more than 80% consistency with the initial clustering results (Supplemental Figure 2B). In particular, we set the clustering parameters to divide into two clusters for Tumor and TME subtyping, and were able to identify clusters with common characteristics in the three datasets. When we set the parameters to divide into more than two clusters, the clusters were divided in a dataset-dependent manner, making it difficult to find common points among the three datasets.

## Other data processing and statistical analysis
- Deconvolution analysis: We integrated information from various studies to compute functional scores based on gene expression. We identified marker genes for regulatory T cells (Tregs), myeloid-derived suppressor cells (MDSCs), immune suppression factors, and tumor-associated macrophages (TAMs), which are listed in Supplemental Table 3. For each patient sample, the sum of expression levels for the genes included in the signature was converted into a z-score, which served as a functional gene expression score. Additionally, we utilized the MCPcounter (R package) (Becht, Giraldo et al. 2016) to estimate microenvironment cell populations from bulk RNA-seq data.
- Single-sample GSEA: Normalized enrichment scores (NES) were obtained for each sample (patient) using the single sample GSEA function in GSEApy (Barbie, Tamayo et al. 2009, Abazeed, Adams et al. 2013, Fang, Liu et al. 2023). Four pathway databases, MSigDB_Hallmark_2020, KEGG_2021_Human, GO_Biological_Process_2021, and Reactome_2016, were used. Genes that showed significant differences (adjusted p-value<0.05) between C1 and C2 or between E1 and E2 were used to calculate the NES for each pathway. Then, the pathways that showed significant differences between subtypes (adjusted p-value<1e-5) were selected.
- Statistical tests: The chi-square test was used to check statistical significance of the associations between the pathologic stage of PDAC patients and TME subtype, whereas the McNemar's test was used for the associations between the tumor subtype and the TME subtype. The t-test was used to obtain the p-values in DEG analysis, ssGSEA and NES differences. Benjamini-Hochberg correction was then applied to obtain the false discovery rate (FDR). The log-rank test was used to check the significance in survival differences.

# Paper Summary Figure: Overview of Research Process and Results
## SNUBH PDAC Cohort
![Fig_](https://github.com/JoonghoLee/PDAC_Subtyping/assets/35910715/47036a01-561f-446b-bf69-0980f0c9e2bc)
![Fig_2png](https://github.com/JoonghoLee/PDAC_Subtyping/assets/35910715/f934c313-0d4b-416b-82ca-4a25aa1f1c44)
## PDAC Subtyping
![Fig1](https://github.com/JoonghoLee/PDAC_Subtyping/assets/35910715/07f39df5-6bfd-4e4b-980b-edcb6663e1bc)

# Result
## Tumor Subtype
![Fig2](https://github.com/JoonghoLee/PDAC_Subtyping/assets/35910715/57f328ce-6265-4624-8f55-12bfcbad3c47)
## TME Subtype
![Fig3](https://github.com/JoonghoLee/PDAC_Subtyping/assets/35910715/2b5f2875-57cc-4a1c-bfa3-d71e9e9b1110)
![Fig4](https://github.com/JoonghoLee/PDAC_Subtyping/assets/35910715/5bd0c092-fa3b-4d2a-8cd2-73c458f568f0)
## Combined Subtype
![Fig5](https://github.com/JoonghoLee/PDAC_Subtyping/assets/35910715/3ded9e69-ae7b-4053-8fa0-6b7bc884f4b5)

## Correlation between deconvolution tools 
![f1](https://github.com/JoonghoLee/PDAC_Subtyping/assets/35910715/d83ea731-3be0-42d2-8d2e-9d8f0666df25)
## FFX response
![f2](https://github.com/JoonghoLee/PDAC_Subtyping/assets/35910715/61e44051-d5e7-4543-a3ec-92d91223dde3)


































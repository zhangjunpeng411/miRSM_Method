# miRSM_Method
Identifying miRNA sponge modules using biclustering and regulatory scores

## Background
MicroRNA (miRNA) sponges with multiple tandem miRNA binding sequences can sequester miRNAs from their endogenous target mRNAs. Therefore, miRNA sponge acting as a decoy is extremely important for long-term lossof-function studies both in vivo and in silico. Recently, a growing number of in silico methods have been used as an effective technique to generate hypotheses for in vivo methods for studying the biological functions and regulatory mechanisms of miRNA sponges. However, most existing in silico methods only focus on studying miRNA sponge interactions or networks in cancer, the module-level properties of miRNA sponges in cancer is still largely unknown. We propose a novel in silico method, called miRSM (miRNA Sponge Module) to infer miRNA sponge modules in breast cancer. We apply miRSM to the breast invasive carcinoma (BRCA) dataset provided by The Cancer Genome Altas (TCGA), and make functional validation of the computational results. The functional validation results demonstrate that miRSM is a promising method to identify miRNA sponge modules and interactions, and may provide new insights for understanding the roles of miRNA sponges in cancer progression and development.

## Description of each file
BRCA_GENE: BRCA-related genes.

BRCA_miRNA: BRCA-related miRNAs

DEG_1E-04: Differentially expressed mRNAs.

DEmiRNA_1E-02.csv: Differentially expressed miRNAs.

Hallmark_GENE.csv: Hallmark-related genes.

HPRD.csv: Putative PPIs.

HTRI.csv: Putative TF-target interactions.

miRNA(453)-mRNA(11157).csv: Matched miRNA and mRNA expression data in BRCA.

miRTarBase_Strong.csv: Experimentally validated miRNA-mRNA inrteractions.

TargetScan_7.1.csv: Putative miRNA-mRNA interactions.

Validated_miRSponge.csv: Experimentally validated miRNA sponge inrteractions.

miRSM_Method.R: Utility functions for identifying miRNA sponge modules.

Main scripts.R: Running scripts for identifying miRNA sponge modules.

## The usage of miRSM_Method
Paste all files into a single folder (set the folder as the directory of R environment), the workflow of miRSM_Method is implemented in miRSM_Method.R. The users can simply run the scripts as follows.

```{r echo=FALSE, results='hide', message=FALSE}
source("Main scripts.R")
```

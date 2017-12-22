# TumorPurity_Prediction
Prediction of tumor purity, stromal and immune infiltration in ovarian cancer samples by implementing library estimate. 

### Details
A custom R script which takes in data from GSE datasets and pre-processes gene expression (.CEL) files to create a gene expression set using package affy was implemented. The script also uses library estimate to predict stromal and immune infiltration in ovarian cancer patient samples along with predicting tumor purity.

### Requirements
This script requires user to pre-install following R packages from bioconductor and R-Forge
1. Affy
2. R.utils
3. Estimate

### Installation Instructions
1. Affy
source("http://bioconductor.org/biocLite.R")
biocLite("affy")

2. R.utils
install.packages("R.utils")

3. estimate
install.packages("estimate", repos="http://R-Forge.R-project.org")

### Usage



### Reference
1. Yoshihara K, Shahmoradgoli M, Martínez E, Vegesna R, Kim H, Torres-Garcia W, Treviño V, Shen H, Laird PW, Levine DA, Carter SL, Getz G, Stemke-Hale K, Mills GB, Verhaak R. (2013)."Inferring tumour purity and stromal and immune cell admixture from expression data."
Nature Communications 4, 2612 doi:10.1038/ncomms3612

2. Zeeberg, Barry R, et al. “Mistaken Identifiers: Gene Name Errors Can Be Introduced Inadvertently When Using Excel in Bioinformatics.” BMC Bioinformatics, BioMed Central, 23 June 2004, bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-5-80.

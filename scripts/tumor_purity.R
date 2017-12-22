#!/usr/bin/Rscript

### Set command line arguments ###
options(warn=1)
args <- commandArgs(trailingOnly=TRUE)
if(length(args) == 2){
stop("No Arguments Provided... Please supply 2 argument...")
}

### Print error if incorrect arguments ###
input_GSE <- args[1]
probe_file <- args[2]
gene_date <- args[3]
extension <- grepl("^GSE", input_GSE)
if(!extension){
stop("Input file not defined, or bad file extension.")
}

### Download GSE dataset and untar to obtain .CEL.gz files ###
var <- gsub(".{3}$", "", input_GSE)
URL <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/xnnn/y/suppl/y_RAW.tar"
URL <- gsub("x", var, URL)
URL <- gsub("y", input_GSE, URL)
input_GSE <- paste(sep="", input_GSE, "_RAW.tar")
download.file(URL, destfile = input_GSE)

GSE <- gsub("_RAW.tar$", "", input_GSE)
dir.create(GSE)
untar(input_GSE, exdir = GSE)
setwd(GSE)
files <- list.files(".", pattern="*.CEL.gz", recursive=TRUE,full.names=TRUE)

### unzip to obtain .CEL files ###
library(R.utils)
#lapply(files, gunzip)
cels = list.files(".", pattern = "CEL", full.names=TRUE)

### Use package affy to preprocess CEL files ###
library(affy)
raw.data <- ReadAffy(filenames = cels, cdfname = "hgu133plus2cdf")
data.rma.norm <- rma(raw.data)
rma <- exprs(data.rma.norm)

### Merge by ProbeID ###
write.table(rma, file = "rma.txt", quote = FALSE, sep = "\t")
tt <- cbind(row.names(rma), rma)
colnames(tt) = c("ProbeID", sub(".cel", "", colnames(rma), ignore.case=TRUE))
rownames(tt)=NULL
probeInfo <- read.table(probe_file, header = TRUE, sep ="\t")
comb <- merge(probeInfo, tt, by.x = "ProbeID", by.y = "ProbeID")
comb <- as.data.frame(comb)

### Aggregate  commom GeneSymbols using Max values ###
mydata.max <- aggregate(. ~ GeneSymbol, data = comb, max)

DateToGene <- read.table(gene_date, header = TRUE, sep = "\t")

### Function to convert from factor to characters ###
Fac_To_Char <- function(x){
f <- sapply(x, is.factor)
x[f] <- lapply(x[f], as.character)
}
 
mydata.max$GeneSymbol <- Fac_To_Char(mydata.max$GeneSymbol)
DateToGene$Gene <- Fac_To_Char(DateToGene$Gene)
DateToGene$Date <- Fac_To_Char(DateToGene$Date)


### For loop to change date format to gene symbols ###
for(i in 1:length(mydata.max$GeneSymbol)){
for(j in 1:length(DateToGene$Date)){
if(as.character(mydata.max$GeneSymbol[i]) != as.character(DateToGene$Date[j])){
next
}
else if(as.character(mydata.max$GeneSymbol[i]) == as.character(DateToGene$Date[j])){
mydata.max$GeneSymbol[i] = DateToGene$Gene[j]
}
}
}

### Convert list attribute to character ###
mydata.max$GeneSymbol <- unlist(mydata.max$GeneSymbol)

### Save Dataframe as a text file ###
write.table(mydata.max, file = "/projects/home/kpatankar7/ESTIMATE_package/estimate/inst/extdata/GSE2.txt", sep = "\t", row.names = FALSE, quote = FALSE)


### Use package estimate to plot estimate scores for individual cel file ###
library(estimate)
OvarianCancerExpr <- "/projects/home/kpatankar7/ESTIMATE_package/estimate/inst/extdata/GSE2.txt"
filterCommonGenes(input.f=OvarianCancerExpr, output.f="OV_10412genes.gct", id="GeneSymbol")
estimateScore("OV_10412genes.gct", "OV_estimate_score.gct", platform="affymetrix")
plotPurity(scores="OV_estimate_score.gct", samples="all_samples", platform="affymetrix")


### Running Time ###
ptm <- proc.time()
proc.time()-ptm

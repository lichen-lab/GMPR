# GMPR
Calculate the normalizing factors for zeroinflated sequencing data such as microbiome sequencing data

GMPR.R - R function for calculating the GMPR size factors, will be slow if sample size is more than 1,000.

GMPR.Example.R - Illustration the use of GMPR size factors.


 **Installation**
---------------

Run following commands in R:

```
library(devtools)
install_github("lichen-lab/GMPR")
```

 **Example**
---------------

```
# Demonstrate the use of GMPR factor
require(GUniFrac)
require(vegan)
require(DESeq2)
library(GMPR)

data(throat.otu.tab)
data(throat.meta)

###########################################################################################################
# Calculate GMPR size factor
# Row - features, column - samples
otu.tab <- t(throat.otu.tab)
gmpr.size.factor <- GMPR(otu.tab)
###########################################################################################################

# Two potential applications of GMPR size factors

###########################################################################################################
# Application 1: Counts are normalized by size factors to reduce the variation due to different library sizes
# The normalized counts are subject to further downstream analysis such as ordination (PCA, PCoA), clustering,
# and other multivariate methods. Note that further data transformation such as VST transformation (DESeq2)
# may be needed in order to reveal patterns. Here shows an example of BC distance based ordination 
otu.tab.norm <- t(t(otu.tab) / gmpr.size.factor)
dist.mat <- vegdist(t(otu.tab.norm))
PCs <- cmdscale(dist.mat, k=2)
plot(PCs[, 1], PCs[, 2], col=factor(throat.meta$SmokingStatus))
###########################################################################################################

###########################################################################################################
# Application 2:  Differential abundance analysis using DESeq2 using GMPR size factors instead of the default
dds <- DESeqDataSetFromMatrix(countData = otu.tab,
		colData = throat.meta,
		design= ~ SmokingStatus)

# Replace size factor
# dds <- estimateSizeFactors(dds)
sizeFactors(dds) <-  gmpr.size.factor
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
results(dds, contrast=c("SmokingStatus", "NonSmoker", "Smoker"))
###########################################################################################################
```





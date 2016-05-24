## Required local file structure:
## (see http://aroma-project.org/setup for details)
## 
## rawData/
##   GSE29172/
##     GenomeWideSNP_6/
##       GSM721892_PA097_mix100_091027_GenomeWideSNP_6_.CEL
##       GSM721893_PA097_mix30_110329_GenomeWideSNP_6_.CEL
##       GSM721894_PA097_mix50_110329_GenomeWideSNP_6_.CEL
##       GSM721895_PA097_mix70_091027_GenomeWideSNP_6_.CEL
##   GSE26302/
##     GenomeWideSNP_6/
##       GSM645856_PA097_NCI-BL1395_091027_GenomeWideSNP_6_.CEL

library("aroma.affymetrix")
log <- verbose <- Arguments$getVerbose(-8, timestamp=TRUE)
options(digits=4) # Don't display too many decimals.

## Setup raw data set
cdf <- AffymetrixCdfFile$byChipType("GenomeWideSNP_6", tags="Full")
print(cdf)

dataSet <- "GSE29172"
csR <- AffymetrixCelSet$byName(dataSet, cdf=cdf)
print(csR)
resR <- doASCRMAv2(csR, verbose=verbose)

dataSet <- "GSE26302"
csN <- AffymetrixCelSet$byName(dataSet, cdf=cdf)
print(csN)
resN <- doASCRMAv2(csN, verbose=verbose)

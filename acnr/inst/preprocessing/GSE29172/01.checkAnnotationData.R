## Required local file structure:
## (see http://aroma-project.org/setup for details)
## 
## annotationData/
##   chipTypes/
##     GenomeWideSNP_6/
##       GenomeWideSNP_6,Full.CDF
##       GenomeWideSNP_6,Full,monocell.CDF
##       GenomeWideSNP_6,Full,na31,hg19,HB20110328.ufl
##       GenomeWideSNP_6,Full,na31,hg19,HB20110328.ugp
##       GenomeWideSNP_6,HB20080710.acs

library("aroma.affymetrix")
log <- verbose <- Arguments$getVerbose(-8, timestamp=TRUE)
options(digits=4) # Don't display too many decimals.

## Verifying annotation data files
cdf <- AffymetrixCdfFile$byChipType("GenomeWideSNP_6", tags="Full")
print(cdf)

gi <- getGenomeInformation(cdf)
print(gi)

si <- getSnpInformation(cdf)
print(si)

acs <- AromaCellSequenceFile$byChipType(getChipType(cdf, fullname=FALSE))
print(acs)

## Do not continue beyond this point if any of the above fails!

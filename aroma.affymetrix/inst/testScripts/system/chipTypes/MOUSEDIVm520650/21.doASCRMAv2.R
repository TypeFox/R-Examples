##########################################################################
# Allele-specific CRMAv2
##########################################################################
library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE27691";
chipType <- "MOUSEDIVm520650";

csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# AS-CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
res <- doASCRMAv2(csR, lengthRange=c(450,2000), drop=FALSE, verbose=verbose);
print(res);

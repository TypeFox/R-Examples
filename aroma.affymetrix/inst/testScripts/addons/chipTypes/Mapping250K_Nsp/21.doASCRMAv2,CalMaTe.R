##########################################################################
# Allele-specific CRMAv2 and CalMaTe
#
# Author: Henrik Bengtsson
# Created on: 2012-06-15
# Last updated: 2016-01-06
##########################################################################
library("aroma.affymetrix");
library("calmate");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE12702";
chipType <- "Mapping250K_Nsp";

csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);
csR <- csR[1:12]
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# AS-CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsList <- doASCRMAv2(csR, verbose=verbose);
print(dsList);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CalMaTe
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cmt <- CalMaTeCalibration(dsList);
print(cmt);

dsCList <- process(cmt, verbose=verbose);
print(dsCList);

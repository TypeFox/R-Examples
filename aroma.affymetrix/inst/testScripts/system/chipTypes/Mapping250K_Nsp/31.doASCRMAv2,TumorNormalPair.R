##########################################################################
# Allele-specific CRMAv2
#
# Author: Henrik Bengtsson
# Created on: 2011-11-11
# Last updated: 2016-01-06
##########################################################################
library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE34754";
chipType <- "Mapping250K_Nsp";

csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extract a single tumor-normal pair
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extract tumor-normal pair
pair <- c(T="GSM854620", N="GSM854615");
csR <- csR[indexOf(csR, pair)];
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# AS-CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsNList <- doASCRMAv2(csR, verbose=verbose);
print(dsNList);

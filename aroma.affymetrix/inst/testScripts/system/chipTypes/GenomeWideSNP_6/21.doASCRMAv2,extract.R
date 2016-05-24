##########################################################################
# Allele-specific CRMAv2
#
# Author: Henrik Bengtsson
# Created on: 2011-11-11
# Last updated: 2012-09-02
##########################################################################
library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE13372,testset";
chipType <- "GenomeWideSNP_6,Full";

csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# AS-CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
res <- doASCRMAv2(csR, drop=FALSE, verbose=verbose);
print(res);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# extractDataFrame();
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cesN <- res$cesN;
units <- c(1:4, 600+1:4, 1000+1:4, 1000000+1:4);  # A mix of unit types
data <- extractDataFrame(cesN, units=units, addNames=TRUE, verbose=verbose);
str(data);
print(data);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# extractTheta()
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
theta <- extractTheta(cesN, units=units, drop=TRUE, verbose=verbose);
str(theta);
print(theta);

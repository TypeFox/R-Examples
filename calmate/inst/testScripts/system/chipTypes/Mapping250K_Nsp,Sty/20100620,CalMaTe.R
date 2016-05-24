###########################################################################
# Author: Henrik Bengtsson
# Created on: 2010-05-18
# Last updated: 2010-09-28
#
# Data:
# rawData/GSE12702/Mapping250K_Nsp/
#
# URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE12702
###########################################################################

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CalMaTe
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
library("calmate");
verbose <- Arguments$getVerbose(-10, timestamp=TRUE);

dataSet <- "Affymetrix_2006-TumorNormal";
tags <- "ACC,-XY,BPN,-XY,RMA,FLN,-XY";
chipType <- "Mapping250K_Nsp";
dsT <- AromaUnitTotalCnBinarySet$byName(dataSet, tags=tags, chipType=chipType);
print(dsT);

dsB <- AromaUnitFracBCnBinarySet$byName(dataSet, tags=tags, chipType=chipType);
print(dsB);

stopifnot(identical(getNames(dsT), getNames(dsB)));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CalMaTe
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsList <- list(total=dsT, fracB=dsB);
asn <- CalMaTeCalibration(dsList);
print(asn);

ugp <- getAromaUgpFile(dsT);
chr <- 17;
units <- getUnitsOnChromosome(ugp, chr);

dsNList <- process(asn, units=units, verbose=verbose);
print(dsNList);

# Array #ii
ii <- 1;

# Extract raw (TCN,BAF)
df <- getFile(dsList$total, ii);
dfR <- getAverageFile(dsList$total, verbose=verbose);
gamma <- extractRawCopyNumbers(df, logBase=NULL, chromosome=chr);
gammaR <- extractRawCopyNumbers(dfR, logBase=NULL, chromosome=chr);
gamma <- 2*divideBy(gamma, gammaR);
df <- getFile(dsList$fracB, ii);
beta <- extractRawAlleleBFractions(df, chromosome=chr);

# Extract calibrated (TCN,BAF)
dfN <- getFile(dsNList$fracB, ii);
betaN <- extractRawAlleleBFractions(dfN, chromosome=chr);
dfN <- getFile(dsNList$total, ii);
gammaN <- extractRawCopyNumbers(dfN, logBase=NULL, chromosome=chr);

# Plot
subplots(4, ncol=2, byrow=FALSE);
plot(beta);
title(sprintf("%s", getName(beta)));
plot(gamma);
plot(betaN);
title(sprintf("%s (CalMaTe)", getName(betaN)));
plot(gammaN);


###########################################################################
# HISTORY:
# 2010-09-28
# o Added total CN track.
# o Using CalMaTeCalibration formely known CalMaTeNormalization.
# 2010-06-20
# o Created.
###########################################################################

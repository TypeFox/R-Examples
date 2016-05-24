###########################################################################
# Author: Henrik Bengtsson
# Created on: 2012-02-06
# Last updated: 2012-02-06
#
# Data:
# rawData/GSE8605/Mapping10K_Xba142/
#
# URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE12702
###########################################################################

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CalMaTe
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
library("calmate");
verbose <- Arguments$getVerbose(-10, timestamp=TRUE);

dataSet <- "GSE8605";
chipType <- "Mapping10K_Xba142";
tags <- "ACC,-XY,BPN,-XY,RMA,FLN,-XY";
dsT <- AromaUnitTotalCnBinarySet$byName(dataSet, tags=tags, chipType=chipType);
dsB <- AromaUnitFracBCnBinarySet$byName(dataSet, tags=tags, chipType=chipType);
dsList <- list(total=dsT, fracB=dsB);
print(dsList);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CalMaTe
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
asn <- CalMaTeCalibration(dsList, fB1=0.33, fB2=0.66);
print(asn);

dsNList <- process(asn, verbose=verbose);
print(dsNList);

stop();

ugp <- getAromaUgpFile(dsT);
chr <- 17;
units <- getUnitsOnChromosome(ugp, chr);


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

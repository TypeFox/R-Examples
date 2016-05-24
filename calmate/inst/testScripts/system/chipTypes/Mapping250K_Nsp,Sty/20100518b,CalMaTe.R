###########################################################################
# Author: Henrik Bengtsson
# Created on: 2010-05-18
# Last updated: 2010-05-18
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
library("aroma.core");
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
# Study array #1 and chromosome 22
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
array <- 1;
chr <- 17;

chrTag <- sprintf("Chr%02d", chr);

# Extract the array
dfT <- getFile(dsT, array);
dfB <- getFile(dsB, array);
stopifnot(identical(getName(dfT), getName(dfB)));

# Extract the chromosome
ugp <- getAromaUgpFile(dsT);
units <- getUnitsOnChromosome(ugp, chr);
pos <- getPositions(ugp, units=units);

# Extract (total, fracB) signals
total <- extractMatrix(dsT, units=units);
fracB <- extractMatrix(dsB, units=units);
dim <- c(nrow(total), 2, ncol(total));
dimnames <- list(rownames(total), c("total", "fracB"), colnames(total));
naValue <- as.double(NA);
data <- array(naValue, dim=dim, dimnames=dimnames);
data[,"total",] <- total;
data[,"fracB",] <- fracB;

dataC <- calmateByTotalAndFracB(data, verbose=-8);

# Plot
ylim <- c(-0.2, 1.2);
subplots(2, ncol=1);

# Array #1
ii <- 1;
name <- dimnames(dataC)[[3]][ii];

fracB <- RawAlleleBFractions(data[,"fracB",ii], x=pos);
plot(fracB, pch=".", ylim=ylim);
stext(side=3, pos=0, name);
stext(side=3, pos=1, chrTag);

fracBC <- RawAlleleBFractions(dataC[,"fracB",ii], x=pos);
plot(fracBC, pch=".", ylim=ylim);
stext(side=3, pos=0, name);
stext(side=3, pos=1, chrTag);


###########################################################################
# HISTORY:
# 2010-05-18
# o Created.
###########################################################################

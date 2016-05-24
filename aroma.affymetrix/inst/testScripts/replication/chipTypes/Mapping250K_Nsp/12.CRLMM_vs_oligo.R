###########################################################################
# Replication test
#
# Description:
# This test verifies that aroma.affymetrix can reproduce the CRLMM
# genotype estimates as estimated by oligo.
# It verifies that they give the same results whether or not one
# is normalizing towards the HapMap reference (as defined by oligo).
#
# Author: Henrik Bengtsson
# Created: 2008-12-04
# Last modified: 2012-09-01
###########################################################################
library("aroma.affymetrix");
library("oligo");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE34754";
chipType <- "Mapping250K_Nsp";

cdf <- AffymetrixCdfFile$byChipType(chipType);
csR <- AffymetrixCelSet$byName(dataSet, cdf=cdf);
print(csR);

# Assert that the correct Platform Design package is installed
pdPkgName <- cleanPlatformName(chipType);
library(pdPkgName, character.only=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# SNPRMA according to aroma.affymetrix
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ces <- justSNPRMA(csR, normalizeToHapmap=TRUE, returnESet=FALSE, verbose=verbose);
print(ces);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CRLMM according to aroma.affymetrix
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
recalibrate <- TRUE;
crlmm <- CrlmmModel(ces, tags="*,oligo", recalibrate=recalibrate);
print(crlmm);

units <- fit(crlmm, ram="oligo", verbose=verbose);
str(units);

callSet <- getCallSet(crlmm);
print(callSet);

confSet <- getConfidenceScoreSet(crlmm);
print(confSet);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CRLMM according to oligo
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create temporary output directory for oligo
path <- file.path("oligoData", dataSet, chipType);
Arguments$getWritablePath(dirname(path)); # Create parent directories
if (!isDirectory(path)) {
  oligo:::justCRLMMv2(getPathnames(csR), tmpdir=path, recalibrate=recalibrate, balance=1.5, verbose=as.logical(verbose));
}
# Files created by oligo
print(list.files(path));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Compare genotype calls
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
calls0 <- readSummaries("calls", path);

unitNames <- rownames(calls0);
units <- indexOf(cdf, names=unitNames);
calls <- extractGenotypes(callSet, units=units, encoding="oligo");

count <- 0;
for (cc in 1:ncol(calls)) {
  idxs <- whichVector(calls[,cc] != calls0[,cc]);
  count <- count + length(idxs);
  cat(sprintf("%s: ", getNames(callSet)[cc]));
  if (length(idxs) > 0) {
    map <- c("AA", "AB", "BB");
    cat(paste(map[calls[idxs,cc]], map[calls0[idxs,cc]], sep="!="), sep=", ");
  }
  cat("\n");
}
cat(sprintf("Averages number of discrepances per array: %.1f\n", count/ncol(calls)));
errorRate <- count/length(calls);
cat(sprintf("Concordance rate: %.5f%%\n", 100*(1-errorRate)));
stopifnot(errorRate < 5e-4);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Compare confidence scores
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
conf <- extractMatrix(confSet, units=units);
conf0 <- readSummaries("conf", path);

res <- all.equal(conf, conf0, check.attributes=FALSE);
print(res);

# Sanity check
stopifnot(all.equal(conf, conf0, check.attributes=FALSE, tolerance=1e-4));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Visual comparison of confidence scores
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
toPNG(getFullName(ces), tags="CRLMM_vs_oligo", aspectRatio=0.8, {
  subplots(ncol(conf));
  par(mar=c(3,2,1,1)+0.1);
  lim <- c(0,1);
  for (cc in 1:ncol(conf)) {
    plot(NA, xlim=lim, ylim=lim);
    abline(a=0, b=1, col="#999999");
    points(conf[,cc], conf0[,cc], pch=".", cex=3);
    rho <- cor(conf[,cc], conf0[,cc]);
    stext(side=3, pos=0, line=-1, sprintf("rho=%.4f", rho));
    stext(side=3, pos=0, getNames(confSet)[cc]);
    cat(sprintf("Array #%d: Correlation: %.4f\n", cc, rho));
  }
});



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot BAFs along genome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ugp <- getAromaUgpFile(cdf);
chr <- 2;
units <- getUnitsOnChromosome(ugp, chromosome=chr);
pos <- getPositions(ugp, units=units)/1e6;

chrTag <- sprintf("Chr%02d", chr);

toPNG(getFullName(ces), tags=c(chrTag, "TCN,BAF"), aspectRatio=0.5*nbrOfArrays(ces), {
  layout(matrix(seq_along(ces), ncol=1));
  par(mar=c(3.5,4,1.5,1), mgp=c(1.8,0.5,0), pch=".");

  for (ii in seq_along(ces)) {
    ce <- ces[[ii]];
    gc <- callSet[[ii]];
    sampleName <- getName(ce);

    data <- extractTotalAndFracB(ce, units=units, drop=TRUE);
    colnames(data) <- c("total", "fracB");
    calls <- extractGenotypes(gc, units=units, drop=TRUE);

    col <- c(AA=1, AB=2, BB=3)[calls];
    plot(pos, data[,"fracB"], col=col, pch=".", cex=4, ylim=c(0,1));
    stext(side=3, pos=0, sampleName);
    stext(side=3, pos=1, chrTag);
  } # for (ii ...)
});

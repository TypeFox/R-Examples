###########################################################################
# Replication test
#
# Description:
# This test verifies that aroma.affymetrix can reproduce the estimates
# of the MAT (Model-based Analysis of Tiling arrays) algorithm.
#
# Author: Mark Robinson and Henrik Bengtsson
# Created: 2008-12-09
# Last modified: 2012-09-02
###########################################################################
library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-20, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE24546";
tags <- "testset";
chipType <- "Hs_PromPR_v02";
sampleNamesMap <- c(
  GSM605951="Prec1_MeDNA_Input1",
  GSM605952="Prec1_MeDNA_IP2",
  GSM605953="Prec1_MeDNA_IP1"
);

cdf <- AffymetrixCdfFile$byChipType(chipType);
print(cdf);

csR <- AffymetrixCelSet$byName(dataSet, tags=tags, cdf=cdf);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Normalize the data using the MAT model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
mn <- MatNormalization(csR);
print(mn);

csM <- process(mn, verbose=more(verbose, 3));
print(csM);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Convert data set such that it maps to the "unique" CDF
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csU <- convertToUnique(csM, verbose=verbose);
print(csU);

# Rename
setFullNamesTranslator(csU, function(names, ...) { sampleNamesMap[names] });


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Validate the "unique" estimates
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csList <- list(csM, csU);
for (ii in seq_along(csU)) {
  yList <- lapply(csList, FUN=function(cs) {
    cf <- cs[[ii]];
    cdf <- getCdf(cf);
    units <- indexOf(cdf, "^chr1F");
    cells <- getCellIndices(cdf, units=units, stratifyBy="pm",
                                            unlist=TRUE, useNames=FALSE);
    y <- extractMatrix(cf, cells=cells, drop=TRUE);
    y;
  });
  stopifnot(all.equal(yList[[1]], yList[[2]]));
  rm(yList);
} # for (ii ...)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot ratios estimates along a single chromosome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calculate the average
cfR <- getAverageFile(csU);

cdfU <- getCdf(csU);
acpU <- AromaCellPositionFile$byChipType(getChipType(cdfU));

# Identify cells on Chr21
chr <- 21;
chrTag <- sprintf("Chr%02d", chr);
cells <- whichVector(acpU[,1,drop=TRUE] == chr);
str(cells);

# Get the chromosomal positions of these cells
pos <- acpU[cells,2,drop=TRUE];

yR <- extractMatrix(cfR, cells=cells, drop=TRUE, verbose=verbose);
for (ii in seq_along(csU)) {
  cf <- csU[[ii]];
  sampleName <- getName(cf);

  toPNG(getFullName(csU), tags=c(sampleName, "tracks"), aspectRatio=1/2, {
    par(mar=c(3,3,1,1)+0.1, mgp=c(1.8,0.5,0));
    xlab <- "Physical position (Mb)";
    ylab <- expression(log2(PM));

    x <- pos / 1e6;
    y <- extractMatrix(cf, cells=cells, drop=TRUE, verbose=verbose);
    M <- log2(y/yR);

    plot(x, M, pch=".", cex=2, xlab=xlab, ylab=ylab);
    stext(side=3, pos=0, cex=0.7, sampleName);
    stext(side=3, pos=1, cex=0.7, chrTag);
    stext(side=4, pos=0, cex=0.7, getFullName(csU));
    stext(side=4, pos=1, cex=0.7, sprintf("n=%d", length(y)));
  });
} # for (ii ...)

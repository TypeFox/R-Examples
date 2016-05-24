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
library("matrixStats"); # colMedians()
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
# Compare to precomputed estimates from external MAT software
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Compare common units with prefix "chr1F".
cdfU <- getCdf(csU);
units <- indexOf(cdfU, "^chr1F");
cells <- getCellIndices(cdfU, units=units, stratifyBy="pm",
                        unlist=TRUE, useNames=FALSE);

# Get the chromosomal positions of these cells
acp <- AromaCellPositionFile$byChipType(getChipType(cdfU));
pos <- acp[cells,2,drop=TRUE];

# Order cells by chromsomal position
o <- order(pos);
pos <- pos[o];
cells <- cells[o];


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Compare to the MAT software
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MAT normalized signals
normFile <- "Prec1_MeDNA_600_IP1-Input.tsv";
tsv <- TabularTextFile(filename=normFile, path=getPath(csR));
print(tsv);

# Sample names ran on external MAT software
sampleNames <- getColumnNames(tsv);
sampleNames <- setdiff(sampleNames, c("Chr", "Pos"));
sampleNames <- gsub("[.]CEL$", "", sampleNames);

# Sanity check
stopifnot(all(is.element(sampleNames, getNames(csU))));

# Corresponding aroma subset
csUt <- csU[sampleNames];


# Read MAT signals
colClasses <- c("Chr"="character", "Pos"="integer", rep("numeric", times=length(sampleNames)));
names(colClasses)[-(1:2)] <- sampleNames;
data <- readDataFrame(tsv, colClasses=colClasses, nrow=435000);
data <- subset(data, (Chr == "chr1" & Pos %in% pos));

# Order as (pos,yN)
o <- match(data$Pos, pos);
# Sanity check
stopifnot(all(is.finite(o)));
data <- data[o,,drop=FALSE];

# Extract signals
YMAT <- as.matrix(data[,-(1:2)]);
rownames(YMAT) <- NULL;
colnames(YMAT) <- sampleNames;

Y <- extractMatrix(csUt, cells=cells);
Y <- log2(Y);

# Sanity check
stopifnot(identical(dim(Y), dim(YMAT)));


# Graphical comparison
for (sampleName in sampleNames) {
  toPNG(getFullName(csU), tags=c(sampleName, "MatNormalization_vs_MAT"), {
    xlab <- expression(log[2](y[MAT]));
    ylab <- expression(log[2](y));
    plot(YMAT[,sampleName], Y[,sampleName], pch=".", xlab=xlab, ylab=ylab);
    abline(a=0, b=1, col="red", lwd=2, lty=3);
    stext(side=3, pos=0, sampleName);
    stext(side=3, pos=1, "Chr01");
    stext(side=4, pos=1, sprintf("n=%d", nrow(Y)));
    stext(side=4, pos=0, getFullName(csU));
  });
}

# Numerical comparison
D2 <- (Y - YMAT)^2;
avgDiff <- colMedians(D2, na.rm=TRUE);
names(avgDiff) <- colnames(D2);
cat("Average differences:\n");
print(avgDiff);

# Sanity check
stopifnot(all(avgDiff < 0.001));

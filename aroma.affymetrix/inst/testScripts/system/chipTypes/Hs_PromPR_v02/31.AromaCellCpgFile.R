###########################################################################/**
# File type test
#
# Description:
# This test verifies that aroma.affymetrix can create AromaCellCpgFile
# and AromaCellPositionFile objects for the (promoter) tiling array.
#
# Author: Mark Robinson
# Created: 2010-02-22
# Last modified: 2012-09-02
#*/###########################################################################
library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup the chip type
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- AffymetrixCdfFile$byChipType("Hs_PromPR_v02");
print(cdf);

cdfU <- getUniqueCdf(cdf);
print(cdfU);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Test allocation, writing and reading of 'acp' object
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allocate
acp <- AromaCellPositionFile$allocateFromCdf(cdfU, tags=c("*", "test"), 
                                            overwrite=TRUE, verbose=verbose);
print(acp);

nbrOfRandomCells <- 20;
cells <- sample(nbrOfCells(cdfU), size=nbrOfRandomCells);
chr <- sample(1:22, size=nbrOfRandomCells);
pos <- sample(1:1e6, size=nbrOfRandomCells);

# Write
acp[cells,1] <- chr;
acp[cells,2] <- pos;

# Read
cp <- acp[cells,];

# Sanity check
stopifnot((cp[,1] == chr) & (cp[,2] == pos));

# Clean up
rm(chr, pos, cp);
file.remove(getPathname(acp));

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Test allocation, writing and reading of 'acc' object
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allocate
acc <- AromaCellCpgFile$allocateFromCdf(cdfU, tags=c("*", "test"), 
                                            overwrite=TRUE, verbose=verbose);
print(acc);

cells <- sample(nbrOfCells(cdfU), size=nbrOfRandomCells);
cpg <- rnorm(nbrOfRandomCells);

# Write
acc[cells,1] <- 2^cpg;

# Read
cpg2 <- acc[cells,1,drop=TRUE];

d <- (log2(cpg2) - cpg)^2;
print(summary(d));

# Sanity check
stopifnot(sum(d) < 1e-8);

# Clean up
file.remove(getPathname(acc));

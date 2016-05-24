library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-4, timestamp=TRUE);

dataSet <- "GSE8605";
chipType <- "Mapping10K_Xba142";
res <- doCRMAv2(dataSet, chipType=chipType, drop=FALSE, verbose=verbose);
print(res);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CBS on a CnChipEffectSet
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cesN <- res$cesN;
segA <- CbsModel(cesN, tags="*,cesN");
print(segA);

# Try to segment
fit(segA, arrays=1:2, chromosomes=19, verbose=verbose);

# Write to file
pathnameA <- writeRegions(segA, arrays=1:2, chromosomes=19, skip=FALSE, verbose=verbose);
print(pathnameA);

segTableA <- read.table(pathnameA, header=TRUE);
print(segTableA);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CBS on an AromaUnitTotalCnBinarySet
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsNList <- res$dsNList;
dsT <- dsNList; # Hmm... should this have a different name?
segB <- CbsModel(dsT);
print(segB);

# Try to segment
fit(segB, arrays=1:2, chromosomes=19, verbose=verbose);

# Write to file
pathnameB <- writeRegions(segB, arrays=1:2, chromosomes=19, skip=FALSE, verbose=verbose);
print(pathnameB);

segTableB <- read.table(pathnameB, header=TRUE);
print(segTableB);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Compare
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Sanity check
stopifnot(all.equal(segTableB, segTableA));

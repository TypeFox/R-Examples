library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-4, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE8605";
chipType <- "Mapping10K_Xba142";
csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType, verbose=verbose);
print(csR);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (a) CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
res <- doCRMAv2(csR, drop=FALSE, verbose=verbose);
print(res);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (b) CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Process only the first six arrays (and in reverse order)
subset <- c(2,5,3);
dsN <- doCRMAv2(csR, arrays=subset, verbose=verbose);
print(dsN);


# Sanity checks
# Assert same output as input arrays (in the same order)
sampleNames <- getNames(csR)[subset];
stopifnot(getNames(dsN) == sampleNames);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (c) CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Process only the first six arrays (and in reverse order)
subset <- c(2,5,3);
dsN2 <- doCRMAv2(dataSet, chipType=chipType, arrays=subset, verbose=verbose);
print(dsN2);

# Sanity checks
stopifnot(equals(dsN2, dsN));

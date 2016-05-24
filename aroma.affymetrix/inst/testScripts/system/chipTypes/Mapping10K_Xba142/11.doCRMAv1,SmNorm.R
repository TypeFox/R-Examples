library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-4, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE8605";
chipType <- "Mapping10K_Xba142";
csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType, verbose=verbose);

# Process only the first six arrays (and in reverse order)
subset <- 6:1;
csR <- csR[subset];
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CRMAv1
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
res <- doCRMAv1(csR, drop=FALSE, verbose=verbose);
print(res);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Scale normalization test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ces <- res$ces;
sn <- ScaleNormalization3(ces);
print(sn);

cesN <- process(sn, verbose=verbose);
print(cesN);

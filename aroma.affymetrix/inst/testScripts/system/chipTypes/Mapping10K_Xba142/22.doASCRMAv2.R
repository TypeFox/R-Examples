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
# (a) AS-CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Process arrays in random order
subset <- sample(seq_along(csR));
res <- doASCRMAv2(csR, arrays=subset, drop=FALSE, verbose=verbose);
print(res);


# Sanity checks
# Assert same output as input arrays (in the same order)
sampleNames <- getNames(csR)[subset];
for (key in names(res)[-1]) {
  ds <- res[[key]];
  if (!inherits(ds, "GenericDataFileSet")) next;
  stopifnot(getNames(ds) == sampleNames);
}



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (b) CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Process only the first six arrays (and in reverse order)
subset <- c(2,5,3);
dsNList <- doASCRMAv2(csR, arrays=subset, verbose=verbose);
print(dsNList);


# Sanity checks
# Assert same output as input arrays (in the same order)
sampleNames <- getNames(csR)[subset];
stopifnot(getNames(dsNList$total) == sampleNames);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (c) CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Process only the first six arrays (and in reverse order)
subset <- c(2,5,3);
dsNList2 <- doASCRMAv2(dataSet, chipType=chipType, arrays=subset, verbose=verbose);
print(dsNList2);

# Sanity checks
for (key in names(dsNList2)) {
  stopifnot(equals(dsNList2[[key]], dsNList[[key]]));
}

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
sampleNames <- getNames(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CRMAv1
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
res <- doCRMAv1(csR, drop=FALSE, verbose=verbose);
print(res);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Sanity checks
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Assert same output as input arrays (in the same order)
for (key in names(res)[-1]) {
  ds <- res[[key]];
  if (!inherits(ds, "GenericDataFileSet")) next;
  stopifnot(getNames(ds) == sampleNames);
}

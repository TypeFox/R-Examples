library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "HapMap,CEU,testset";
chipType <- "Mapping50K_Hind240";

cdf <- AffymetrixCdfFile$byChipType(chipType);
csR <- AffymetrixCelSet$byName(dataSet, cdf=cdf);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# AS-CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
res <- doASCRMAv2(csR, drop=FALSE, verbose=verbose);
print(res);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CRLMM
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ces <- res$cesN;
print(ces);

# This should throw "Exception: Cannot fit CRLMM model: CRLMM 
# requires that the probe-level model was fitted without
# merging the strands (mergeStrands=FALSE)."
tryCatch({
  crlmm <- CrlmmModel(ces, tags="*,oligo", recalibrate=TRUE);
}, error = function(ex) {
  print(ex);
})

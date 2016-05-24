library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-10, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE37861";
chipType <- "MoGene-1_0-st-v1";

cdf <- AffymetrixCdfFile$byChipType(chipType, tags="r3");
print(cdf);

csR <- AffymetrixCelSet$byName(dataSet, cdf=cdf);
csR <- csR[1:6];
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# gcRMA-style background correction
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Currently, you must use the standard CDF file.
bc <- GcRmaBackgroundCorrection(csR, type="affinities");
print(bc);

# This should throw "Error: Cannot perform GCRMA background
# (type="affinities") correction: The number (0) of negative
# control is too small."
tryCatch({
  csB <- process(bc, verbose=verbose);
}, error = function(ex) {
  print(ex);
})

############################################################################
# ROBUSTNESS TEST
############################################################################
library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-4, timestamp=TRUE);

dataSet <- "GSE8605";
chipType <- "Mapping10K_Xba142";

csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);
csR <- csR[1:6];
print(csR);

cesN <- doASCRMAv1(csR, verbose=verbose, drop=FALSE)$cesN;

# This should throw "Exception: Unsupported chip effects. ..."
cbs <- NULL;
tryCatch({
  cbs <- CbsModel(cesN);
}, error = function(ex) {
  print(ex);
})
print(cbs);

# Sanity check
stopifnot(is.null(cbs));


############################################################################
# HISTORY:
# 2012-08-30
# o Simplified test script to utilize doASCRMAv1().
# 2009-11-13
# o Created.  Thanks to Pierre Neuvial for the report.
############################################################################

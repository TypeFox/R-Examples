###########################################################################
# Replication test
#
# Description:
# This test verifies that aroma.affymetrix can reproduce the SNPRMA 
# chip-effect estimates as estimated by oligo.
# It verifies that they give the same results whether or not one
# is normalizing towards the HapMap reference (as defined by oligo).
#
# Author: Henrik Bengtsson
# Created: 2008-12-04
# Last modified: 2010-09-01
###########################################################################
library("aroma.affymetrix");
library("oligo");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
compareESets <- function(eSet1, eSet2, tolerance=1e-4) {
  # Default
  res <- all.equal(eSet1, eSet2, tolerance=tolerance);
  if (!isTRUE(res)) return(res);

  # Compare assayData()
  data1 <- assayData(eSet1);
  data2 <- assayData(eSet2);
  data1 <- as.list(data1);
  data2 <- as.list(data2);
  res <- all.equal(data1, data2, tolerance=tolerance);
  if (!isTRUE(res)) return(res);
  TRUE;
} # compareESets()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "HapMap,CEU,testset";
chipType <- "Mapping50K_Hind240";

cdf <- AffymetrixCdfFile$byChipType(chipType);
csR <- AffymetrixCelSet$byName(dataSet, cdf=cdf);
print(csR);

# Assert that the correct Platform Design package is installed
pdPkgName <- cleanPlatformName(chipType);
library(pdPkgName, character.only=TRUE);

for (normalizeToHapmap in c(FALSE, TRUE)) {
  verbose && enter(verbose, sprintf("justSNPRMA(..., normalizeToHapMap=%s)", normalizeToHapmap));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # SNPRMA according to aroma.affymetrix
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  eSet <- justSNPRMA(csR, normalizeToHapmap=normalizeToHapmap, verbose=verbose);
  print(eSet);
  
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # SNPRMA according to oligo
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  eSet0 <- justSNPRMA(getPathnames(csR), normalizeToHapmap=normalizeToHapmap, verbose=as.logical(verbose));
  print(eSet0);
  
  if (!normalizeToHapmap) {
    # CLEAN UP: justSNPRMA() stores a target distribution file
    # in the working directory that we don't need
    filename <- sprintf("%s.quantileReference.rda", pdPkgName);
    if (isFile(filename)) {
      file.remove(filename);
    }
  }

 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Compare
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- compareESets(eSet, eSet0);
  verbose && print(verbose, res);

  # Sanity check
  stopifnot(compareESets(eSet, eSet0));

  verbose && exit(verbose);
} # for (normalizeToHapmap ...)

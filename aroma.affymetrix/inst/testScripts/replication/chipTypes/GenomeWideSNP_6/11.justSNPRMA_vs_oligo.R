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
dataSet <- "GSE13372,testset";
chipType <- "GenomeWideSNP_6";

cdf <- AffymetrixCdfFile$byChipType(chipType);
csR <- AffymetrixCelSet$byName(dataSet, cdf=cdf);

# Process only a subset of the arrays.  Since this data set
# contains many replicates (cf. GEO), they need to be for
# different samples.
sampleNamesMap <- c(
  GSM337641="HCC1143_GLEYS_A02",
#  GSM337646="HCC1143_TRIBE_H11",
  GSM337662="HCC1143BL_GLEYS_A01",
#  GSM337666="HCC1143BL_TRIBE_D02",
#  GSM337668="HCC1143BL_GHATS_H04",
#  GSM337674="HCC1143BL_TRIGS_G07",
  GSM337683="HCC1954_GLEYS_B02",
#  GSM337688="HCC1954_TRIBE_G12",
  GSM337696="HCC1954BL_GLEYS_B01",
#  GSM337700="HCC1954BL_TRIBE_B01",
#  GSM337702="HCC1954BL_GHATS_G10",
#  GSM337703="HCC1954BL_TRIGS_G11",
  GSM337707="NCI-H2347",
  GSM337708="NCI-H2347BL"
);
sampleNames <- names(sampleNamesMap);

sampleNames <- rev(sampleNames);
csR <- csR[sampleNames];
setFullName(csR, sprintf("%s,crlmmSubset", dataSet));
print(csR);

# Assert that the file header of the first CEL file in
# truly a GenomeWideSNP_6 (and not GenomeWideEx_6 because
# oligo::justSNPRMA() cannot handle chip type aliases)
stopifnot(getHeader(csR[[1]])$chiptype == chipType);

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
  res <- compareESets(eSet, eSet0, tolerance=0.02);
  verbose && print(verbose, res);

  # Sanity check
  stopifnot(res);

  verbose && exit(verbose);
} # for (normalizeToHapmap ...)

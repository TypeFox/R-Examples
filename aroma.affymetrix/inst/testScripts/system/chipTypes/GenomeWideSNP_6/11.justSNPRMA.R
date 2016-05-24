library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

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


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# SNPRMA according to aroma.affymetrix
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
eSet <- justSNPRMA(csR, normalizeToHapmap=TRUE, verbose=verbose);
print(eSet);

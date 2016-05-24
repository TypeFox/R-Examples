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
# RMA background correction
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Currently, you must use the standard CDF file.
bc <- RmaBackgroundCorrection(csR);
print(bc);

csB <- process(bc, verbose=verbose);
print(csB);


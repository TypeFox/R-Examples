library("aroma.affymetrix");
library("ACNE");
verbose <- Arguments$getVerbose(-4, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "HapMap,CEU,testset";
chipType <- "Mapping50K_Hind240";

csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# ACNE
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsNList <- doACNE(csR, fln=TRUE, verbose=verbose);
print(dsNList);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segment TCN
dsT <- dsNList$total;
seg <- CbsModel(dsT);
print(seg);

fit(seg, arrays=1, chromosomes=19, verbose=verbose);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# ChromosomeExplorer
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ce <- ChromosomeExplorer(seg);
print(ce);
process(ce, arrays=1:2, chromosomes=c(19,23), verbose=verbose);

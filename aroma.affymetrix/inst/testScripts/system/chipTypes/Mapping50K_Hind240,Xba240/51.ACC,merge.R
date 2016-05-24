library("aroma.affymetrix");
log <- Arguments$getVerbose(-4, timestamp=TRUE);


dataSet <- "HapMap,CEU,testset";
chipType <- "Mapping50K_Hind240";


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tests for setting up CEL sets and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allelic cross-talk calibration
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
acc <- AllelicCrosstalkCalibration(csR, mergeShifts=FALSE);
print(acc);
csC <- process(acc, verbose=log);
print(csC);

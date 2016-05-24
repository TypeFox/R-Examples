library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-50, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE8605";
chipType <- "Mapping10K_Xba142";

csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allelic-crosstalk calibration with and without mergeShifts
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
acc <- AllelicCrosstalkCalibration(csR, mergeShifts=TRUE);
print(acc);

csC <- process(acc, verbose=verbose);
print(csC);


acc <- AllelicCrosstalkCalibration(csR, mergeShifts=FALSE);
print(acc);

csC <- process(acc, verbose=verbose);
print(csC);

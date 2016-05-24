if (interactive()) savehistory();
library("aroma.affymetrix");
library("ACNE");

# - - - - - - - - - - - - - - - - - - - - - - -
# setup dataset and chip names
# - - - - - - - - - - - - - - - - - - - - - - -
log <- Arguments$getVerbose(-10, timestamp=TRUE);

dataSetName <- "HapMap270,100K,CEU,5trios"
chipType <- "Mapping50K_Hind240"

# - - - - - - - - - - - - - - - - - - - - - - -
# Setup annotation data
# - - - - - - - - - - - - - - - - - - - - - - -
cdf <- AffymetrixCdfFile$byChipType(chipType);
print(cdf);
gi <- getGenomeInformation(cdf);
print(gi);
si <- getSnpInformation(cdf);
print(si);


# - - - - - - - - - - - - - - - - - - - - - - -
# Setup data set
# - - - - - - - - - - - - - - - - - - - - - - -
csR <- AffymetrixCelSet$byName(dataSetName, cdf=cdf);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - -
# Calibrate and normalize
# - - - - - - - - - - - - - - - - - - - - - - -
acc <- AllelicCrosstalkCalibration(csR, model="CRMAv2");
print(acc);
csC <- process(acc, verbose=log);
print(csC);

bpn <- BasePositionNormalization(csC, target="zero");
print(bpn);
csN <- process(bpn, verbose=log);
print(csN);


# - - - - - - - - - - - - - - - - - - - - - - -
# Summarize replicated probes
# - - - - - - - - - - - - - - - - - - - - - - -
plm <- NmfSnpPlm(csN, mergeStrands=TRUE);
print(plm);
fit(plm, verbose=log);

ces <- getChipEffectSet(plm);
print(ces);


# - - - - - - - - - - - - - - - - - - - - - - -
# Fragment-length normalization
# - - - - - - - - - - - - - - - - - - - - - - -
fln <- FragmentLengthNormalization(ces, target="zero");
print(fln);
cesN <- process(fln, verbose=log);
print(cesN);



# - - - - - - - - - - - - - - - - - - - - - - -
# Segmentation
# - - - - - - - - - - - - - - - - - - - - - - -
# CBS needs theta = thetaA + thetaB
as <- AlleleSummation(cesN);
cesT <- process(as, verbose=log);
cbs <- CbsModel(cesT);
print(cbs);
ce <- ChromosomeExplorer(cbs);
print(ce);

process(ce, chromosomes=1:5, verbose=log);

library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-50, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "Affymetrix-CytoSampleData";
chipType <- "Cytogenetics_Array";

cdf <- AffymetrixCdfFile$byChipType(chipType);
print(cdf);

csR <- AffymetrixCelSet$byName(dataSet, cdf=cdf);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allelic-crosstalk calibration
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
acc <- AllelicCrosstalkCalibration(csR);
print(acc);
csC <- process(acc, verbose=verbose);
print(csC);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot allele pairs before and after calibration
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for (what in c("input", "output")) {
  toPNG(getFullName(acc), tags=c("allelePairs", what), aspectRatio=0.7, {
    plotAllelePairs(acc, array=1, what=what, verbose=verbose);
  });
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Probe-level summarization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plm <- AvgSnpPlm(csC, mergeStrands=TRUE);
print(plm);
if (length(findUnitsTodo(plm)) > 0) {
  # Fit CN probes quickly (~5-10s/array + some overhead)
  units <- fitCnProbes(plm, verbose=verbose);
  str(units);
  # int [1:2377527] 401698 401699 401700 401701 401702 401703 ...

  # Fit remaining units, i.e. SNPs (~5-10min/array)
  units <- fit(plm, verbose=verbose);
  str(units);
  # int [1:418181] 1 2 3 4 5 6 7 8 9 10 ...
}

ces <- getChipEffectSet(plm);
print(ces);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Export (theta, beta) for all arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsList <- exportTotalAndFracB(ces, verbose=verbose);
print(dsList);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# TCN Segmentation and plotting
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cns <- CbsModel(dsList$total);
print(cns);

ce <- ChromosomeExplorer(cns, zooms=2^(0:5));
print(ce);
process(ce, chromosomes=c(19, 22, 23), verbose=verbose);

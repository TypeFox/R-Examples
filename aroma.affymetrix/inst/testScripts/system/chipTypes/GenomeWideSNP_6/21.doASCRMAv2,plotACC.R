library("aroma.affymetrix")
verbose <- Arguments$getVerbose(-4, timestamp=TRUE);

dataSet <- "GSE13372,testset";
chipType <- "GenomeWideSNP_6,Full";
csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# AS-CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
res <- doASCRMAv2(csR, drop=FALSE, verbose=verbose);
print(res);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot allele pairs before and after calibration
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
acc <- res$acc;
for (what in c("input", "output")) {
  toPNG(getFullName(acc), tags=c("allelePairs", what), aspectRatio=0.7, {
    plotAllelePairs(acc, array=1, what=what, verbose=verbose);
  });
}

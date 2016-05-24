library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-4, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tests for setting up CEL sets and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "HapMap,CEU,testset";
chipType <- "Mapping50K_Hind240";

cdf <- AffymetrixCdfFile$byChipType(chipType);
csR <- AffymetrixCelSet$byName(dataSet, cdf=cdf);
setFullNamesTranslator(csR, function(names, ...) {
  # Turn into comma-separated tags
  names <- gsub("_", ",", names);
  # Drop any Hind/Xba tags
  names <- gsub(",(Hind|Xba)", "", names);
  names;
})
print(csR);
print(getFullNames(csR));

# Check SAF attributes
nXY <- t(sapply(csR, function(cf) getAttributes(cf)[c("n23", "n24")]));
rownames(nXY) <- getNames(csR);
print(verbose, nXY);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allelic cross-talk calibration tests
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
acc <- AllelicCrosstalkCalibration(csR);
print(acc);
csC <- process(acc, verbose=verbose);
print(csC);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Normalization for nucleotide-position probe sequence effects
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bpn <- BasePositionNormalization(csC, target="zero");
print(bpn);
csN <- process(bpn, verbose=verbose);
print(csN);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Probe-level modelling test (for total CN analysis)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plm <- RmaCnPlm(csN, mergeStrands=TRUE, combineAlleles=TRUE, shift=300);
print(plm);
fit(plm, verbose=verbose);
ces <- getChipEffectSet(plm);
print(ces);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Normalize for PCR fragment-length effects
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fln <- FragmentLengthNormalization(ces, target="zero");
print(fln);
cesN <- process(fln, verbose=verbose);
print(cesN);

# The attributes applies to chip effects as well
nXY0 <- nXY;
nXY <- t(sapply(cesN, function(cf) getAttributes(cf)[c("n23", "n24")]));
rownames(nXY) <- getNames(cesN);
print(nXY);
stopifnot(identical(nXY, nXY0));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calculate sex-chromosome bias-corrected reference signals
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# For Chr1-22 and ChrX the copy-neutral ploidy should be two.
# If the ploidy of a sample is unknown, assume the default is two.
ceR <- calculateBaseline(cesN, chromosomes=1:23, ploidy=2,
                                       defaultPloidy=2, verbose=verbose);
print(ceR);

# Compare to regular median average
ceRb <- getAverageFile(cesN, verbose=verbose);
print(ceRb);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Verify that the ChrX CNs are bias corrected
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cbs <- CbsModel(cesN, ceR);
print(cbs);

M <- NULL;
for (kk in seq_along(cbs)) {
  rawCNs <- extractRawCopyNumbers(cbs, array=kk, chromosome=23);
  rawCNs <- getSignals(rawCNs);
  M <- cbind(M, rawCNs);
}
colnames(M) <- getArrays(cbs);

n23 <- sapply(cesN, getAttribute, "n23");
col <- c("blue", "red")[n23];
Mlab <- expression(log[2](theta/theta[R]));
Mlim <- c(-5,2);

boxplot(M, col=col, ylim=Mlim, ylab=Mlab, las=2);
abline(h=0, lty=4);
title("Copy numbers on ChrX\n(bias corrected)");


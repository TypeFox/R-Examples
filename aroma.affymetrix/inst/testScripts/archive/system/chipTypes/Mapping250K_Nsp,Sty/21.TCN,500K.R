library("aroma.affymetrix");
log <- Arguments$getVerbose(-4, timestamp=TRUE);


dataSetName <- "Affymetrix_2006-TumorNormal";
chipTypes <- c("Mapping250K_Nsp", "Mapping250K_Sty");

pairs <- matrix(c(
  "CRL-2325D", "CRL-2324D",
  "CRL-5957D", "CRL-5868D",
  "CCL-256.1D", "CCL-256D",
  "CRL-2319D", "CRL-2320D",
  "CRL-2362D", "CRL-2321D",
  "CRL-2337D", "CRL-2336D",
  "CRL-2339D", "CRL-2338D",
  "CRL-2341D", "CRL-2340D",
  "CRL-2346D", "CRL-2314D"
), ncol=2, byrow=TRUE);
colnames(pairs) <- c("normal", "tumor");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tests for setting up CEL sets and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csRList <- list();
for (chipType in chipTypes) {
  cs <- AffymetrixCelSet$byName(dataSetName, chipType=chipType, verbose=log);
  print(cs);
  stopifnot(all(getNames(cs) %in% pairs));
  csRList[[chipType]] <- cs;
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allelic cross-talk calibration
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csList <- csRList;
csCList <- list();
for (chipType in names(csList)) {
  cs <- csList[[chipType]];
  acc <- AllelicCrosstalkCalibration(cs);
  print(acc);
  csC <- process(acc, verbose=log);
  print(csC);
  stopifnot(identical(getNames(csC), getNames(cs)));
  csCList[[chipType]] <- csC;
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Probe-level modelling test (for CN analysis)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csList <- csCList;
cesCnList <- list();
for (chipType in names(csList)) {
  cs <- csList[[chipType]];
  plm <- RmaCnPlm(cs, mergeStrands=TRUE, combineAlleles=TRUE, shift=300);
  print(plm);
  fit(plm, verbose=log);
  ces <- getChipEffectSet(plm);
  print(ces);
  stopifnot(identical(getNames(ces), getNames(cs)));
  cesCnList[[chipType]] <- ces;
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fragment-length normalization test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cesCnList <- cesCnList;
cesNList <- list();
for (chipType in names(csList)) {
  ces <- cesCnList[[chipType]];
  fln <- FragmentLengthNormalization(ces);
  print(fln);
  cesN <- process(fln, verbose=log);
  print(cesN);
  stopifnot(identical(getNames(cesN), getNames(ces)));
  cesNList[[chipType]] <- cesN;
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup a paired CBS model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Split data set in (tumor, normal) pairs
sets <- list(tumor=list(), normal=list());
for (chipType in names(cesNList)) {
  ces <- cesNList[[chipType]];
  for (type in colnames(pairs)) {
    idxs <- match(pairs[,type], getNames(ces));
    sets[[type]][[chipType]] <- ces[idxs];
  }
}
cns <- CbsModel(sets$tumor, sets$normal, maxNAFraction=1/5);
print(cns);

# Link the ChromosomeExplorer to the segmentation model
ce <- ChromosomeExplorer(cns);
print(ce);

# Fit the model for a few chromosomes
process(ce, arrays=1:2, chromosomes=c(2, 19), verbose=log);

# The X chromosome is very noisy and generates quite a few missing values
process(ce, arrays=1:2, chromosomes=23, verbose=log);

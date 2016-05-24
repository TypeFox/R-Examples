library("aroma.affymetrix");
log <- verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

dataSetName <- "Affymetrix_2006-TumorNormal";
chipType <- "Mapping250K_Nsp";

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
# Setting up data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- AffymetrixCdfFile$byChipType(chipType);
csR <- AffymetrixCelSet$byName(dataSetName, cdf=cdf);
print(csR);

# Reorder arrays according to 'pairs' matrix
csR <- csR[indexOf(csR, pairs)];

acc <- AllelicCrosstalkCalibration(csR, model="CRMAv2");
print(acc);

csC <- process(acc, verbose=log);
print(csC);

bpn <- BasePositionNormalization(csC, target="zero");
print(bpn);

csN <- process(bpn, verbose=log);
print(csN);

plm <- RmaCnPlm(csN, mergeStrands=TRUE, combineAlleles=FALSE);
print(plm);

fit(plm, verbose=log);

ces <- getChipEffectSet(plm);
print(ces);

fln <- FragmentLengthNormalization(ces, target="zero");
print(fln);

cesN <- process(fln, verbose=log);
print(cesN);

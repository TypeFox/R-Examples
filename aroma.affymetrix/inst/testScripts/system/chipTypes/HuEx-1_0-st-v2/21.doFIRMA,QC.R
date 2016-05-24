library("aroma.affymetrix");
library("matrixStats"); # rowMedians()
verbose <- Arguments$getVerbose(-3, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "Affymetrix-HeartBrain";
chipType <- "HuEx-1_0-st-v2";
cdf <- AffymetrixCdfFile$byChipType(chipType, tags="coreR3,A20071112,EP");
print(cdf);

# Setup CEL set using the core CDF.
csR <- AffymetrixCelSet$byName(dataSet, cdf=cdf);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Process only cerebellum and heart
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
types <- c("cerebellum", "heart");
csR <- csR[indexOf(csR, patterns=types)];

setFullName(csR, sprintf("%s,%s", dataSet, paste(types, collapse="+")));
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# FIRMA
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
res <- doFIRMA(csR, drop=FALSE, verbose=verbose);
print(res);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# PLMs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plmList <- list(
  merge   = res$plm,                               # all exons together
  noMerge = ExonRmaPlm(res$csN, mergeGroups=FALSE) # each exon separately
);
print(plmList);

# Fit the per-exon PLM
fit(plmList$noMerge, verbose=verbose);

# Chip effects
cesList <- lapply(plmList, FUN=getChipEffectSet);
print(cesList);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# RLE/NUSE for first 100 units
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# The reason for observing a small difference is scores below,
# is that when saving to file the estimates are rounded of to
# floats, whereas the estimates calculated in memory are kept
# in full precision.
tol <- 1e-5;

units <- 1:100;

for (key in names(cesList)) {
  verbose && enter(verbose, key);

  ces <- cesList[[key]];

  # Assert correctness of unit-specific RLE scores
  theta <- extractMatrix(ces, field="theta", units=units);
  thetaR <- 2^rowMedians(log2(theta), na.rm=TRUE);
  rle0 <- log2(theta/thetaR);
  rle1 <- extractMatrix(ces, field="RLE", units=units, verbose=verbose);
  stopifnot(all.equal(rle1, rle0, tolerance=tol));

  # Assert correctness of boxplot statistics of RLE scores
  stats0 <- boxplot.stats(rle0[,1]);
  stats1 <- boxplotStats(ces, type="RLE", arrays=1:2, subset=units);
  stopifnot(all.equal(stats1[[1]], stats0, tolerance=tol));

  # Assert correctness of unit-specific NUSE scores
  se <- extractMatrix(ces, field="sdTheta", units=units);
  seR <- 2^rowMedians(log2(se), na.rm=TRUE);
  nuse0 <- log2(se)/log2(seR);
  nuse1 <- extractMatrix(ces, field="NUSE", units=units, verbose=verbose);
  stopifnot(all.equal(nuse1, nuse0, tolerance=tol));

  # Assert correctness of boxplot statistics of NUSE scores
  stats0 <- boxplot.stats(nuse0[,1]);
  stats1 <- boxplotStats(ces, type="NUSE", arrays=1:2, subset=units);
  stopifnot(all.equal(stats1[[1]], stats0, tolerance=tol));

  verbose && exit(verbose);
} # for (key in ...)

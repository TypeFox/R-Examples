library("aroma.affymetrix");

log <- Arguments$getVerbose(-4, timestamp=TRUE);

dataSet <- "HapMap,CEU,testset";
chipType <- "Mapping50K_Hind240";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Assert existence of probe-sequence annotation files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
acs <- AromaCellSequenceFile$byChipType(chipType);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tests for setting up CEL sets and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType);
print(csR);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Allelic cross-talk calibration tests
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
acc <- AllelicCrosstalkCalibration(csR);
print(acc);
csC <- process(acc, verbose=log);
print(csC);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Base-position normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bpn <- BasePositionNormalization(csC, shift=+300);
print(bpn);

csN <- process(bpn, verbose=log);
print(csN);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot probe sequence base count effects
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ff <- 3;

counts <- countBases(acs, verbose=log);

Ms <- list();
whats <- c("raw", "acc", "bcn");
for (what in whats) {
  if (what == "raw") {
    cs <- csR;
  } else if (what == "acc") {
    cs <- csC;
  } else if (what == "bcn") {
    cs <- csN;
  }
  cf <- cs[[ff]];
  cfR <- getAverageFile(cs, verbose=log);
  yR <- readRawData(cfR, fields="intensities", drop=TRUE);
  y <- readRawData(cf, fields="intensities", drop=TRUE);
  M <- log2(y/yR);
  Ms[[what]] <- M;
  rm(y, yR, cfR, cs, M);
} # for (what ...)


# Estimate standard deviations for log-ratios
print(sapply(Ms, FUN=function(M) { mad(diff(M), na.rm=TRUE) }));

Mlim <- c(-1,1);
Mlab <- expression(M == log[2](y/y[R]));
xlim <- c(0, 25);

layout(matrix(1:(length(whats)*4), ncol=4, byrow=TRUE));
par(mar=c(5,4,2,1)+0.1);
for (what in names(Ms)) {
  M <- Ms[[what]];

  for (bb in c("A", "C", "G", "T")) {
    xlab <- sprintf("Number of %s:s", bb);
    boxplot(M ~ counts[,bb], outline=FALSE,
            ylim=Mlim, xlim=xlim, ylab=Mlab, xlab=xlab);
    stext(side=3, pos=0, getFullName(cf));
    stext(side=3, pos=0.98, line=-1, cex=2, bb);
    stextChipType(chipType);
  }
  rm(M);
} # for (what ...)
devDone();

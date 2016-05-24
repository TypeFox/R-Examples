library("PSCBS")
subplots <- R.utils::subplots
stext <- R.utils::stext

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load SNP microarray data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
data <- PSCBS::exampleData("paired.chr01")
str(data)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Paired PSCBS segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Drop single-locus outliers
dataS <- dropSegmentationOutliers(data)

# Run light-weight tests by default
if (Sys.getenv("_R_CHECK_FULL_") == "") {
  # Use only every 5th data point
  dataS <- dataS[seq(from=1, to=nrow(data), by=5),]
  # Number of segments (for assertion)
  nSegs <- 3L
  # Number of bootstrap samples (see below)
  B <- 100L
} else {
  # Full tests
  nSegs <- 12L
  B <- 1000L
}

str(dataS)

R.oo::attachLocally(dataS)


gaps <- findLargeGaps(dataS, minLength=2e6);
knownSegments <- gapsToSegments(gaps, dropGaps=TRUE);

# Paired PSCBS segmentation
fit <- segmentByPairedPSCBS(dataS, knownSegments=knownSegments,
                            seed=0xBEEF, verbose=-10);
print(fit);

fit1 <- fit;
fit2 <- renameChromosomes(fit1, from=1, to=2);
fit <- append(fit1, fit2);
knownSegments <- tileChromosomes(fit)$params$knownSegments;

segList <- seqOfSegmentsByDP(fit, verbose=-10);
K <- length(segList);
ks <- seq(from=1, to=K, length.out=min(5,K));
subplots(length(ks), ncol=1, byrow=TRUE);
par(mar=c(2,1,1,1));
for (kk in ks) {
  knownSegmentsKK <- segList[[kk]];
  fitKK <- resegment(fit, knownSegments=knownSegmentsKK, undoTCN=+Inf, undoDH=+Inf);
  plotTracks(fitKK, tracks="tcn,c1,c2", Clim=c(0,5), add=TRUE);
  abline(v=c(knownSegments$start, knownSegments$end)/1e6, lty=3);
  stext(side=3, pos=0, sprintf("Number of segments: %d", nrow(knownSegmentsKK)));
} # for (kk ...)

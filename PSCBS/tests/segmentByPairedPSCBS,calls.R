library("PSCBS")

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

# Find centromere
gaps <- findLargeGaps(dataS, minLength=2e6)
knownSegments <- gapsToSegments(gaps)


# Run light-weight tests by default
if (Sys.getenv("_R_CHECK_FULL_") == "") {
  # Use only every 5th data point
  dataS <- dataS[seq(from=1, to=nrow(data), by=5),]
  # Number of segments (for assertion)
  nSegs <- 4L
  # Number of bootstrap samples (see below)
  B <- 100L
} else {
  # Full tests
  nSegs <- 11L
  B <- 1000L
}

str(dataS)

# Paired PSCBS segmentation
fit <- segmentByPairedPSCBS(dataS, knownSegments=knownSegments,
                            seed=0xBEEF, verbose=-10)
print(fit)

# Plot results
plotTracks(fit)

# Sanity check
stopifnot(nbrOfSegments(fit) == nSegs)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Bootstrap segment level estimates
# (used by the AB caller, which, if skipped here,
#  will do it automatically)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit <- bootstrapTCNandDHByRegion(fit, B=B, verbose=-10)
print(fit)
plotTracks(fit)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calling segments with run of homozygosity (ROH)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit <- callROH(fit, verbose=-10)
print(fit)
plotTracks(fit)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Estimate background
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
kappa <- estimateKappa(fit, verbose=-10)
print(kappa)
## [1] 0.226011


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calling segments in allelic balance (AB)
# NOTE: Ideally, this should be done on whole-genome data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Explicitly estimate the threshold in DH for calling AB
# (which be done by default by the caller, if skipped here)
deltaAB <- estimateDeltaAB(fit, flavor="qq(DH)", verbose=-10)
if (Sys.getenv("_R_CHECK_FULL_") == "") {
  # Ad hoc workaround for not utilizing all of the data
  # in the test, which results in a poor estimate
  deltaAB <- 0.165
}
print(deltaAB)
## [1] 0.1657131

fit <- callAB(fit, delta=deltaAB, verbose=-10)
print(fit)
plotTracks(fit)

# Even if not explicitly specified, the estimated
# threshold parameter is returned by the caller
stopifnot(fit$params$deltaAB == deltaAB)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calling segments in loss-of-heterozygosity (LOH)
# NOTE: Ideally, this should be done on whole-genome data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Explicitly estimate the threshold in C1 for calling LOH
# (which be done by default by the caller, if skipped here)
deltaLOH <- estimateDeltaLOH(fit, flavor="minC1|nonAB", verbose=-10)
print(deltaLOH)
## [1] 0.625175

fit <- callLOH(fit, delta=deltaLOH, verbose=-10)
print(fit)
plotTracks(fit)

# Even if not explicitly specified, the estimated
# threshold parameter is returned by the caller
stopifnot(fit$params$deltaLOH == deltaLOH)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calling segments that are gained, copy neutral, and lost
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit <- callGNL(fit, verbose=-10)
print(fit)
plotTracks(fit)

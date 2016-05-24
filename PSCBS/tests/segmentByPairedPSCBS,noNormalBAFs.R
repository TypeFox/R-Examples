library("PSCBS")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load SNP microarray data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
data <- PSCBS::exampleData("paired.chr01")
str(data)

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
  nSegs <- 8L
  B <- 1000L
}

str(dataS)

R.oo::attachLocally(dataS)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulate that genotypes are known by other means
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
library("aroma.light")
muN <- aroma.light::callNaiveGenotypes(betaN, censorAt=c(0,1))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Paired PSCBS segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit <- segmentByPairedPSCBS(CT, betaT=betaT, muN=muN, tbn=FALSE,
                            chromosome=chromosome, x=x,
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
# Calling segments in allelic balance (AB) and
# in loss-of-heterozygosity (LOH)
# NOTE: Ideally, this should be done on whole-genome data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit <- callAB(fit, verbose=-10)
fit <- callLOH(fit, verbose=-10)
print(fit)
plotTracks(fit)

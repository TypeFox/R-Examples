library("PSCBS")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load SNP microarray data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
data <- PSCBS::exampleData("paired.chr01")
str(data)

# Drop single-locus outliers
dataS <- dropSegmentationOutliers(data)

# Run light-weight tests
# Use only every 5th data point
dataS <- dataS[seq(from=1, to=nrow(data), by=5),]
# Number of segments (for assertion)
nSegs <- 3L
# Number of bootstrap samples (see below)
B <- 100L

str(dataS)
R.oo::attachLocally(dataS)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calculate DH
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
muN <- aroma.light::callNaiveGenotypes(betaN, censorAt=c(0,1))
# SNPs are identifies as those loci that have non-missing 'betaT' & 'muN'
isSnp <- (!is.na(betaT) & !is.na(muN))
isHet <- isSnp & (muN == 1/2)
rho <- rep(NA_real_, length=length(muN))
rho[isHet] <- 2*abs(betaT[isHet]-1/2)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Paired PSCBS segmentation using TCN and DH only
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit <- segmentByPairedPSCBS(CT, rho=rho,
                            chromosome=chromosome, x=x,
                            seed=0xBEEF, verbose=-10)
print(fit)

# Plot results
plotTracks(fit)


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

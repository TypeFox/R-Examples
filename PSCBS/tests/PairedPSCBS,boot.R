###########################################################
# This tests:
# - Bootstrapping for PairedPSCBS objects
###########################################################
library("PSCBS")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load SNP microarray data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
data <- PSCBS::exampleData("paired.chr01")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Paired PSCBS segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Drop single-locus outliers
dataS <- dropSegmentationOutliers(data)
dataS <- dataS[seq(from=1, to=nrow(data), by=5),]
nSegs <- 4L
str(dataS)
# Segment known regions
knownSegments <- data.frame(
  chromosome = c(        1,  1,         1),
  start      = c(     -Inf, NA, 141510003),
  end        = c(120992603, NA,      +Inf)
)
fit <- segmentByPairedPSCBS(dataS, knownSegments=knownSegments, avgDH="median", seed=0xBEEF)
print(fit)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Bootstrap
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
B <- 1L
seed <- 0xBEEF
probs <- c(0.025, 0.05, 0.95, 0.975)

sets <- getBootstrapLocusSets(fit, B=B, seed=seed)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Subset by first segment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ss <- 1L

fitT <- extractSegment(fit, ss)
dataT <- getLocusData(fitT)
segsT <- getSegments(fitT)

# Truth
bootT <- bootstrapSegmentsAndChangepoints(fitT, B=B, seed=seed)
bootT1 <- bootT$segments[1L,,,drop=FALSE]
types <- dimnames(bootT1)[[3L]]
dim(bootT1) <- dim(bootT1)[-1L]
colnames(bootT1) <- types
sumsT <- apply(bootT1, MARGIN=2L, FUN=quantile, probs=probs)
print(sumsT)

fitTB <- bootstrapTCNandDHByRegion(fitT, B=B, seed=seed)
segsTB <- getSegments(fitTB)
segsTB <- unlist(segsTB[,grep("_", colnames(segsTB))])
dim(segsTB) <- dim(sumsT)
dimnames(segsTB) <- dimnames(sumsT)
print(segsTB)

# Sanity check
stopifnot(all.equal(segsTB, sumsT))

# Calculate summaries using the existing bootstrap samples
fitTBp <- bootstrapTCNandDHByRegion(fitT, .boot=bootT)
# Sanity check
all.equal(fitTBp, fitTB)


# Bootstrap from scratch
setsT <- getBootstrapLocusSets(fitT, B=B, seed=seed)
lociT <- setsT$locusSet[[1L]]$bootstrap$loci
idxs <- lociT$tcn
tcnT <- array(dataT$CT[idxs], dim=dim(idxs))
tcnT <- apply(tcnT, MARGIN=2L, FUN=mean, na.rm=TRUE)
idxs <- lociT$dh
dhT <- array(dataT$rho[idxs], dim=dim(idxs))
dhT <- apply(dhT, MARGIN=2L, FUN=median, na.rm=TRUE)
c1T <- (1-dhT) * tcnT / 2
c2T <- tcnT - c1T
bootT2 <- array(c(tcnT, dhT, c1T, c2T), dim=c(1L, 4L))
colnames(bootT2) <- colnames(bootT1)
print(bootT2)

# This comparison is only valid if B == 1L
if (B == 1L) {
  # Sanity check
  stopifnot(all.equal(bootT2, bootT1))
}

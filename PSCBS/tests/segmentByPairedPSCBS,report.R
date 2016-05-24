# This test script calls a report generator which requires
# the 'ggplot2' package, which in turn will require packages
# 'colorspace', 'dichromat', 'munsell', 'reshape2' and 'scales'.

# Only run this test in full testing mode
if (Sys.getenv("_R_CHECK_FULL_") != "") {
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

# Speed up example by segmenting fewer loci
dataS <- dataS[seq(from=1, to=nrow(data), by=5),]

str(dataS)

gaps <- findLargeGaps(dataS, minLength=2e6)
knownSegments <- gapsToSegments(gaps)

# Paired PSCBS segmentation
fit <- segmentByPairedPSCBS(dataS, knownSegments=knownSegments,
                            seed=0xBEEF, verbose=-10)

# Fake a multi-chromosome segmentation
fit1 <- fit
fit2 <- renameChromosomes(fit, from=1, to=2)
fit <- append(fit1, fit2)

report(fit, sampleName="PairedPSCBS", studyName="PSCBS-Ex", verbose=-10)

} # if (Sys.getenv("_R_CHECK_FULL_"))

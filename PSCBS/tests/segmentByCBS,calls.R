# This test script calls a report generator which requires
# the 'ggplot2' package, which in turn will require packages
# 'colorspace', 'dichromat', 'munsell', 'reshape2' and 'scales'.

# Only run this test in full testing mode
if (Sys.getenv("_R_CHECK_FULL_") != "") {
library("PSCBS")
stext <- R.utils::stext

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load SNP microarray data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
data <- PSCBS::exampleData("paired.chr01")
str(data)

data <- data.frame(chromosome=data$chromosome, x=data$x, y=data$CT)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CBS segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Drop single-locus outliers
dataS <- dropSegmentationOutliers(data)

# Speed up example by segmenting fewer loci
dataS <- dataS[seq(from=1, to=nrow(data), by=5),]

str(dataS)

gaps <- findLargeGaps(dataS, minLength=2e6)
knownSegments <- gapsToSegments(gaps)

# CBS segmentation
fit <- segmentByCBS(dataS, knownSegments=knownSegments,
                            seed=0xBEEF, verbose=-10)
signalType(fit) <- "ratio"
plotTracks(fit)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Call using the UCSF MAD caller (not recommended)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fitC <- callGainsAndLosses(fit)
plotTracks(fitC)
pars <- fitC$params$callGainsAndLosses
stext(side=3, pos=1/2, line=-1, substitute(sigma==x, list(x=sprintf("%.2f", pars$sigmaMAD))))
mu <- pars$muR
tau <- unlist(pars[c("tauLoss", "tauGain")], use.names=FALSE)
abline(h=mu, lty=2, lwd=2)
abline(h=tau, lwd=2)
mtext(side=4, at=tau[1], expression(Delta[LOSS]), adj=-0.2, cex=0.7, las=2, xpd=TRUE)
mtext(side=4, at=tau[2], expression(Delta[GAIN]), adj=-0.2, cex=0.7, las=2, xpd=TRUE)
title(main="CN caller: \"ucsf-mad\"")


# Caller to be implemented
deltaCN <- estimateDeltaCN(fit)
tau <- mu + 1/2*c(-1,+1)*deltaCN
abline(h=tau, lty=2, lwd=1, col="red")



} # if (Sys.getenv("_R_CHECK_FULL_"))

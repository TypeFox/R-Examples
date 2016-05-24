###########################################################################
## This 'tangle' R script was created from an RSP document.
## RSP source document: './PairedPSCBS.tex.rsp'
## Metadata 'title': 'Parent-specific copy-number segmentation using Paired PSCBS'
## Metadata 'author': 'Henrik Bengtsson'
## Metadata 'engine': 'R.rsp::rsp'
## Metadata 'keywords': 'copy numbers, allele specific, parent specific, genomic aberrations'
###########################################################################

t0 <- Sys.time()
R.utils::use("R.utils")

# RSP specific
R.rsp <- R.oo::Package("R.rsp")
withCapture <- R.utils::withCapture
options("withCapture/newline"=FALSE)
options(str=strOptions(strict.width="cut"))
options(width=85)
options(digits=3)

# Graphics
use("R.devices")
options("devEval/args/field"="fullname") # Preferred for LaTeX
devOptions("png", width=840)

# Analysis
use("PSCBS")
PSCBS <- R.oo::Package("PSCBS")
fixLocations <- function(fit, ...) {
  for (key in grep("(end|start)$", colnames(fit$output))) {
    fit$output[[key]] <- as.integer(fit$output[[key]])
  }
  fit
} # fixLocations()
R.rsp$version
R.rsp$author
format(as.Date(PSCBS$date), format="%B %d, %Y")
fullname <- "PairedPSCBS,exData,chr01"
withCapture({
data <- PSCBS::exampleData("paired.chr01")
str(data)
})
withCapture({
data <- dropSegmentationOutliers(data)
})
withCapture({
gaps <- findLargeGaps(data, minLength=1e6)
gaps
})
withCapture({
knownSegments <- gapsToSegments(gaps)
knownSegments
})
withCapture({
fit <- segmentByPairedPSCBS(data, knownSegments=knownSegments, preserveScale=FALSE, seed=0xBEEF, verbose=-10)
})
nbrOfSegments(fit)
fit <- fixLocations(fit)
withCapture({
getSegments(fit, simplify=TRUE)
})
segs <- getSegments(fit, simplify=TRUE)
which(segs$tcnNbrOfLoci == 0)
which(segs$dhNbrOfLoci == 0)
toPNG(fullname, tags=c(class(fit)[1L], "tracks"), aspectRatio=0.6, {
    plotTracks(fit)
  })
withCapture({
fit <- callROH(fit, verbose=-10)
})
withCapture({
fit <- callAB(fit, verbose=-10)
})
withCapture({
fit <- callLOH(fit, verbose=-10)
})
withCapture({
fit <- callNTCN(fit, verbose=-10)
})
withCapture({
getSegments(fit, simplify=TRUE)
})
segs <- getSegments(fit, simplify=TRUE)
withCapture({
fitP <- pruneByHClust(fit, h=0.25, verbose=-10)
})
toPNG(fullname, tags=c(class(fitP)[1L], "pruned", "tracks"), aspectRatio=0.6, {
    plotTracks(fitP)
  })
toLatex(sessionInfo())
dt <- round(Sys.time()-t0, digits=2)
attr(dt, "units")

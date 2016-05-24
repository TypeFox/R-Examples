library("PSCBS")
subplots <- R.utils::subplots

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulating copy-number data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
set.seed(0xBEEF)

# Number of loci
J <- 1000

mu <- double(J)
eps <- rnorm(J, sd=1/2)
y <- mu + eps
x <- sort(runif(length(y), max=length(y)))

idxs <- which(200 <= x & x < 300)
y[idxs] <- y[idxs] + 1
idxs <- which(350 <= x & x < 400)
y[idxs] <- NA # centromere
x[idxs] <- NA # centromere
idxs <- which(650 <= x & x < 800)
y[idxs] <- y[idxs] - 1
x <- x*1e5

keep <- is.finite(x)
x <- x[keep]
y <- y[keep]

data <- list()
for (chr in 1:2) {
  data[[chr]] <- data.frame(chromosome=chr, y=y, x=x)
}
data <- Reduce(rbind, data)


subplots(7, ncol=1)
par(mar=c(1.7,1,0.2,1)+0.1)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit <- segmentByCBS(data)
print(fit)

Clim <- c(-3,3) + c(0,10)
plotTracks(fit, Clim=Clim)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Shifting every other chromosome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fitList <- list()
chrs <- getChromosomes(fit)
for (kk in seq(along=chrs)) {
  chr <- chrs[kk]
  fitKK <- extractChromosome(fit, chr)
  if (kk %% 2 == 0) {
    fitKK <- shiftTCN(fitKK, shift=+10)
  }
  fitList[[kk]] <- fitKK
} # for (kk ...)
fitT <- Reduce(append, fitList)
# Sanity check
stopifnot(nbrOfSegments(fitT) == nbrOfSegments(fit))

plotTracks(fitT, Clim=Clim)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Shifting every other known segment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
gaps <- findLargeGaps(data, minLength=40e5)
knownSegments <- gapsToSegments(gaps, dropGaps=TRUE)
fit <- segmentByCBS(data, knownSegments=knownSegments)

subplots(2, ncol=1)
plotTracks(fit, Clim=Clim)
abline(v=c(knownSegments$start, knownSegments$end)/1e6, lty=3)

fitList <- list()
for (kk in seq(length=nrow(knownSegments))) {
  seg <- knownSegments[kk,]
  start <- seg$start
  end <- seg$end
  fitKK <- extractChromosome(fit, seg$chromosome)
  segsKK <- getSegments(fitKK)
  idxStart <- min(which(segsKK$start >= start))
  idxEnd <- max(which(segsKK$end <= end))
  idxs <- idxStart:idxEnd
  fitKK <- extractSegments(fitKK, idxs)
  if (kk %% 2 == 0) {
    fitKK <- shiftTCN(fitKK, shift=+10)
  }
  fitList[[kk]] <- fitKK
} # for (kk ...)
fitT <- Reduce(append, fitList)
# Sanity check
stopifnot(nbrOfSegments(fitT) == nbrOfSegments(fit))

plotTracks(fitT, Clim=Clim)
abline(v=c(knownSegments$start, knownSegments$end)/1e6, lty=3)


segList <- seqOfSegmentsByDP(fit);
K <- length(segList);
subplots(K, ncol=2, byrow=FALSE);
par(mar=c(2,1,1,1));
for (kk in 1:K) {
  knownSegments <- segList[[kk]];
  fitKK <- resegment(fit, knownSegments=knownSegments, undo=+Inf);
  plotTracks(fitKK, Clim=c(-3,3));
} # for (kk ...)

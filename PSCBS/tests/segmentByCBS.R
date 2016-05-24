###########################################################
# This tests:
# - segmentByCBS(...)
# - segmentByCBS(..., knownSegments)
# - tileChromosomes()
# - plotTracks()
###########################################################
library("PSCBS")
subplots <- R.utils::subplots

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulating copy-number data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
set.seed(0xBEEF)

# Number of loci
J <- 1000

mu <- double(J)
mu[200:300] <- mu[200:300] + 1
mu[350:400] <- NA # centromere
mu[650:800] <- mu[650:800] - 1
eps <- rnorm(J, sd=1/2)
y <- mu + eps
x <- sort(runif(length(y), max=length(y))) * 1e5
w <- runif(J)
w[650:800] <- 0.001


subplots(8, ncol=1L)
par(mar=c(1.7,1,0.2,1)+0.1)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit <- segmentByCBS(y, x=x)
sampleName(fit) <- "CBS_Example"
print(fit)
plotTracks(fit)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segmentation with some known change points
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
knownSegments <- data.frame(
  chromosome=c(    0,   0),
  start     =x[c(  1, 401)],
  end       =x[c(349,   J)]
)
fit2 <- segmentByCBS(y, x=x, knownSegments=knownSegments, verbose=TRUE)
sampleName(fit2) <- "CBS_Example_2"
print(fit2)
plotTracks(fit2)
abline(v=c(knownSegments$start, knownSegments$end)/1e6, lty=3)


# Chromosome boundaries can be specified as -Inf and +Inf
knownSegments <- data.frame(
  chromosome=c(     0,      0),
  start     =c(  -Inf, x[401]),
  end       =c(x[349],   +Inf)
)
fit2b <- segmentByCBS(y, x=x, knownSegments=knownSegments, verbose=TRUE)
sampleName(fit2b) <- "CBS_Example_2b"
print(fit2b)
plotTracks(fit2b)
abline(v=c(knownSegments$start, knownSegments$end)/1e6, lty=3)


# As a proof of concept, it is possible to segment just the centromere,
# which contains no data.  All statistics will be NAs.
knownSegments <- data.frame(
  chromosome=c(    0),
  start     =x[c(350)],
  end       =x[c(400)]
)
fit3 <- segmentByCBS(y, x=x, knownSegments=knownSegments, verbose=TRUE)
sampleName(fit3) <- "CBS_Example_3"
print(fit3)
plotTracks(fit3, Clim=c(0,5), xlim=c(0,100))
abline(v=c(knownSegments$start, knownSegments$end)/1e6, lty=3)



# If one specify the (empty) centromere as a segment, then its
# estimated statistics will be NAs, which becomes a natural
# separator between the two "independent" arms.
knownSegments <- data.frame(
  chromosome=c(    0,   0,   0),
  start     =x[c(  1, 350, 401)],
  end       =x[c(349, 400,   J)]
)
fit4 <- segmentByCBS(y, x=x, knownSegments=knownSegments, verbose=TRUE)
sampleName(fit4) <- "CBS_Example_4"
print(fit4)
plotTracks(fit4)
abline(v=c(knownSegments$start, knownSegments$end)/1e6, lty=3)



fit5 <- segmentByCBS(y, x=x, knownSegments=knownSegments, undo=Inf, verbose=TRUE)
sampleName(fit5) <- "CBS_Example_5"
print(fit5)
plotTracks(fit5)
abline(v=c(knownSegments$start, knownSegments$end)/1e6, lty=3)
stopifnot(nbrOfSegments(fit5) == nrow(knownSegments))


# One can also force a separator between two segments by setting
# 'start' and 'end' to NAs ('chromosome' has to be given)
knownSegments <- data.frame(
  chromosome=c(    0,  0,   0),
  start     =x[c(  1, NA, 401)],
  end       =x[c(349, NA,   J)]
)
fit6 <- segmentByCBS(y, x=x, knownSegments=knownSegments, verbose=TRUE)
sampleName(fit6) <- "CBS_Example_6"
print(fit6)
plotTracks(fit6)
abline(v=c(knownSegments$start, knownSegments$end)/1e6, lty=3)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segment multiple chromosomes
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulate multiple chromosomes
fit1 <- renameChromosomes(fit, from=0, to=1)
fit2 <- renameChromosomes(fit, from=0, to=2)
fitM <- append(fit1, fit2)
fitM <- segmentByCBS(fitM)
sampleName(fitM) <- "CBS_Example_M"
print(fitM)
plotTracks(fitM, Clim=c(-3,3))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tiling multiple chromosomes
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tile chromosomes
fitT <- tileChromosomes(fitM)
fitTb <- tileChromosomes(fitT)
stopifnot(identical(fitTb, fitT))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Write segmentation to file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pathT <- tempdir()

## Tab-delimited file
pathname <- writeSegments(fitM, path=pathT)
print(pathname)

## WIG file
pathname <- writeWIG(fitM, path=pathT)
print(pathname)

unlink(pathT, recursive=TRUE)

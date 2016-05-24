library("PSCBS")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulating copy-number data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
set.seed(0xBEEF)

# Number of loci
J <- 1000

x <- sort(runif(J, max=J)) * 1e5

mu <- double(J)
mu[200:300] <- mu[200:300] + 1
mu[350:400] <- NA # centromere
mu[650:800] <- mu[650:800] - 1
eps <- rnorm(J, sd=1/2)
y <- mu + eps

outliers <- seq(from=1L, to=J, length.out=0.2*J)
y[outliers] <- y[outliers] + 1.5

w <- rep(1.0, times=J)
w[outliers] <- 0.01

data <- data.frame(chromosome=1L, x=x, y=y)
dataW <- cbind(data, w=w)


par(mar=c(2,3,0.2,1)+0.1)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Single-chromosome segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segment without weights
fit <- segmentByCBS(data)
sampleName(fit) <- "CBS_Example"
print(fit)
plotTracks(fit)
## Highlight outliers (they pull up the mean levels)
points(x[outliers]/1e6, y[outliers], col="purple")

# Segment with weights
fitW <- segmentByCBS(dataW)
sampleName(fitW) <- "CBS_Example (weighted)"
print(fitW)
drawLevels(fitW, col="red")

legend("topright", bg="white", legend=c("outliers", "non-weighted CBS", "weighted CBS"), col=c("purple", "purple", "red"), lwd=c(NA,3,3), pch=c(1,NA,NA))

## Assert that weighted segment means are less biased
dmean <- getSegments(fit)$mean - getSegments(fitW)$mean
cat("Segment mean differences:\n")
print(dmean)
stopifnot(all(dmean > 0, na.rm=TRUE))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segmentation with some known change points
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
knownSegments <- data.frame(
  chromosome=c(    1,   1),
  start     =x[c(  1, 401)],
  end       =x[c(349,   J)]
)
fit2 <- segmentByCBS(dataW, knownSegments=knownSegments, verbose=TRUE)
sampleName(fit2) <- "CBS_Example_2 (weighted)"
print(fit2)
plotTracks(fit2)
abline(v=c(knownSegments$start, knownSegments$end)/1e6, lty=3)


# Chromosome boundaries can be specified as -Inf and +Inf
knownSegments <- data.frame(
  chromosome=c(     1,      1),
  start     =c(  -Inf, x[401]),
  end       =c(x[349],   +Inf)
)
fit2b <- segmentByCBS(dataW, knownSegments=knownSegments, verbose=TRUE)
sampleName(fit2b) <- "CBS_Example_2b (weighted)"
print(fit2b)
plotTracks(fit2b)
abline(v=c(knownSegments$start, knownSegments$end)/1e6, lty=3)


# As a proof of concept, it is possible to segment just the centromere,
# which contains no data.  All statistics will be NAs.
knownSegments <- data.frame(
  chromosome=c(    1),
  start     =x[c(350)],
  end       =x[c(400)]
)
fit3 <- segmentByCBS(dataW, knownSegments=knownSegments, verbose=TRUE)
sampleName(fit3) <- "CBS_Example_3"
print(fit3)
plotTracks(fit3, Clim=c(0,5), xlim=c(0,100))
abline(v=c(knownSegments$start, knownSegments$end)/1e6, lty=3)


# If one specify the (empty) centromere as a segment, then its
# estimated statistics will be NAs, which becomes a natural
# separator between the two "independent" arms.
knownSegments <- data.frame(
  chromosome=c(    1,   1,   1),
  start     =x[c(  1, 350, 401)],
  end       =x[c(349, 400,   J)]
)
fit4 <- segmentByCBS(dataW, knownSegments=knownSegments, verbose=TRUE)
sampleName(fit4) <- "CBS_Example_4"
print(fit4)
plotTracks(fit4)
abline(v=c(knownSegments$start, knownSegments$end)/1e6, lty=3)


fit5 <- segmentByCBS(dataW, knownSegments=knownSegments, undo=Inf, verbose=TRUE)
sampleName(fit5) <- "CBS_Example_5"
print(fit5)
plotTracks(fit5)
abline(v=c(knownSegments$start, knownSegments$end)/1e6, lty=3)
stopifnot(nbrOfSegments(fit5) == nrow(knownSegments))


# One can also force a separator between two segments by setting
# 'start' and 'end' to NAs ('chromosome' has to be given)
knownSegments <- data.frame(
  chromosome=c(    1,  1,   1),
  start     =x[c(  1, NA, 401)],
  end       =x[c(349, NA,   J)]
)
fit6 <- segmentByCBS(dataW, knownSegments=knownSegments, verbose=TRUE)
sampleName(fit6) <- "CBS_Example_6"
print(fit6)
plotTracks(fit6)
abline(v=c(knownSegments$start, knownSegments$end)/1e6, lty=3)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Multi-chromosome segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
data2 <- data
data2$chromosome <- 2L
data <- rbind(data, data2)
dataW <- cbind(data, w=w)

par(mar=c(2,3,0.2,1)+0.1)
# Segment without weights
fit <- segmentByCBS(data)
sampleName(fit) <- "CBS_Example"
print(fit)
plotTracks(fit, Clim=c(-3,3))

# Segment with weights
fitW <- segmentByCBS(dataW)
sampleName(fitW) <- "CBS_Example (weighted)"
print(fitW)
drawLevels(fitW, col="red")

legend("topright", bg="white", legend=c("outliers", "non-weighted CBS", "weighted CBS"), col=c("purple", "purple", "red"), lwd=c(NA,3,3), pch=c(1,NA,NA))

## Assert that weighted segment means are less biased
dmean <- getSegments(fit)$mean - getSegments(fitW)$mean
cat("Segment mean differences:\n")
print(dmean)
stopifnot(all(dmean > 0, na.rm=TRUE))

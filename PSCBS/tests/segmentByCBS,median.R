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


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Single-chromosome segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
par(mar=c(2,3,0.2,1)+0.1)
# Segment without weights
fit <- segmentByCBS(data)
sampleName(fit) <- "CBS_Example"
print(fit)
plotTracks(fit)
## Highlight outliers (they pull up the mean levels)
points(x[outliers]/1e6, y[outliers], col="purple")

# Segment without weights but with median
fitM <- segmentByCBS(data, avg="median")
sampleName(fitM) <- "CBS_Example (median)"
print(fitM)
drawLevels(fitM, col="magenta", lty=3)

# Segment with weights
fitW <- segmentByCBS(dataW, avg="median")
sampleName(fitW) <- "CBS_Example (weighted)"
print(fitW)
drawLevels(fitW, col="red")

# Segment with weights and median
fitWM <- segmentByCBS(dataW, avg="median")
sampleName(fitWM) <- "CBS_Example (weighted median)"
print(fitWM)
drawLevels(fitWM, col="orange", lty=3)

legend("topright", bg="white", legend=c("outliers", "non-weighted CBS (mean)", "non-weighted CBS (median)", "weighted CBS (mean)", "weighted CBS (median)"), col=c("purple", "purple", "magenta", "red", "orange"), lwd=c(NA,3,3,3,3), lty=c(NA,1,3,1,3), pch=c(1,NA,NA,NA,NA))

## Assert that weighted segment means are less biased
dmean <- getSegments(fit)$mean - getSegments(fitW)$mean
cat("Segment mean differences:\n")
print(dmean)
stopifnot(all(dmean > 0, na.rm=TRUE))

dmean <- getSegments(fitM)$mean - getSegments(fitWM)$mean
cat("Segment median differences:\n")
print(dmean)
stopifnot(all(dmean > 0, na.rm=TRUE))


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

# Segment without weights but with median
fitM <- segmentByCBS(data, avg="median")
sampleName(fitM) <- "CBS_Example (median)"
print(fitM)
drawLevels(fitM, col="magenta", lty=3)

# Segment with weights
fitW <- segmentByCBS(dataW, avg="median")
sampleName(fitW) <- "CBS_Example (weighted)"
print(fitW)
drawLevels(fitW, col="red")

# Segment with weights and median
fitWM <- segmentByCBS(dataW, avg="median")
sampleName(fitWM) <- "CBS_Example (weighted median)"
print(fitWM)
drawLevels(fitWM, col="orange", lty=3)

legend("topright", bg="white", legend=c("outliers", "non-weighted CBS (mean)", "non-weighted CBS (median)", "weighted CBS (mean)", "weighted CBS (median)"), col=c("purple", "purple", "magenta", "red", "orange"), lwd=c(NA,3,3,3,3), lty=c(NA,1,3,1,3), pch=c(1,NA,NA,NA,NA))

## Assert that weighted segment means are less biased
dmean <- getSegments(fit)$mean - getSegments(fitW)$mean
cat("Segment mean differences:\n")
print(dmean)
stopifnot(all(dmean > 0, na.rm=TRUE))

dmean <- getSegments(fitM)$mean - getSegments(fitWM)$mean
cat("Segment median differences:\n")
print(dmean)
stopifnot(all(dmean > 0, na.rm=TRUE))

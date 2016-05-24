library("calmate")

# Load example (thetaA,thetaB) signals
path <- system.file("exData", package="calmate")
theta <- loadObject("thetaAB,100x2x40.Rbin", path=path)

# Transform to (total,fracB) signals
data <- thetaAB2TotalAndFracB(theta)

# Calibrate (total,fracB) by CalMaTe
dataC <- calmateByTotalAndFracB(data)

# Calculate copy-number ratios
theta <- data[,"total",]
thetaR <- matrixStats::rowMedians(theta, na.rm=TRUE)
data[,"total",] <- 2*theta/thetaR

# Plot two "random" arrays
Clim <- c(0,4)
Blim <- c(0,1)
subplots(4, ncol=2, byrow=FALSE)
for (ii in c(1,5)) {
  sampleName <- dimnames(data)[[3]][ii]
  sampleLabel <- sprintf("Sample #%d ('%s')", ii, sampleName)
  plot(data[,,ii], xlim=Clim, ylim=Blim)
  title(main=sampleLabel)
  plot(dataC[,,ii], xlim=Clim, ylim=Blim)
  title(main=sprintf("%s\ncalibrated", sampleLabel))
}


# Assert that it also works with a single unit
dummy <- calmateByTotalAndFracB(data[1,,,drop=FALSE])
stopifnot(length(dim(dummy)) == 3)

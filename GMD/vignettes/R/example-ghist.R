## load library
require("GMD")

## create two normally-distributed samples
## with unequal means and unequal variances
set.seed(2012)
v1 <- rnorm(1000,mean=-5, sd=10)
v2 <- rnorm(1000,mean=10, sd=5)

## create common bins
n <- 20 # desired number of bins
breaks <- gbreaks(c(v1,v2),n) # bin boundaries
x <-
  list(ghist(v1,breaks=breaks,digits=0),
       ghist(v2,breaks=breaks,digits=0))
mhist.obj <- as.mhist(x)

## plot histograms side-by-side
plot(mhist.obj,mar=c(1.5,1,1,0),main="Histograms of simulated normal distributions")

## plot histograms as subplots,
## with corresponding bins aligned
plot(mhist.obj,beside=FALSE,mar=c(1.5,1,1,0),
     main="Histograms of simulated normal distributions")

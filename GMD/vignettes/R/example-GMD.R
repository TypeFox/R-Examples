
## ------------------------------------------------------------------------
## chunk1: Load library and simulate data 
## ------------------------------------------------------------------------
require("GMD") # load library

## create two normally-distributed samples
## with unequal means and unequal variances
set.seed(2012)
x1 <- rnorm(1000,mean=-5, sd=10)
x2 <- rnorm(1000,mean=10, sd=5)


## ------------------------------------------------------------------------
## chunk2: Construct histograms
## ------------------------------------------------------------------------
## create common bins
n <- 20 # desired number of bins
breaks <- gbreaks(c(x1,x2),n) # bin boundaries

## make two histograms
v1 <- ghist(x1,breaks=breaks,digits=0)
v2 <- ghist(x2,breaks=breaks,digits=0)


## ------------------------------------------------------------------------
## chunk3: Save histograms as multiple-histogram (`mhist') object
## ------------------------------------------------------------------------
x <- list(v1,v2)
mhist.obj <- as.mhist(x)


## ------------------------------------------------------------------------
## chunk4: Visualize a `mhist' object
## ------------------------------------------------------------------------
## plot histograms side-by-side
plot(mhist.obj,mar=c(1.5,1,1,0),main="Histograms of simulated normal distributions")

## plot histograms as subplots, with corresponding bins aligned
plot(mhist.obj,beside=FALSE,mar=c(1.5,1,1,0),
     main="Histograms of simulated normal distributions")


## ------------------------------------------------------------------------
## chunk5: Measure the pairwise distance betwwen two histograms by GMD
## ------------------------------------------------------------------------
gmdp.obj <- gmdp(v1,v2)
print(gmdp.obj)                       # print a brief version by default
print(gmdp.obj,mode="detailed") # print a detailed version
print(gmdp.obj,mode="full")     # print a full version


## ------------------------------------------------------------------------
## chunk6: Show alignment
## ------------------------------------------------------------------------
plot(gmdp.obj,beside=FALSE)
plot(gmdp.obj,labels=c("Distribution1","Distribution2"),beside=FALSE) # add labels



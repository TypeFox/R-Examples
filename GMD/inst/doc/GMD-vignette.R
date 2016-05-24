### R code from vignette source 'GMD-vignette.Rnw'

###################################################
### code chunk number 1: GMD-vignette.Rnw:23-24
###################################################
  cat(as.character(packageVersion('GMD')))


###################################################
### code chunk number 2: GMD-vignette.Rnw:28-29
###################################################
  cat(unlist(strsplit(packageDescription('GMD')[['Date']],' '))[1])


###################################################
### code chunk number 3: GMD-vignette.Rnw:167-174
###################################################
require("GMD") # load library

## create two normally-distributed samples
## with unequal means and unequal variances
set.seed(2012)
x1 <- rnorm(1000,mean=-5, sd=10)
x2 <- rnorm(1000,mean=10, sd=5)


###################################################
### code chunk number 4: GMD-vignette.Rnw:179-187
###################################################
## create common bins

n <- 20 # desired number of bins
breaks <- gbreaks(c(x1,x2),n) # bin boundaries

## make two histograms
v1 <- ghist(x1,breaks=breaks,digits=0)
v2 <- ghist(x2,breaks=breaks,digits=0)


###################################################
### code chunk number 5: GMD-vignette.Rnw:193-195
###################################################
x <- list(v1,v2)
mhist.obj <- as.mhist(x)


###################################################
### code chunk number 6: GMD-vignette.Rnw:199-201
###################################################
## plot histograms side-by-side
plot(mhist.obj,mar=c(1.5,1,1,0),main="Histograms of simulated normal distributions")


###################################################
### code chunk number 7: GMD-vignette.Rnw:204-207
###################################################
## plot histograms as subplots, with corresponding bins aligned
plot(mhist.obj,beside=FALSE,mar=c(1.5,1,1,0),
main="Histograms of simulated normal distributions")


###################################################
### code chunk number 8: GMD-vignette.Rnw:218-223
###################################################

gmdp.obj <- gmdp(v1,v2,sliding=TRUE)
print(gmdp.obj)                       # print a brief version by default
print(gmdp.obj,mode="detailed") # print a detailed version
print(gmdp.obj,mode="full")     # print a full version


###################################################
### code chunk number 9: GMD-vignette.Rnw:232-233
###################################################
plot(gmdp.obj,beside=FALSE)


###################################################
### code chunk number 10: GMD-vignette.Rnw:372-373 (eval = FALSE)
###################################################
## data(package="GMD")


###################################################
### code chunk number 11: GMD-vignette.Rnw:392-393 (eval = FALSE)
###################################################
## help(cage)


###################################################
### code chunk number 12: GMD-vignette.Rnw:396-405
###################################################
require(GMD)

data(cage)
class(cage)
length(cage)
names(cage)

data(cagel)
names(cagel)


###################################################
### code chunk number 13: GMD-vignette.Rnw:412-413 (eval = FALSE)
###################################################
## help(chipseq)


###################################################
### code chunk number 14: GMD-vignette.Rnw:416-423
###################################################
data(chipseq_mES)
class(chipseq_mES)
length(chipseq_mES)
names(chipseq_mES)

data(chipseq_hCD4T)
names(chipseq_hCD4T)



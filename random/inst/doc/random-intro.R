### R code from vignette source 'random-intro.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
library(random)
options(SweaveHooks=list(twofig=function() {par(mfrow=c(1,2))},
                         twofig2=function() {par(mfrow=c(2,1))},
                         onefig=function() {par(mfrow=c(1,1))}))


###################################################
### code chunk number 2: data
###################################################
## cached data to not depend on a) a network connection 
## and b) data at random.org
if ( !(file.exists("random.Rdata")) ) {
  randomOrg <- randomNumbers(n=5000, min=1, max=1e6, col=2)/1e6
  save(randomOrg, file="random.Rdata")
} else {
  load("random.Rdata")
}
#


###################################################
### code chunk number 3: <stats
###################################################
summary(randomOrg)
apply(randomOrg, 2, function(X) quantile(X,c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99), digits=4))
#


###################################################
### code chunk number 4: plotunif
###################################################
getOption("SweaveHooks")[["twofig"]]()
plot(randomOrg, ylab="", xlab="", main="5000 random,org U(0,1) draws", pch='.')
hist(matrix(randomOrg, ncol=1), xlab = "", ylab = "", main="Histogram")



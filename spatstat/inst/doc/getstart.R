### R code from vignette source 'getstart.Rnw'

###################################################
### code chunk number 1: getstart.Rnw:5-6
###################################################
options(SweaveHooks=list(fig=function() par(mar=c(1,1,1,1))))


###################################################
### code chunk number 2: getstart.Rnw:25-32
###################################################
library(spatstat)
spatstat.options(image.colfun=function(n) { grey(seq(0,1,length=n)) })
sdate <- read.dcf(file = system.file("DESCRIPTION", package = "spatstat"),
         fields = "Date")
sversion <- read.dcf(file = system.file("DESCRIPTION", package = "spatstat"),
         fields = "Version")
options(useFancyQuotes=FALSE)


###################################################
### code chunk number 3: getstart.Rnw:56-58
###################################################
getOption("SweaveHooks")[["fig"]]()
data(redwood)
plot(redwood, pch=16, main="")


###################################################
### code chunk number 4: getstart.Rnw:79-81
###################################################
getOption("SweaveHooks")[["fig"]]()
data(longleaf)
plot(longleaf, main="")


###################################################
### code chunk number 5: getstart.Rnw:138-141
###################################################
data(finpines)
mypattern <- unmark(finpines)
mydata <- round(as.data.frame(finpines), 2)


###################################################
### code chunk number 6: getstart.Rnw:156-157 (eval = FALSE)
###################################################
## mydata <- read.csv("myfile.csv")


###################################################
### code chunk number 7: getstart.Rnw:167-168
###################################################
head(mydata)


###################################################
### code chunk number 8: getstart.Rnw:183-184 (eval = FALSE)
###################################################
##   mypattern <- ppp(mydata[,3], mydata[,7], c(100,200), c(10,90))


###################################################
### code chunk number 9: getstart.Rnw:187-188 (eval = FALSE)
###################################################
## ppp(x.coordinates, y.coordinates, x.range, y.range)


###################################################
### code chunk number 10: getstart.Rnw:197-198
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(mypattern)


###################################################
### code chunk number 11: getstart.Rnw:205-206 (eval = FALSE)
###################################################
## summary(mypattern)


###################################################
### code chunk number 12: getstart.Rnw:210-211
###################################################
options(SweaveHooks=list(fig=function() par(mar=rep(4,4)+0.1)))


###################################################
### code chunk number 13: getstart.Rnw:213-214
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(Kest(mypattern))


###################################################
### code chunk number 14: getstart.Rnw:220-221 (eval = FALSE)
###################################################
## plot(envelope(mypattern,Kest))


###################################################
### code chunk number 15: getstart.Rnw:223-224
###################################################
env <- envelope(mypattern,Kest, nsim=39)


###################################################
### code chunk number 16: getstart.Rnw:226-227
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(env, main="envelope(mypattern, Kest)")


###################################################
### code chunk number 17: getstart.Rnw:229-230
###################################################
options(SweaveHooks=list(fig=function() par(mar=c(1,1,1,1))))


###################################################
### code chunk number 18: getstart.Rnw:236-237
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(density(mypattern))


###################################################
### code chunk number 19: getstart.Rnw:247-248 (eval = FALSE)
###################################################
## marks(mypattern) <- mydata[, c(5,9)]


###################################################
### code chunk number 20: getstart.Rnw:250-251
###################################################
mypattern <-finpines


###################################################
### code chunk number 21: getstart.Rnw:254-255 (eval = FALSE)
###################################################
## plot(Smooth(mypattern))


###################################################
### code chunk number 22: getstart.Rnw:258-259
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(Smooth(mypattern, sigma=1.2), main="Smooth(mypattern)")



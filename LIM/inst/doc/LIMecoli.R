### R code from vignette source 'LIMecoli.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: preliminaries
###################################################
library("LIM")
options(prompt = "> ")
options(width=70)


###################################################
### code chunk number 2: LIMecoli.Rnw:361-362
###################################################
pars <- Ldei(LIMEcoli)


###################################################
### code chunk number 3: LIMecoli.Rnw:371-372
###################################################
LP <- Linp(LIMEcoli)


###################################################
### code chunk number 4: LIMecoli.Rnw:375-376
###################################################
LP$X <- LP$X[1,]


###################################################
### code chunk number 5: LIMecoli.Rnw:379-380
###################################################
xr <- Xranges(LIMEcoli)


###################################################
### code chunk number 6: LIMecoli.Rnw:385-386
###################################################
data.frame(simplest = pars$X, optimal = LP$X, xr)


###################################################
### code chunk number 7: range
###################################################
par(mfrow = c(1, 2))
nr <- LIMEcoli$NUnknowns
ii <- 1:(nr/2)
dotchart(LP$X[ii], xlim = range(xr), pch = 16, cex = 0.8)
segments(xr[ii,1], 1:nr, xr[ii,2], 1:nr)
ii <- (nr/2+1):nr
dotchart(LP$X[ii], xlim = range(xr), pch = 16, cex = 0.8)
segments(xr[ii,1], 1:nr, xr[ii,2], 1:nr)
mtext(side =  3, cex = 1.5, outer = TRUE, line = -1.5,
      "E coli Core Metabolism, optimal solution and ranges")


###################################################
### code chunk number 8: figrange
###################################################
par(mfrow = c(1, 2))
nr <- LIMEcoli$NUnknowns
ii <- 1:(nr/2)
dotchart(LP$X[ii], xlim = range(xr), pch = 16, cex = 0.8)
segments(xr[ii,1], 1:nr, xr[ii,2], 1:nr)
ii <- (nr/2+1):nr
dotchart(LP$X[ii], xlim = range(xr), pch = 16, cex = 0.8)
segments(xr[ii,1], 1:nr, xr[ii,2], 1:nr)
mtext(side =  3, cex = 1.5, outer = TRUE, line = -1.5,
      "E coli Core Metabolism, optimal solution and ranges")


###################################################
### code chunk number 9: LIMecoli.Rnw:429-433
###################################################
print(system.time(
  xs <- Xsample(LIMEcoli, iter = 500, type = "mirror", test = TRUE)  #))
))
  


###################################################
### code chunk number 10: sample
###################################################
  pairs(xs[, 1:12], pch = ".", cex = 2, gap = 0, upper.panel = NULL)


###################################################
### code chunk number 11: figsample
###################################################
  pairs(xs[, 1:12], pch = ".", cex = 2, gap = 0, upper.panel = NULL)



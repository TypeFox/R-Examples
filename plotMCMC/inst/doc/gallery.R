### R code from vignette source 'gallery.Rnw'

###################################################
### code chunk number 1: require
###################################################
require(plotMCMC)


###################################################
### code chunk number 2: Auto1
###################################################
plotAuto(xpar$R0)


###################################################
### code chunk number 3: Auto2
###################################################
plotAuto(xpar$R0, thin=10)


###################################################
### code chunk number 4: Auto3
###################################################
plotAuto(xpar, lag.max=50, ann=FALSE, axes=FALSE)


###################################################
### code chunk number 5: Cumu1
###################################################
plotCumu(xpar$R0, main="R0")


###################################################
### code chunk number 6: Cumu2
###################################################
plotCumu(xpar$cSfull, main="cSfull")


###################################################
### code chunk number 7: Cumu3
###################################################
plotCumu(xpar, probs=c(0.25,0.75), ann=FALSE, axes=FALSE)


###################################################
### code chunk number 8: Dens1
###################################################
plotDens(xbio$"2004", points=TRUE, div=1000, main="2004\n",
         xlab="Biomass age 4+ (1000 t)", tick.number=6, strip=FALSE)


###################################################
### code chunk number 9: Dens2
###################################################
plotDens(xpar, xlab="Parameter value", ylab="Posterior density\n")


###################################################
### code chunk number 10: Quant1
###################################################
plotQuant(xrec, names=substring(names(xrec),3), div=1000, xlab="Year",
          ylab="Recruitment (million one-year-olds)")


###################################################
### code chunk number 11: Quant2
###################################################
plotQuant(xbio, div=1000, xlab="Year", ylab="Biomass age 4+ (kt)")


###################################################
### code chunk number 12: Quant3
###################################################
plotQuant(xbio, style="bars", div=1000, sfrac=0, xlab="Year",
          ylab="Biomass age 4+ (kt)")


###################################################
### code chunk number 13: Quant4
###################################################
plotQuant(xbio, style="lines", div=1000, xlab="Year",
          ylab="Biomass age 4+ (kt)")


###################################################
### code chunk number 14: Quant5
###################################################
plotQuant(xpro, axes=1:2, div=1000, xlab="Year",
          ylab="Biomass age 4+ (kt)")


###################################################
### code chunk number 15: Splom1
###################################################
plotSplom(xpar, pch=".")


###################################################
### code chunk number 16: Splom2
###################################################
plotSplom(xpro, axes=TRUE, between=1, div=1000, main="Future biomass",
          cex.labels=1.5, pch=".", cex=3)


###################################################
### code chunk number 17: Trace1
###################################################
plotTrace(xpar, xlab="Iterations", ylab="Parameter value",
          layout=c(2,4))


###################################################
### code chunk number 18: Trace2
###################################################
plotTrace(xpar$R0, axes=TRUE, div=1000)



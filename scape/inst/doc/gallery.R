### R code from vignette source 'gallery.Rnw'

###################################################
### code chunk number 1: require
###################################################
require(scape)


###################################################
### code chunk number 2: B1
###################################################
plotB(x.ling, series=c("VB.1","VB.2","Y"), div=1000, xlab="Year\n",
      ylab="Biomass and landings (1000 t)")


###################################################
### code chunk number 3: B2
###################################################
plotB(x.ling, "s", div=1000, xlab="Biomass age 4+ (1000 t)",
      ylab="Recruitment (million one-year-olds)")


###################################################
### code chunk number 4: CA1
###################################################
plotCA(x.sbw, fit=FALSE, strip=FALSE, xlab="Age", ylab="Year",
       tick.number=10)


###################################################
### code chunk number 5: CA2
###################################################
plotCA(x.cod, xlab="Age", ylab="Proportion in catch", cex.strip=0.7,
       cex.axis=0.7, col.lines="brown", layout=c(8,4))


###################################################
### code chunk number 6: CA3
###################################################
plotCA(x.cod, xlab="Age", ylab="Proportion in catch", cex.strip=0.7,
       cex.axis=0.7, col.lines="brown", layout=c(2,4), swap=TRUE,
       ages=3:10, same.limits=FALSE)


###################################################
### code chunk number 7: CA4
###################################################
plotCA(x.ling, "s", col.points=c("red","blue"), lty.lines=0, xlab="Age",
       ylab="Observed proportion in survey", tck=0.5, cex.strip=0.7,
       cex.axis=0.7)


###################################################
### code chunk number 8: CA5
###################################################
plotCA(x.ling, "s", xlab="Age", ylab="Observed proportion in survey",
       fit=FALSE, cex.strip=0.7, cex.axis=0.7, tck=0.5, layout=c(5,2))


###################################################
### code chunk number 9: CA6
###################################################
plotCA(x.ling, "s", xlab="Age", ylab="Observed proportion in survey",
       fit=FALSE, cex.strip=0.7, cex.axis=0.7, tck=0.5, layout=c(5,6),
       swap=TRUE)


###################################################
### code chunk number 10: CL1
###################################################
plotCL(x.ling, fit=FALSE, strip=FALSE, series="1", sex="Female",
       xlab="Length (cm)", ylab="Year")


###################################################
### code chunk number 11: CL2
###################################################
plotCL(x.oreo, xlab="Length (cm)", ylab="Proportion in catch")


###################################################
### code chunk number 12: CL3
###################################################
plotCL(x.oreo, "s", layout=c(2,1), xlab="Length (cm)",
       ylab="Observed proportion in survey", cex.points=0.8,
       col.points=c("red","blue"), lty.lines=0)


###################################################
### code chunk number 13: CL4
###################################################
plotCL(x.ling, fit=FALSE, series="2", xlab="Length (cm)",
       ylab="Observed proportion in trawl catch", tck=0.5)


###################################################
### code chunk number 14: CL5
###################################################
plotCL(x.ling, series="2", swap=TRUE, lengths=70:150, lty.grid=0)


###################################################
### code chunk number 15: Index1
###################################################
plotIndex(x.cod, xlab="Year", ylab="Survey abundance index",
          strip=FALSE)


###################################################
### code chunk number 16: Index2
###################################################
plotIndex(x.oreo, "c", series="Series 1-1", xlim=c(1981,1990))


###################################################
### code chunk number 17: Index3
###################################################
plotIndex(x.oreo, "c", xlim=list(c(1981,1990),c(1992,2002)),
          xlab="Year", ylab="Observed CPUE",
          col.points=c("salmon","seagreen"), lty.lines=0)


###################################################
### code chunk number 18: LA1
###################################################
plotLA(x.oreo, xlab="Age", ylab="Length (cm)")


###################################################
### code chunk number 19: LA2
###################################################
mykey <- list(text=list(lab=c("Female","Male")), space="right",
              lines=list(lwd=4,col=c("red","blue")))
plotLA(x.oreo, together=TRUE, xlab="Age", ylab="Length (cm)", pch=NA,
       key=mykey)


###################################################
### code chunk number 20: LA3
###################################################
mykey <- list(text=list(lab=c("Female","Male")), space="right",
              points=list(pch=16,cex=0.5,col=c("red","blue")))
plotLA(x.oreo, together=TRUE, xlab="Age", ylab="Length (cm)",
       col.points=c("red","blue"), lty.lines=0, key=mykey)


###################################################
### code chunk number 21: N1
###################################################
plotN(x.cod, div=1000, xlab=c("Age (years)","Year"),
      ylab="Individuals (million)")


###################################################
### code chunk number 22: N2
###################################################
plotN(x.cod, "l", div=1000, xlab="Age", ylab="Individuals (million)")


###################################################
### code chunk number 23: N3
###################################################
plotN(x.cod, "r", age=3, div=1000, xlim=c(1967,2002))


###################################################
### code chunk number 24: N4
###################################################
plotN(x.cod, "p", div=1000, ages=3:10, xlim=c(2,11), xlab="Age",
      ylab="Individuals (million)", cex.strip=0.7, cex.axis=0.7,
      tck=0.5)


###################################################
### code chunk number 25: N5
###################################################
plotN(x.cod, "b", xlab="Age (years)", ylab="Year", cex.points=0.7)


###################################################
### code chunk number 26: Sel1
###################################################
plotSel(x.ling, xlab="Age", ylab="Selectivity and maturity")


###################################################
### code chunk number 27: Sel2
###################################################
plotSel(x.cod, together=TRUE, xlab="Age\n", ylab="Selectivity",
        pch=NA, col.lines=c("coral","navyblue"), strip=FALSE)



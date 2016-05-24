### R code from vignette source 'fishervig.Rnw'

###################################################
### code chunk number 1: setup
###################################################
ignore <- require(aylmer,quietly=TRUE)
ignore <- require(partitions,quietly=TRUE)


###################################################
### code chunk number 2: load_iqd
###################################################
data(iqd)  # Table 1


###################################################
### code chunk number 3: test_iqd
###################################################
jj <- aylmer.test(iqd)


###################################################
### code chunk number 4: print_iqd
###################################################
jj


###################################################
### code chunk number 5: test_iqd_alternative
###################################################
jj1 <- aylmer.test(iqd, alternative=function(x){x[1,1]})


###################################################
### code chunk number 6: print_iqd_test
###################################################
jj1


###################################################
### code chunk number 7: aylmerTestA
###################################################
shifts
aylmer.test(shifts)


###################################################
### code chunk number 8: fisherA
###################################################
fisher.test(shifts[-1,  ])$p.value
fisher.test(shifts[  ,-3])$p.value


###################################################
### code chunk number 9: loadIcons
###################################################
set.seed(0)
data(icons)


###################################################
### code chunk number 10: useGood
###################################################
good(icons)


###################################################
### code chunk number 11: calculateIconStats
###################################################
iconstats <- aylmer.test(icons, simulate.p.value=TRUE)


###################################################
### code chunk number 12: printIconStats
###################################################
iconstats


###################################################
### code chunk number 13: dataPurum
###################################################
data(purum)


###################################################
### code chunk number 14: randomPurum
###################################################
set.seed(1)
rP <- randomprobs(purum,B=800)


###################################################
### code chunk number 15: fishervig.Rnw:598-599
###################################################
g <- function(x) max(abs(x-t(x)),na.rm=TRUE)


###################################################
### code chunk number 16: purum_gender
###################################################
set.seed(1)
jj7 <- aylmer.test(purum, alternative=g, simulate.p.value=TRUE,B=2000)


###################################################
### code chunk number 17: print_iqd_test
###################################################
jj7


###################################################
### code chunk number 18: plotPurum
###################################################
plot(rP, xlab="index", ylab="log(Prob)", type="o", pch=16, cex=0.4, axes=FALSE)
axis(side=1,pos=min(rP)-1)
axis(side=2,pos= -20)
val <-  -175.5185
segments(x0=0,x1=800,y0=val,y1=val,col="gray",lwd=2)
points(x = 1, y=rP[1],pch=16,cex=1.5)
text(x = 1, y=rP[1],"observation",pos=4)


###################################################
### code chunk number 19: fishervig.Rnw:769-770
###################################################
data(frogs)


###################################################
### code chunk number 20: CalculateFrogStats
###################################################
frogResults <- aylmer.test(frogs, simulate.p.value=TRUE,B=2000)


###################################################
### code chunk number 21: printFrogStats
###################################################
frogResults


###################################################
### code chunk number 22: printFrogsMatrix
###################################################
frogs.matrix


###################################################
### code chunk number 23: fishervig.Rnw:916-918
###################################################
data(gear) # Table 6
aylmer.test(gear)


###################################################
### code chunk number 24: calcGear
###################################################
jj <- aylmer.test(gear,alternative="less")$p.value
jj <- round(jj*1e8)/1e3



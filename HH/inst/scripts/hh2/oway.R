### R code from vignette source '~/WindowsC/HOME/rmh/hh.e2/hh2/oway.tex'

###################################################
### code chunk number 1: oway.tex:8-9
###################################################
library(HH)


###################################################
### code chunk number 2: oway.tex:12-17
###################################################
## the standard lattice color 2 is difficult for people with color deficient vision
data(col3x2)
## These colors look like a 3x2 color array when run through
## the vischeck simulator to see how they look for the three most
## common color vision deficiencies: Protanope, Deuteranope, Tritanope.


###################################################
### code chunk number 3: oway.tex:70-77
###################################################
## hhpdf("catalystm1.pdf", width=5, height=3, col=col3x2) ## col is not an argument for grDevices:::pdf
data(catalystm)
bwplot(concent ~ catalyst, data=catalystm,
       panel=panel.bwplot.superpose, groups=catalyst,
       ylab=list("concentration"),
       xlab=list("catalyst"))
## hhdev.off()


###################################################
### code chunk number 4: oway.tex:106-111
###################################################
## hhcapture("catalystm-aov1.Rout", '
catalystm1.aov <- aov(concent ~ catalyst, data=catalystm)
anova(catalystm1.aov)
model.tables(catalystm1.aov, "means")
## ')


###################################################
### code chunk number 5: oway.tex:480-485
###################################################
## hhcapture("catalystm-aov2.Rout", '
catalystm.mmc <-
   mmc(catalystm1.aov, linfct = mcp(catalyst = "Tukey"))
catalystm.mmc
## ')


###################################################
### code chunk number 6: oway.tex:526-529
###################################################
## hhpdf("catalystm-mmc-mmc.pdf", width=6, height=6)
mmcplot(catalystm.mmc, style="both")
## hhdev.off()


###################################################
### code chunk number 7: oway.tex:723-728
###################################################
## hhcapture("batch.Rout", '
data(batch)
batch1.aov <- aov(Calcium ~ Batch, data=batch)
anova(batch1.aov)
## ')


###################################################
### code chunk number 8: oway.tex:745-750
###################################################
old.par <- options(width=72, digits=6)
## hhcapture("batchhov.Rout", '
hovBF(Calcium ~ Batch, data=batch)
## ')
options(old.par)


###################################################
### code chunk number 9: oway.tex:765-768
###################################################
## hhpdf("batchhov.pdf", width=7, height=4, col=col3x2) ## col is not an argument for grDevices:::pdf
hovplotBF(Calcium ~ Batch, data=batch)
## hhdev.off()


###################################################
### code chunk number 10: oway.tex:785-788
###################################################
## hhpdf("batchoway.pdf", width=5, height=4, col=col3x2) ## col is not an argument for grDevices:::pdf
OneWayVarPlot(Calcium ~ Batch, data = batch)
## hhdev.off()


###################################################
### code chunk number 11: oway.tex:913-918
###################################################
## hhpdf("turkey-f1.pdf", width=6, height=4, col=col3x2) ## col is not an argument for grDevices:::pdf)
data(turkey)
bwplot(wt.gain ~ diet, data=turkey, groups=diet,
       panel=panel.bwplot.superpose, xlab="diet", ylab="Weight Gain")
## hhdev.off()


###################################################
### code chunk number 12: oway.tex:929-934
###################################################
## hhcapture("turkey-aov1.Rout", '
turkey.aov <- aov(wt.gain ~ diet, data=turkey)
summary(turkey.aov)
model.tables(turkey.aov, type="means", se=TRUE)
## ')


###################################################
### code chunk number 13: oway.tex:972-983
###################################################
## hhcapture("turkey-contrasts.Rout", '
contrasts(turkey$diet)
contrasts(turkey$diet) <-
  cbind(control.vs.treatment=c(1,-.25,-.25,-.25,-.25),
        A.vs.B              =c(0, .5,  .5, -.5, -.5 ),
        amount              =c(0, .5, -.5,  .5, -.5 ),
        A.vs.B.by.amount    =c(0, .5, -.5, -.5,  .5 ))
contrasts(turkey$diet)
tapply(turkey$wt.gain, turkey$diet, mean) %*%
   contrasts(turkey$diet)
## ')


###################################################
### code chunk number 14: oway.tex:1001-1013
###################################################
## hhcapture("turkey-anova-contrasts.Rout", '
turkey2.aov <- aov(wt.gain ~ diet, data=turkey)
summary(turkey2.aov)
old.width <- options(width=67)
summary(turkey2.aov,
        split=list(diet=list(
                     control.vs.treatment=1,
                     A.vs.B=2,
                     amount=3,
                     A.vs.B.by.amount=4)))
options(old.width)
## ')


###################################################
### code chunk number 15: oway.tex:1417-1421
###################################################
## hhpdf("catalystm-hov.pdf", width=7, height=4, col=col3x2) ## col is not an argument for grDevices:::pdf
hovBF(concent ~ catalyst, data=catalystm)
hovplotBF(concent ~ catalyst, data=catalystm)
## hhdev.off()


###################################################
### code chunk number 16: oway.tex:1721-1731
###################################################
## hhcapture("SumsSquareIdentities.Rout", '
data(catalystm)
catalystm.aov <- aov(concent ~ catalyst, data=catalystm)
anova(catalystm.aov)
model.tables(catalystm.aov)
Proj <- proj(catalystm.aov)
Proj <-cbind(Proj, Sum=apply(Proj, 1, sum))
Proj
apply(Proj, 2, function(x) sum(x^2))
## ')


###################################################
### code chunk number 17: oway.tex:1757-1766
###################################################
## hhcapture("ANOVAbyRegression.Rout", '
contrasts(catalystm$catalyst)
X <- model.matrix(catalystm.aov)[,2:4]
X
catalystm.lm <-
    lm(concent ~ X[,"catalystB"] + X[,"catalystC"] + X[,"catalystD"],
       data=catalystm)
anova(catalystm.lm)
## ')


###################################################
### code chunk number 18: oway.tex:1784-1796
###################################################
## hhcapture("oopplot.Rout", '
tmp <- data.frame(AA=c(5,6,8,7,8),
                  BB=factor(letters[c(5,6,8,7,8)]),
                  CC=ts(c(5,6,8,7,8)),
                  stringsAsFactors=FALSE)
tmp
sapply(tmp, class)
is.numeric(tmp$A)
plot(tmp$AA)
plot(tmp$BB)
plot(tmp$CC)
## ')


###################################################
### code chunk number 19: oway.tex:1798-1805
###################################################
## hhpdf("oopplot.pdf", width=9, height=3.25)
old.par <- par(mfrow=c(1,3), cex=1.25, mar=c(4,2,1,5))
plot(tmp$AA, ylab=" ", xlim=c(.5, 5.5), ylim=c(4.5, 8.5))
plot(tmp$BB, ylab=NULL)
plot(tmp$CC, ylab=NULL)
par(old.par)
## hhdev.off()


###################################################
### code chunk number 20: oway.tex:1831-1838
###################################################
## hhcapture("oopaov.Rout", '
class(catalystm.aov)
summary(catalystm.aov)
old.par <- par(mfrow=c(1,4))
plot(catalystm.aov)
par(old.par)
## ')


###################################################
### code chunk number 21: oway.tex:1840-1845
###################################################
## hhpdf("oopaov.pdf", width=12, height=4)
old.par <- par(mfrow=c(1,4), cex=1.25)
plot(catalystm.aov)
par(old.par)
## hhdev.off()



### R code from vignette source '~/WindowsC/HOME/rmh/hh.e2/hh2/mcomp.tex'

###################################################
### code chunk number 1: mcomp.tex:8-9
###################################################
library(HH)


###################################################
### code chunk number 2: mcomp.tex:12-17
###################################################
## the standard lattice color 2 is difficult for people with color deficient vision
data(col3x2)
## These colors look like a 3x2 color array when run through
## the vischeck simulator to see how they look for the three most
## common color vision deficiencies: Protanope, Deuteranope, Tritanope.


###################################################
### code chunk number 3: mcomp.tex:179-188
###################################################
## hhpdf("weightloss-data.pdf", width=7, height=4, col=likertColor(2)[2]) ## col is not an argument for grDevices:::pdf
data(weightloss)
bwplot(loss ~ group, data=weightloss,
       scales=list(cex=1.5),
       ylab=list("Weight Loss", cex=1.5),
       xlab=list("group",cex=1.5),
       par.settings=list(box.dot=list(
          col=trellis.par.get()$superpose.symbol$col[1])))
## hhdev.off()


###################################################
### code chunk number 4: mcomp.tex:205-209
###################################################
## hhcapture("weightloss.Rout", '
weightloss.aov <- aov(loss ~ group, data=weightloss)
summary(weightloss.aov)
## ')


###################################################
### code chunk number 5: mcomp.tex:289-297
###################################################
## hhcapture("weightloss-b.Rout", '
weightloss.dunnett <-
glht(weightloss.aov,
     linfct=mcp(group=
                contrMat(table(weightloss$group), base=4)),
     alternative = "greater")
confint(weightloss.dunnett)
## ')


###################################################
### code chunk number 6: mcomp.tex:312-315
###################################################
## hhpdf("weightloss-dunnett.pdf", width=7, height=2.75)
mmcplot(weightloss.dunnett, focus="group")
## hhdev.off()


###################################################
### code chunk number 7: mcomp.tex:330-338
###################################################
## hhcapture("weightloss-dunnet-mmc.Rout", '
weightloss.mmc <-
  mmc(weightloss.aov,
      linfct=mcp(group=
                 contrMat(table(weightloss$group), base=4)),
      alternative = "greater")
weightloss.mmc
## ')


###################################################
### code chunk number 8: mcomp.tex:352-355
###################################################
## hhpdf("weightloss-dunnet-mmc.pdf", width=7, height=7)
mmcplot(weightloss.mmc, h=c(.80, .20), style="both")
## hhdev.off()


###################################################
### code chunk number 9: mcomp.tex:494-509
###################################################
## hhcapture("turkey-contrasts2a.Rout", '
data(turkey)
turkey.aov <- aov(wt.gain ~ diet, data=turkey)
scheffe.quantile <- sqrt(4*qf(.95, 4, 25))
turkey.lmat <-
  cbind(control.vs.treatment=c(1,-.25,-.25,-.25,-.25),
        A.vs.B              =c(0, .5,  .5, -.5, -.5 ),
        amount              =c(0, .5, -.5,  .5, -.5 ),
        A.vs.B.by.amount    =c(0, .5, -.5, -.5,  .5 ))
row.names(turkey.lmat) <- row.names(contrasts(turkey$diet))
turkey.mmc <- mmc(turkey.aov, calpha=scheffe.quantile, focus="diet",
                  focus.lmat=turkey.lmat,
                  estimate.sign=0, order.contrasts=FALSE)
turkey.mmc$lmat
## ')


###################################################
### code chunk number 10: mcomp.tex:543-546
###################################################
## hhpdf("turkey-scheffelmat.pdf", width=6.5, height=2.5)
mmcplot(turkey.mmc, type="lmat", style="confint", axis.right=2.3)
## hhdev.off()


###################################################
### code chunk number 11: mcomp.tex:710-713
###################################################
## hhpdf("turkey-scheffe.pdf", width=7.5, height=7.5)
mmcplot(turkey.mmc, style="both")
## hhdev.off()


###################################################
### code chunk number 12: mcomp.tex:760-773
###################################################
## hhpdf("turkey-lmat.pdf", width=8, height=7.5)
turkey.lmat <-
  cbind(control.vs.treatment=c(1,-.25,-.25,-.25,-.25),
        A.vs.B              =c(0, .5,  .5, -.5, -.5 ),
        amount              =c(0, .5, -.5,  .5, -.5 ),
        A.vs.B.by.amount    =c(0, .5, -.5, -.5,  .5 ))
row.names(turkey.lmat) <- row.names(contrasts(turkey$diet))
turkey.lmat ## these are the constructed contrasts
turkey.mmc <- mmc(turkey.aov, calpha=scheffe.quantile, focus="diet",
                  focus.lmat=turkey.lmat)
mmcplot(turkey.mmc, type="lmat", style="both")
turkey.mmc
## hhdev.off()


###################################################
### code chunk number 13: mcomp.tex:881-888
###################################################
data(catalystm)
catalystm1.aov <- aov(concent ~ catalyst, data=catalystm)
## hhcapture("catalystm-glht.Rout", '
catalystm.glht <-
   glht(catalystm1.aov, linfct = mcp(catalyst = "Tukey"))
confint(catalystm.glht)
## ')


###################################################
### code chunk number 14: mcomp.tex:899-902
###################################################
## hhpdf("catalystm-glht.pdf", width=7, height=3)
mmcplot(catalystm.glht, order.contrasts=FALSE, estimate.sign=0, focus="catalyst")
## hhdev.off()


###################################################
### code chunk number 15: mcomp.tex:952-957
###################################################
## hhpdf("catalystm-tukey-lines.pdf", width=7, height=3)
plot(cld(catalystm.glht))
axis(side="3", at=1:4, labels= signif(model.tables(catalystm1.aov, "means")[[1]]$catalyst, 4) )
segments(c(.5, 2.5), c(62.1, 63.7), c(3.5, 4.5), c(62.1, 63.7), xpd=NA, lty=3, lwd=2)
## hhdev.off()


###################################################
### code chunk number 16: mcomp.tex:1050-1063
###################################################
## hhcapture("inconsistent.Rout", '
group <- factor(LETTERS[1:4])
n <- c(5,100,100,5)
ybar <- c(2, 2.1, 2.8, 3)

inconsistent.aov <- aovSufficient(ybar ~ group, weights=n, sd=.8)
anova(inconsistent.aov)
inconsistent.glht <-
   glht(inconsistent.aov, linfct=mcp(group="Tukey"),
        vcov.=vcovSufficient, df=inconsistent.aov$df.residual)
crit.point <- qtukey(.95, 4, 206)/sqrt(2)
confint(inconsistent.glht, calpha=crit.point)
## ')


###################################################
### code chunk number 17: mcomp.tex:1082-1087
###################################################
## hhpdf("inconsistent.pdf", width=7, height=3)
plot(cld(inconsistent.glht))
axis(side="3", at=1:4, labels=format(ybar,1))
segments(c(.5, .5), c(3.43, 3.61), c(4.5, 4.5), c(3.43, 3.61), xpd=NA, lty=3, lwd=2)
## hhdev.off()


###################################################
### code chunk number 18: mcomp.tex:1101-1104
###################################################
## hhpdf("inconsistent-glht.pdf", width=7, height=3)
mmcplot(inconsistent.glht, focus="group")
## hhdev.off()


###################################################
### code chunk number 19: mcomp.tex:1123-1131
###################################################
## hhpdf("inconsistent-none.pdf", width=7, height=2.5)
inconsistent.mmc <- mmc(inconsistent.aov,
                        linfct=mcp(group="Tukey"),
                        vcov.=vcovSufficient,
                        df=inconsistent.aov$df.residual,
                        calpha=qtukey(.95, 4, 206)/sqrt(2))
mmcplot(inconsistent.mmc$none$glht, focus="group")
## hhdev.off()


###################################################
### code chunk number 20: mcomp.tex:1148-1151
###################################################
## hhpdf("inconsistent-mmc.pdf", width=7, height=5)
mmcplot(inconsistent.mmc)
## hhdev.off()


###################################################
### code chunk number 21: mcomp.tex:1199-1204
###################################################
## hhpdf("catalystm-mmc-mca.pdf", width=6, height=4)
catalystm.mmc <-
   mmc(catalystm1.aov, linfct = mcp(catalyst = "Tukey"))
mmcplot(catalystm.mmc)
## hhdev.off()


###################################################
### code chunk number 22: mcomp.tex:1256-1276
###################################################
## just the B-D contrast
`lmat.B-D` <- catalystm.mmc$mca$lmat[,"B-D", drop=FALSE]
dimnames(`lmat.B-D`)[[1]] <- levels(catalystm$catalyst)
catalystmBD.mmc <- mmc(catalystm1.aov, focus.lmat=`lmat.B-D`)
print(catalystmBD.mmc)

group <- levels(catalystm$catalyst)
n <- c(4,4,4,4)
ybar <- tapply(catalystm$concent, catalystm$catalyst, mean)
ms.5 <- summary.aov(catalystm1.aov)[[1]][2,"Mean Sq"]^.5
crit.point <- catalystm.mmc$mca$crit.point
xlim.explain <- c(45.5,62.5)

## hhpdf("mmc1-a.pdf", width=7, height=7)
HH:::mmc.explain(group, n, ybar, ms.5, crit.point,
                 ylabel="concent", factor.label="catalyst",
                 xlim=xlim.explain,
                 exit=1,
                 col=c("gray30", "navyblue", "SlateBlue", "black", "royalblue", "red"))
## hhdev.off()


###################################################
### code chunk number 23: mcomp.tex:1322-1336
###################################################
## This figure showing steps 7 and 8 is not displayed in the book.
HH:::mmc.explain(group, n, ybar, ms.5, crit.point,
                 ylabel="concent", factor.label="catalyst",
                 xlim=xlim.explain,
                 exit=2)
## col=c( "gray20", "navyblue", "SlateBlue", "gray50", "gray20", "red") ## exit 2, 3

## hhpdf("mmc1-b.pdf", width=7, height=7)
HH:::mmc.explain(group, n, ybar, ms.5, crit.point,
                 ylabel="concent", factor.label="catalyst",
                 xlim=xlim.explain,
                 exit=3)
## col=c( "gray20", "navyblue", "SlateBlue", "gray50", "gray20", "red") ## exit 2, 3
## hhdev.off()


###################################################
### code chunk number 24: mcomp.tex:1385-1396
###################################################
## This figure showing step 11 (actually steps 1 through 9, and 11)
## is not displayed in the book.
HH:::mmc.explain(group, n, ybar, ms.5, crit.point,
                 ylabel="concent", factor.label="catalyst",
                 xlim=xlim.explain,
                 exit=4)
## col=c("gray20", "gray10", "gray70", "gray50", "gray20", "red")

## hhpdf("mmc2.pdf", width=6, height=4)
mmcplot(catalystmBD.mmc, type="lmat")
## hhdev.off()


###################################################
### code chunk number 25: mcomp.tex:1695-1712
###################################################
## hhcapture("pairwiseContrasts.Rout", '
## aov contrast matrix for catalyst factor.  The columns are
## constructed by contr.treatment with the default base=1
contrasts(catalystm$catalyst)

## Linear function used internally by glht for pairwise contrasts.
## The rows of linfct are the differences of the columns
## of the contrast matrix.
catalystm.mmc$mca$glht$linfct

## Contrasts in lmat format, each column sums to zero.
## The last three rows are the transpose of the last three columns
## of the linfct matrix.
## The first row is prepended to make the column sum be zero.
catalyst.pairwise <- lmatPairwise(catalystm.mmc)
catalyst.pairwise
## ')


###################################################
### code chunk number 26: mcomp.tex:1745-1757
###################################################
## hhcapture("catalystm-mmc.orth-matrix.Rout", '
## An orthogonal set of ($4-1$) contrasts for the catalyst factor.
## user-specified contrasts       A  B  C  D
catalystm.lmat <- cbind("AB-D" =c(1, 1, 0,-2),
                        "A-B"  =c(1,-1, 0, 0),
                        "ABD-C"=c(1, 1,-3, 1))
dimnames(catalystm.lmat)[[1]] <- levels(catalystm$catalyst)
catalystm.lmat
crossprod(catalystm.lmat)
catalyst.pairwise
resid(lm(catalystm.lmat ~ catalyst.pairwise))
## ')


###################################################
### code chunk number 27: mcomp.tex:1791-1798
###################################################
## hhcapture("mmc-complete.Rout", '
catalystm.mmc <-
   mmc(catalystm1.aov,
       linfct = mcp(catalyst = "Tukey"),
       focus.lmat=catalystm.lmat)
catalystm.mmc
## ')


###################################################
### code chunk number 28: mcomp.tex:1812-1815
###################################################
## hhpdf("catalystm-mmc-lmat.pdf", width=6, height=6)
mmcplot(catalystm.mmc, type="lmat", style="both")
## hhdev.off()


###################################################
### code chunk number 29: mcomp.tex:1865-1873
###################################################
## hhcapture("pulmonary.Rout", '
data(pulmonary)
pulmonary
pulmonary.aov <-
 aovSufficient(FVC ~ smoker, data=pulmonary,
               weights=pulmonary$n, sd=pulmonary$s)
summary(pulmonary.aov)
## ')


###################################################
### code chunk number 30: mcomp.tex:1942-1957
###################################################
pulm.lmat <- cbind("npnl-mh"=c( 1, 1, 1, 1,-2,-2), ## not.much vs lots
                   "n-pnl"  =c( 3,-1,-1,-1, 0, 0), ## none vs light
                   "p-nl"   =c( 0, 2,-1,-1, 0, 0), ## {} arbitrary 2 df
                   "n-l"    =c( 0, 0, 1,-1, 0, 0), ## {} for 3 types of light
                   "m-h"    =c( 0, 0, 0, 0, 1,-1)) ## moderate vs heavy
dimnames(pulm.lmat)[[1]] <- row.names(pulmonary)
pulmonary.mmc <-
       mmc(pulmonary.aov,
           linfct=mcp(smoker="Tukey"),
           df=pulmonary.aov$df.residual,
           vcov.=vcovSufficient,
           focus.lmat=pulm.lmat)
## hhpdf("pulmonary-mmc-mca.pdf", width=8, height=10)
mmcplot(pulmonary.mmc, style="both")
## hhdev.off()


###################################################
### code chunk number 31: mcomp.tex:1972-1975
###################################################
## hhpdf("pulmonary-mmc-lmat.pdf", width=8, height=10)
mmcplot(pulmonary.mmc, type="lmat", style="both")
## hhdev.off()



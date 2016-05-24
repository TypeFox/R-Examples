### R code from vignette source '~/WindowsC/HOME/rmh/hh.e2/hh2/regb.tex'

###################################################
### code chunk number 1: regb.tex:9-10
###################################################
library(HH)


###################################################
### code chunk number 2: regb.tex:94-129
###################################################
## rgl graphics doesn't automate.  It needs manual intervention for screen shots.
data(fat)
## car::scatter3d(bodyfat~abdomin+biceps, data=fat, fit="linear",
##                residuals="squares",
##                bg="white", axis.scales=TRUE, grid=TRUE, ellipsoid=FALSE)

fat2.resid <- resid(lm(bodyfat ~ abdomin + biceps, data=fat))
car::scatter3d(bodyfat ~ abdomin + biceps, data=fat, fit="linear",
               residuals="squares",
               bg="white", axis.scales=TRUE, grid=TRUE, ellipsoid=FALSE,
               square.color = "gray80", surface.col="#a6cafe",
               surface.alpha=.3, sphere.size=.7,
               point.col=c("red","green")[1+(fat2.resid >= 0)])
rgl::par3d(windowRect=c(1803, 57, 2646, 930)) ## Macintosh

rgl::par3d(userMatrix=structure(
             c(0.5226189494133, 0.0857846289873123, -0.848239600658417, 0,
               0, 0.994925022125244, 0.100619301199913, 0,
               0.852566361427307, -0.0525855533778667, 0.519966661930084, 0,
               0, 0, 0, 1),
             .Dim = c(4L, 4L)))
## Take a screen shot with Preview, save as fat3d-right.pdf
## rgl.snapshot("tmp.png") ## not good enough

rgl::par3d(userMatrix=structure(c(0.487184554338455, 0.133282631635666, 0.863068342208862,
0, 0, 0.988285005092621, -0.152619689702988, 0, -0.873299062252045,
0.0743539556860924, 0.481477200984955, 0, 0, 0, 0, 1), .Dim = c(4L,
4L)))
## Take a screen shot with Preview, save as fat3d-left.pdf

rgl::par3d(userMatrix=structure(c(0.983415901660919, 0.0405623130500317, 0.17677067220211,
0, 0, 0.974669396877289, -0.223650485277176, 0, -0.181364744901657,
0.219941437244415, 0.958505392074585, 0, 0, 0, 0, 1), .Dim = c(4L,
4L)))
## Take a screen shot with Preview, save as fat3d-center.pdf


###################################################
### code chunk number 3: regb.tex:190-195
###################################################
## hhcapture("ls2.Rout", '
fat2.lm <- lm(bodyfat ~ abdomin + biceps, data=fat)
anova(fat2.lm)
summary(fat2.lm)
## ')


###################################################
### code chunk number 4: regb.tex:273-276
###################################################
## hhpdf("f7.pdf", height=7, width=9, col=likertColor(2)[2:1]) ## col is not an argument for grDevices:::pdf
lmplot(fat2.lm)
## hhdev.off()


###################################################
### code chunk number 5: regb.tex:717-729
###################################################
## hhpdf("regb_f_hpr.pdf", height=7, width=7, col=likertColor(2)[2]) ## col is not an argument for grDevices:::pdf
data(houseprice)
houseprice$customf <- factor(houseprice$custom,
                             levels=c(0,1),
                             labels=c("regular","custom"))
houseprice$cornerf <- factor(houseprice$corner,
                             levels=c(0,1),
                             labels=c("middle","corner"))

splom(~houseprice[,c(1,2,5,6,7)], pch=16, cex=.5, main="houseprice",
      axis.text.cex=.7, pscales=3, xlab=NULL)
## hhdev.off()


###################################################
### code chunk number 6: regb.tex:745-759
###################################################
## hhpdf("regb_f_hprcc.pdf", height=7, width=10, col=likertColor(2)[2:1]) ## col is not an argument for grDevices:::pdf
splom(~houseprice[c(1,2,5)] | cornerf,
      group=customf,
      data=houseprice,
      layout=c(2,1),
      auto.key=list(space="top", border=TRUE, title="custom"),
      main="houseprice by custom | corner",
      par.strip.text=list(cex=1.5),
      subpanel.scales=list(cex=.8), pscales=4,
      axis.text.cex=.7,
      panel.cex=1.2,
      par.settings=list(superpose.symbol=list(pch=c(17,16))),
      xlab=NULL)
## hhdev.off()


###################################################
### code chunk number 7: regb.tex:785-791
###################################################
## hhcapture("houseprice3.Rout", '
houseprice.lm2 <- lm(price ~ sqft + taxes + custom + corner,
                     data=houseprice)
anova(houseprice.lm2)
summary(houseprice.lm2)
## ')


###################################################
### code chunk number 8: regb.tex:819-823
###################################################
## hhcapture("houseprice2.Rout", '
houseprice.lm1 <- lm(price ~ sqft + taxes, data=houseprice)
anova(houseprice.lm1, houseprice.lm2)
## ')


###################################################
### code chunk number 9: regb.tex:963-975
###################################################
## hhcapture("hardness-lm.Rout", '
data(hardness)

hardness.lin.lm  <- lm(hardness ~ density,
                       data=hardness)
anova(hardness.lin.lm)

hardness.quad.lm <- lm(hardness ~ density + I(density^2),
                       data=hardness)
anova(hardness.quad.lm)
coef(summary.lm(hardness.quad.lm))
## ')


###################################################
### code chunk number 10: regb.tex:1027-1035
###################################################
## hhpdf("hardness-ls.pdf", height=4, width=7, col=likertColor(2)[2]) ## col is not an argument for grDevices:::pdf
SQ <- regrresidplot(hardness$density, hardness$hardness, xlim=c(20, 85),
                    resid.plot="square")
QU <- regrresidplot(hardness$density, hardness$hardness, xlim=c(20, 85),
                    resid.plot="square", lm.object=hardness.quad.lm)
update(c(SQ, QU, layout=c(2, 1)), xlab="density", ylab="hardness",
       between=list(x=1), scales=list(alternating=FALSE))
## hhdev.off()


###################################################
### code chunk number 11: regb.tex:1050-1063
###################################################
## hhcapture("hardness-lm-orth.Rout", '
data(hardness)
hardness.lin.lm <- lm(hardness ~ density,
                      data=hardness)
anova(hardness.lin.lm)
coef(summary.lm(hardness.lin.lm))
h2 <- data.frame(density=hardness$density, poly(hardness$density, 2))
xyplot(X1 + X2 ~ density, data=h2)  ## graph not shown in book
hardness.quad.orth.lm <- lm(hardness ~ density + h2$X2,
                            data=hardness)
anova(hardness.quad.orth.lm)
coef(summary.lm(hardness.quad.orth.lm))
## ')


###################################################
### code chunk number 12: regb.tex:1115-1122
###################################################
## hhcapture("nointercept.Rout", '
data(fat)
## usual model with intercept
xy.int.lm <- lm(bodyfat ~ biceps, data=fat)
summary(xy.int.lm)
anova(xy.int.lm)
## ')


###################################################
### code chunk number 13: regb.tex:1142-1149
###################################################
## hhcapture("nointerceptB.Rout", '
data(fat)
## model without a constant term
xy.noint.lm <- lm(bodyfat ~ biceps - 1, data=fat)
summary(xy.noint.lm)
anova(xy.noint.lm)
## ')


###################################################
### code chunk number 14: regb.tex:1174-1192
###################################################
## hhpdf("nointercept.pdf", height=4.5, width=9)
A <-
xyplot(bodyfat ~ biceps, data=fat, pch=19, col=likertColor(2)[2],
       key=list(title="model",
                space="right",
                text=list(c("intercept", "no intercept")),
                lines=list(lty=c(1,2), lwd=2),
                col=c("black","red"),
                border=TRUE)) +
  layer(panel.abline(xy.int.lm, col="black", lty=1)) +
  layer(panel.abline(xy.noint.lm, col="red", lty=2))

B <-
update(A, xlim=c(-5,46), ylim=c(-22,35)) +
  layer(panel.abline(h=0, v=0, col="gray40", lty=2))

c(A, B, layout=c(2,1))
## hhdev.off()


###################################################
### code chunk number 15: regb.tex:1302-1315
###################################################
## hhcapture("fat2.Rout", '
fat2.lm <- lm(bodyfat ~ abdomin + biceps, data=fat)
pi.fit <- predict(fat2.lm,
                  newdata=data.frame(abdomin=93:94, biceps=33:34),
                  se.fit=TRUE, interval="prediction")

ci.fit <- predict(fat2.lm,
                  newdata=data.frame(abdomin=93:94,
                  biceps=33:34),
                  se.fit=TRUE, interval="confidence")
pi.fit
ci.fit$fit
## ')


###################################################
### code chunk number 16: regb.tex:1412-1416
###################################################
## hhpdf("longley.pdf", height=14, width=14)
data(longley)  ## from the datasets package
splom(longley, pch=19, cex=1.5, xlab=NULL, col=likertColor(2)[2])
## hhdev.off()


###################################################
### code chunk number 17: regb.tex:1451-1457
###################################################
## hhcapture("longley.Rout", '
longley.lm <- lm( Employed ~ . , data=longley)
summary(longley.lm)
anova(longley.lm)
vif(longley.lm)
## ')


###################################################
### code chunk number 18: regb.tex:1509-1527
###################################################
## rgl graphics doesn't automate.  It needs manual intervention for screen shots.
longley2.lm <- lm(Employed ~ Year + GNP, data=longley)
longley2.resid <- resid(longley2.lm)
car::scatter3d(Employed ~ Year + GNP, data=longley, fit="linear",
               residuals="squares",
               bg="white", axis.scales=TRUE, grid=TRUE, ellipsoid=FALSE,
               square.color = "gray80", surface.col="#a6cafe",
               surface.alpha=.3, sphere.size=.7,
               point.col=c("red","green")[1+(longley2.resid >= 0)])
rgl::par3d(userMatrix=
           structure(c(0.808914661407471, -0.262622833251953, 0.526009798049927,
                       0, 0, 0.894686996936798, 0.446693629026413, 0, -0.587926089763641,
                       -0.361337035894394, 0.723725438117981, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)))
rgl::par3d(windowRect=c(20,  576, 1220, 1964))
rgl::par3d(windowRect=c(20,  100, 500, 2000))
rgl::par3d(zoom=.3)
## screen capture as longley-collinear.pdf
## rgl.snapshot("tmp.png") ## not good enough


###################################################
### code chunk number 19: regb.tex:1693-1701
###################################################
## hhcapture("longley3.Rout", '
longley3.lm <- lm( Employed ~
         GNP.deflator + GNP + Unemployed + Armed.Forces + Year,
         data=longley)
summary(longley3.lm)
anova(longley3.lm)
vif(longley3.lm)
## ')


###################################################
### code chunk number 20: regb.tex:1719-1727
###################################################
## hhcapture("longley4.Rout", '
longley4.lm <- lm(Employed ~
                  GNP + Unemployed + Armed.Forces + Year,
                  data=longley)
summary(longley4.lm)
anova(longley4.lm)
vif(longley4.lm)
## ')


###################################################
### code chunk number 21: regb.tex:1751-1760
###################################################
## hhcapture("longley5.Rout", '
longley5.lm <- lm(Employed ~
                  Unemployed + Armed.Forces + Year,
                  data=longley)

summary(longley5.lm)
anova(longley5.lm)
vif(longley5.lm)
## ')


###################################################
### code chunk number 22: regb.tex:1983-1994
###################################################
## hhcapture("longley6.Rout", '
longley.subsets <-
  leaps::regsubsets(Employed ~ GNP.deflator + GNP +
                    Unemployed +
                    Armed.Forces + Population + Year,
                    data=longley, nbest=2)
longley.subsets.Summary <- summaryHH(longley.subsets)
## longley.subsets.Summary
tmp <- (longley.subsets.Summary$cp <= 10)
longley.subsets.Summary[tmp,]
## ')


###################################################
### code chunk number 23: regb.tex:2018-2022
###################################################
## hhpdf("regb-f4-longley.pdf", height=7, width=7)
plot(longley.subsets.Summary[tmp,], statistic='cp', legend=FALSE,
     ylim=c(3,7.5))
## hhdev.off()


###################################################
### code chunk number 24: regb.tex:2037-2041
###################################################
## hhcapture("longley7.Rout", '
longley.lm.7 <- lm.regsubsets(longley.subsets, 7)
summary(longley.lm.7)
## ')


###################################################
### code chunk number 25: regb.tex:2114-2120
###################################################
## hhpdf("longley-resid.pdf", height=9, width=14)
print(A4.left=.0125, panel.width=list(5,"cm"),
      residual.plots.lattice(longley.lm, par.strip.text=list(cex=1.1),
                             pch=19, col=likertColor(2)[2], cex=1.2)
)
## hhdev.off()


###################################################
### code chunk number 26: regb.tex:2367-2379
###################################################
## hhpdf("regb-fa-usair.pdf", height=7.5, width=7)
data(usair)
splom( ~ usair,
      main=expression("U.S. Air Pollution Data with" ~ SO[2] ~ "response variable"),
      xlab="Original Scaling",
      pch=19,                ## solid circles
      cex=.7,                ## size of points
      col=likertColor(2)[2], ## color of points
      pscales=3,             ## fewer tick labels
      axis.text.cex=.5,      ## smaller tick labels
      varname.cex=.7)        ## smaller variable name
## hhdev.off()


###################################################
### code chunk number 27: regb.tex:2390-2400
###################################################
## hhpdf("regb-fb-usair.pdf", height=7.5, width=7)
usair$lnSO2 <- log(usair$SO2)
usair$lnmfg <- log(usair$mfgfirms)
usair$lnpopn <- log(usair$popn)
splom( ~ usair[, c(8,2,9,10,5,6,7)],
              main=expression("U.S. Air Pollution Data with ln"*(SO[2])*" response variable"),
              xlab="Three log-transformed variables",
              pch=19, cex=.7, col=likertColor(2)[2],
              pscales=3, axis.text.cex=.5, varname.cex=.7)
## hhdev.off()


###################################################
### code chunk number 28: regb.tex:2440-2448
###################################################
## hhcapture("usair2.Rout", '
usair.regsubset <- leaps::regsubsets(
     lnSO2 ~ lnmfg + lnpopn + precip + raindays + temp + wind,
     data=usair, nbest=2)
usair.subsets.Summary <- summaryHH(usair.regsubset)
tmp <- (usair.subsets.Summary$cp <= 10)
usair.subsets.Summary[tmp,]
## ')


###################################################
### code chunk number 29: regb.tex:2463-2466
###################################################
## hhpdf("usair3.pdf", height=5.5, width=7)
plot(usair.subsets.Summary[tmp,], statistic='cp')
## hhdev.off()


###################################################
### code chunk number 30: regb.tex:2480-2486
###################################################
## hhcapture("usair5.Rout", '
usair.lm7 <- lm.regsubsets(usair.regsubset, 7)
anova(usair.lm7)
summary(usair.lm7)
vif(usair.lm7)
## ')



library(HH)


## rega.f4.s
data(fat)
par(mfrow=c(1,1))
regr2.plot(fat[,"abdomin"], xlab="abdomin",
           fat[,"biceps"],  ylab="biceps",
           fat[,"bodyfat"], zlab="bodyfat",
           resid.plot="square",
           eye=c(335.5, 115.65, 171.9),   ## used only in S-Plus
           theta=140, phi=35, r=sqrt(15), ## used only in R
           box=is.R(),
           plot.back.planes=FALSE,
           main="Least-squares with two X-variables")
## export.eps(hh("regb/figure/f4.eps"))


## ls2.s
fat2.lm <- lm(bodyfat ~ abdomin + biceps, data=fat)
summary(fat2.lm, corr=FALSE)
anova(fat2.lm)


## rega.f7.s
## follows ls2.s

old.par <- par(pch=16, mar=c(5,6,4,2)+.1)
par(mfrow=c(2,3))
plot(fat2.lm, cex=.8)
mtext("diagnostics for lm(bodyfat ~ abdomin + biceps, data=fat)",
      outer=TRUE, cex=1.4)
par(mfrow=c(1,1))
par(old.par)
## export.eps(hh("regb/figure/f7.eps"))


## splus.library/regr2.plot.s
## functions in the HH package

## splus.library/resid.squares.s
## functions in the HH package


## housepriceread.s
data(houseprice)
houseprice <- houseprice  ## local copy
houseprice$customf <- factor(houseprice$custom,
                             levels=c(0,1),
                             labels=c("regular","custom"))
houseprice$cornerf <- factor(houseprice$corner,
                             levels=c(0,1),
                             labels=c("middle","corner"))

## houseprice.f1.s
## houseprice.f1.r
if.R(s=
splom(~houseprice[,c(1,2,5,6,7)], pch=16, cex=.35, main="houseprice",
      superpanel=panel.pairs.hh, subpanel.scales=list(cex=.7), pscales=3,
      panel.cex=1)
,r=
splom(~houseprice[,c(1,2,5,6,7)], pch=16, cex=.35, main="houseprice",
      superpanel=panel.pairs, axis.text.cex=.7, pscales=3,
      panel.cex=1)
)
## export.eps(hh("regb/figure/regb.f.hpr.eps.gz"))


tpg <- trellis.par.get("superpose.symbol")
tpg <- lapply(tpg, function(x) x[1:2])
if.R(s={}, r={tpg$pch=c(15,16)})
if.R(s=
splom(~houseprice[c(1,2,5)],
      cex=.65,
      panel=panel.superpose,
      key=list(space="right",
        text=list(levels(houseprice$customf)),
        points=tpg,
        border=1,
        title="custom",
        cex=1),
      group=houseprice$customf,
      data=houseprice,
      main="houseprice by custom",
      superpanel=panel.pairs.hh, subpanel.scales=list(cex=.8), pscales=4,
      panel.cex=1.5)
,r=
splom(~houseprice[c(1,2,5)],
      cex=.65,
      panel=panel.superpose,
      key=list(space="right",
        text=list(levels(houseprice$customf)),
        points=tpg,
        border=1,
        title="custom",
        cex=1),
      group=houseprice$customf,
      data=houseprice,
      main="houseprice by custom",
      superpanel=panel.pairs, subpanel.scales=list(cex=.8), pscales=4,
      panel.cex=1.5,
      pch=c(15,16))
     )
## export.eps(hh("regb/figure/regb.f.hprc.eps.gz"))


tpg <- trellis.par.get("superpose.symbol")
tpg <- lapply(tpg, function(x) x[1:2])
if.R(s={}, r={tpg$pch=c(15,16)})
if.R(s=
splom(~houseprice[c(1,2,5)] | cornerf,
      panel=panel.superpose,
      layout=c(2,1),
      key=list(space="top",
        text=list(levels(houseprice$customf)),
        points=tpg,
        border=1,
        title="custom",
        cex=1),
      group=houseprice$customf,
      data=houseprice,
      main="houseprice by custom | corner",
      par.strip.text=list(cex=1.5),
      superpanel=panel.pairs.hh, subpanel.scales=list(cex=.8), pscales=4,
      panel.cex=1.2)
,r=
splom(~houseprice[c(1,2,5)] | cornerf,
      panel=panel.superpose,
      layout=c(2,1),
      key=list(space="top",
        text=list(levels(houseprice$customf)),
        points=tpg,
        border=1,
        title="custom",
        cex=1),
      group=houseprice$customf,
      data=houseprice,
      main="houseprice by custom | corner",
      par.strip.text=list(cex=1.5),
      superpanel=panel.pairs, subpanel.scales=list(cex=.8), pscales=4,
      axis.text.cex=.7,
      panel.cex=1.2,
      pch=c(15,16))
)
## export.eps(hh("regb/figure/regb.f.hprcc.eps.gz"))

## houseprice2.s
houseprice.lm1 <- lm(price ~ sqft + taxes, data=houseprice)
houseprice.lm2 <- lm(price ~ sqft + custom + corner + taxes,
                     data=houseprice)
anova(houseprice.lm1, houseprice.lm2)


## hardness.s
## hardness data
data(hardness)

## data
plot(hardness ~ density, data=hardness, pch=16)

## linear and quadratic regressions
hardness.lin.lm  <- lm(hardness ~ density,                data=hardness)
hardness.quad.lm <- lm(hardness ~ density + I(density^2), data=hardness)

anova(hardness.lin.lm)
anova(hardness.quad.lm)

summary.lm(hardness.lin.lm, corr=FALSE)
summary.lm(hardness.quad.lm, corr=FALSE)
coef(summary.lm(hardness.quad.lm, corr=FALSE))


## plots of both models
plot(hardness ~ density, data=hardness,
     main="linear and quadratic fit")
abline(hardness.lin.lm)
lines(x=hardness$density, y=predict(hardness.quad.lm))


plot(hardness ~ density, data=hardness, xlim=c(-10,70), ylim=c(-1000,3300),
     main="linear and quadratic fit with\nextended range to show 0 intercept for quadratic")
abline(hardness.lin.lm)
lines(x=hardness$density, y=predict(hardness.quad.lm))
x.density <- seq(-13,73, length=87)
lines(x=x.density,
      y=cbind(1, x.density, x.density^2) %*% coef(hardness.quad.lm),
      lty=4)
abline(h=0,v=0,lty=2)


## The squares of the least squares fit
## There are two options for display of the squares.
## 1. Skip both par() commands to make both plots full size
##    and then toggle between them on the screen.
## 2. Use both par() commands to put both plots on the same page so
##    they can looked at simultaneously.


par(mfrow=c(1,2))

plot(hardness ~ density, data=hardness,
     main="squared residuals for quadratic fit",
     xlim=c(20,83), ylim=c(0,3900))
hardness.hat <-
  cbind(1, hardness$density, hardness$density^2) %*% coef(hardness.quad.lm)
resid.squares(x=hardness$density,
              y=hardness$hardness,
              y.hat=hardness.hat)
lines(x=hardness$density, y=hardness.hat)

x.density <- seq(18,80, length=63)
y.hat <-
  cbind(1, x.density, x.density^2) %*% coef(hardness.quad.lm)
lines(x=x.density, y=y.hat)


regr1.plot(hardness$density, hardness$hardness,
           resid.plot="square",
           main="squared residuals for linear fit",
           xlab="density", ylab="hardness",
           points.yhat=FALSE,
           xlim=c(20,83), ylim=c(0,3900))

par(mfrow=c(1,1))
## export.eps(hh("regb/figure/hardness-ls.eps"))


## nointercept.s
## usual model with intercept
xy.int.lm <- lm(bodyfat ~ biceps, data=fat)
summary(xy.int.lm, corr=FALSE)

## model without a constant term
xy.noint.lm <- lm(bodyfat ~ biceps -1, data=fat)
summary(xy.noint.lm, corr=FALSE)


## plot the points and both lines
old.par <- par(cex=1.2, mar=par()$mar+c(0,.5,0,0), pch=16)
plot(bodyfat ~ biceps, data=fat)
abline(xy.int.lm, lty=1)
abline(xy.noint.lm, lty=4)
legend(c(38,45.5), c(8,15), c("xy.int.lm", "xy.noint.lm"), lty=c(1,4))
par(old.par)
## export.eps(hh("regb/figure/nointercept.eps"))


## splus.library/conf-lm.s
## functions in the HH package


## fat2.s
if.R(s=
     result <- predict(fat2.lm,
             newdata=data.frame(abdomin=93:94, biceps=33:34),
             se.fit=TRUE, pi=TRUE, ci=TRUE)
     ,r={
       result <- predict(fat2.lm,
                         newdata=data.frame(abdomin=93:94, biceps=33:34),
                         se.fit=TRUE, interval="prediction")
       result$pi.fit <- result$fit[,c(2,3)]
       result$fit <- result$fit[,1]
       result$ci.fit <- predict(fat2.lm,
                                newdata=data.frame(abdomin=93:94,
                                  biceps=33:34),
                                se.fit=TRUE, interval="confidence")$fit[,c(2,3)]
     })
result


## longley.s
## longley regression example.

## data is included with S-Plus and R

data(longley)

if.R(s=
     splom( ~ longley, pch=16, cex=.55,
           superpanel=panel.pairs.hh, subpanel.scales=list(cex=.8),
           pscales=2,
           panel.cex=.8)
     ,r=
     splom( ~ longley, pch=16,
           pscales=2,
           varname.cex=.8,
           axis.text.cex=.5)
   )
## export.eps(hh("regb/figure/longley.eps"))

longley.lm <- lm( Employed ~ . , data=longley)
summary(longley.lm, corr=FALSE)
anova(longley.lm)

vif( Employed ~ . , data=longley)


## partial residual plots
## standard S-Plus function
if.R(s={
  par(mfrow=c(2,3))
  partial.plot(longley.lm)
  par(mfrow=c(1,1))
},r=
     print("partial residual plots not included in R")
     )


if (FALSE) { ## prevent non-intentional running of this code.
## To make this illustration look good in the book,
## we are creating a custom-size square 11in by 11in postscript page.
## To see it in ghostview, set the media to "11x17" instead of
## the standard "letter" size.
##
psprint.800.800 <- function(filename) {
  trellis.device(postscript, horizontal=FALSE,
                 region=c(0,0,800,800), height=11, width=11,
                 file=filename)
}

print.four.rows <- function(tmp) {
  print(position=c(-.025,-.050, 1.025,  .250), more=TRUE, tmp[[4]])
  print(position=c(-.025, .200, 1.025,  .500), more=TRUE, tmp[[3]])
  print(position=c(-.025, .450, 1.025,  .750), more=TRUE, tmp[[2]])
  print(position=c(-.025, .700, 1.025, 1.000), more=FALSE, tmp[[1]])
  invisible(tmp)
}

## trellis.device.hh.color(orientation="portrait")

## HH constructed display
## free y-scale in partial residual and added variable plots

tmp.free <- residual.plots(longley.lm)
print.four.rows(tmp.free)
## The strip labels are badly placed on the S-Plus gui graphsheet.
##
## The strip labels are properly placed when sent directly to the
## postscript driver.
## psprint.800.800(hh("regb/figure/longley.resid.eps"))
## print.four.rows(tmp.free)
## dev.off()


## same y-scale in partial residual and added variable plots
  tmp.same <- residual.plots(longley.lm, y.relation="same")
  print.four.rows(tmp.same)
  ## The strip labels are badly placed on the S-Plus gui graphsheet.
  ##
  ## The strip labels are properly placed when sent directly to the
  ## postscript driver.
  ## psprint.800.800(hh("regb/figure/longley.resid.same.eps"))
  ## print.four.rows(tmp.same)
  ## dev.off()


print.two.rows <- function(tmp) {
  print(position=c(-.025, .050, 1.025,  .500), more=TRUE, tmp[[2]])
  print(position=c(-.025, .450, 1.025,  .900), more=FALSE, tmp[[1]])
  invisible(tmp)
}

print.two.rows(tmp.same[1:2])
print.two.rows(tmp.same[3:4])
## The strip labels are badly placed on the S-Plus gui graphsheet.
##
## The strip labels are properly placed when sent directly to the
## postscript driver.
## trellis.device(postscript, horizontal=TRUE,
##                file=hh("regb/figure/longley.resid.same.2yr.eps"))
## print.two.rows(tmp.same[1:2])
## dev.off()
##
## trellis.device(postscript, horizontal=TRUE,
##                file=hh("regb/figure/longley.resid.same.2pa.eps"))
## print.two.rows(tmp.same[3:4])
## dev.off()
}

## collinear.s
## following longley.s

regr2.plot(longley[,"Year"],  xlab="Year",
           longley[,"GNP"], ylab="GNP",
           longley[,"Employed"], zlab="Employed",
           resid.plot="square",
           theta=-40, phi=25, r=sqrt(3), ## used only in R
           box=is.R(),
           plot.back.planes=FALSE,
           plot.base.plane=FALSE,
           main="Least squares with two highly collinear X-variables")
## export.eps(hh("regb/figure/longley.collinear.eps"))


## longley.back.s
## hands-on and mechanical approach to thinning the full Longley model

## drop Population from complete model in longley.lm
longley3.lm <- lm( Employed ~
                  GNP.deflator + GNP + Unemployed + Armed.Forces + Year,
                  data=longley)
summary(longley3.lm, corr=FALSE)
anova(longley3.lm)
vif( Employed ~   GNP.deflator + GNP + Unemployed + Armed.Forces + Year,
    data=longley)


## drop GNP.deflator, highest p-value
longley4.lm <- lm(Employed ~
                                 GNP + Unemployed + Armed.Forces + Year,
                  data=longley)
summary(longley4.lm, corr=FALSE)
anova(longley4.lm)
vif( Employed ~                  GNP + Unemployed + Armed.Forces + Year,
    data=longley)



## drop GNP, large VIF and largest p-value
longley5.lm <- lm( Employed ~
                                       Unemployed + Armed.Forces + Year,
                  data=longley)

summary(longley5.lm, corr=FALSE)
anova(longley5.lm)
vif( Employed ~                        Unemployed + Armed.Forces + Year,
    data=longley)


## longley6.s
## longley6.r
## stepwise regression analysis of longley data

if.R(s={
longley.step <- stepwise(y=longley$Employed,
                         x=longley[,c(1:6)],
                         method="exhaustive",
                         plot=FALSE, nbest=2)
## longley.step  ## no need to print this.  longley.cp is more legible.

longley.cp <- cp.calc(longley.step, longley, "Employed")
tmp <- (longley.cp$cp <= 10)
longley.cp[tmp,]

old.par <- par(mar=par()$mar+c(0,1,0,0))
plot(cp ~ p, data=longley.cp[tmp,], ylim=c(0,10), type="n", cex=1.3)
abline(b=1)
text(x=longley.cp$p[tmp], y=longley.cp$cp[tmp],
     row.names(longley.cp)[tmp], cex=1.3)
title(main="Cp plot for longley.dat, Cp<10")
par(old.par)
## export.eps(hh("regb/figure/regb.f4.longley.eps"))
}
     ,r={
       longley.subsets <-
         leaps::regsubsets(Employed ~ GNP.deflator + GNP + Unemployed + Armed.Forces + Population + Year,
                    data=longley, nbest=2)
       longley.subsets.Summary <- summaryHH(longley.subsets)
       longley.subsets.Summary
       tmp <- (longley.subsets.Summary$cp <= 10)
       longley.subsets.Summary[tmp,]
       plot(longley.subsets.Summary[tmp,], statistic='cp', legend=FALSE)
       longley.lm.7 <- lm.regsubsets(longley.subsets, 7)  ## subset 7 has largest adjr2
       summary(longley.lm.7)
     })


## splus.library/residual.plots.s
## functions in the HH package


## splus.library/partial.corr.s
## functions in the HH package


## usair.s
## usair.r
data(usair)
usair <- usair  ## local copy
splom(~usair, main="U.S. Air Pollution Data with SO2 response", cex=.5)
## export.eps(hh("regb/figure/regb.f1.usair.eps"))

usair$lnSO2 <- log(usair$SO2)
usair$lnmfg <- log(usair$mfgfirms)
usair$lnpopn <- log(usair$popn)

usair[1:3,]   ## lnSO2 is in position 8, SO2 is in position 1
              ## lnmfg is in position 9, lnpopn is in position 10

splom(~usair[, c(8,2,9,10,5:7)],
              main="U.S. Air Pollution Data with 3 log-transformed variables",
              cex=.5)
## export.eps(hh("regb/figure/regb.f2.usair.eps"))

if.R(s={
  usair.step <- stepwise(y=usair$lnSO2,
                         x=usair[, c(2,9,10,5:7)],
                         method="exhaustive",
                         plot=FALSE, nbest=2)
  ## print for pedagogical purposes only.  The plot of cp ~ p is more useful.
  ## The line with rss=1e35 is a stepwise() bug, that we reported to S-Plus.
  print(usair.step, digits=4)
  usair.cp <- cp.calc(usair.step, usair, "lnSO2")
  ## print for pedagogical purposes only.  The plot of cp ~ p is more useful.
  usair.cp
  tmp <- (usair.cp$cp <= 10)
  usair.cp[tmp,]

  old.par <- par(mar=par()$mar+c(0,1,0,0))
  tmp <- (usair.cp$cp <= 10)
  plot(cp ~ p, data=usair.cp[tmp,], ylim=c(0,10), type="n", cex=1.3)
  abline(b=1)
  text(x=usair.cp$p[tmp], y=usair.cp$cp[tmp],
       row.names(usair.cp)[tmp], cex=1.3)
  title(main="Cp plot for usair.dat, Cp<10")
  par(old.par)
## export.eps(hh("regb/figure/regb.f3.usair.eps"))
},r={
  usair.regsubset <- leaps::regsubsets(lnSO2~lnmfg+lnpopn+precip+raindays+temp+wind, data=usair, nbest=2)
  usair.subsets.Summary <- summaryHH(usair.regsubset)
  tmp <- (usair.subsets.Summary$cp <= 10)
  usair.subsets.Summary[tmp,]
  plot(usair.subsets.Summary[tmp,], statistic='cp', legend=FALSE)

  usair.lm7 <- lm.regsubsets(usair.regsubset, 7)
  anova(usair.lm7)
  summary(usair.lm7)
})

vif(lnSO2 ~ temp + lnmfg + lnpopn + wind + precip + raindays, data=usair)

vif(lnSO2 ~ temp + lnmfg + wind + precip, data=usair)

usair.lm <- lm(lnSO2 ~ temp + lnmfg + wind + precip, data=usair)
anova(usair.lm)
summary(usair.lm, corr=FALSE)

## dsgn/code/cc176.s  ## Figure 13.1
## please see the code in file Ch13-dsgn.r

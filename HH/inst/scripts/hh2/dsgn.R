### R code from vignette source '~/WindowsC/HOME/rmh/hh.e2/hh2/dsgn.tex'

###################################################
### code chunk number 1: dsgn.tex:17-18
###################################################
library(HH)


###################################################
### code chunk number 2: dsgn.tex:21-26
###################################################
## the standard lattice color 2 is difficult for people with color deficient vision
data(col3x2)
## These colors look like a 3x2 color array when run through
## the vischeck simulator to see how they look for the three most
## common color vision deficiencies: Protanope, Deuteranope, Tritanope.


###################################################
### code chunk number 3: dsgn.tex:86-112
###################################################
## hhpdf("cc176-1.pdf", col=col3x2) ## col is not an argument for grDevices:::pdf
data(cc176)
## the positions need to be moved into the data object
useOuterStrips(
xyplot(wt.d ~ wt.n | n.treats*current, data=cc176,
       group=minutes,
       panel=function(x, y, ...) {
         panel.superpose(x, y, ...)
         panel.abline(lm(cc176$wt.d ~ cc176$wt.n), lty=3)
       },
       cex=1.6, pch=levels(cc176$minutes),
       par.strip.text=list(cex=1),
       strip=function(...) strip.default(strip.names=c(TRUE, TRUE), ...),
       scales=list(
         cex=1,
         x=list(alternating=1, at=seq(75, 200, 25)),
         y=list(alternating=2)),
       xlab=list("wt.n, weight of untreated other side", cex=1.2),
       ylab.right=list("wt.d, weight of treated muscle", cex=1.2),
       ylab=list("current", cex=1.2),
       xlab.top=list("n.treats, number of treatments", cex=1.2),
       between=list(x=1, y=1),
       sub=list("Plotting symbol is duration of the treatment in minutes", cex=.8)
       )
)
## hhdev.off()


###################################################
### code chunk number 4: dsgn.tex:154-175
###################################################
## hhcapture("cc176-1.Rout", '
## y=wt.d with x=wt.n as covariate
## (get essentially the same ANOVA as the approximate (y-bx)^2
## ANOVA table in Cochran and Cox)
cc176.aov <- aov(wt.d ~ rep + wt.n + n.treats*minutes*current,
                 data=cc176)
## summary(cc176.aov)
summary(cc176.aov,
        split=list(n.treats=list(n.treats.lin=1,
                                 n.treats.quad=2)),
        expand.split=FALSE)
##
## adjust y for x
cc176$y.adj <- cc176$wt.d  -
  (cc176$wt.n - mean(cc176$wt.n))*coef(cc176.aov)["wt.n"]
## duplicate CC Table 5.17
cc176.means <- tapply(cc176$y.adj,
                      cc176[,c("current","n.treats")], mean)
cc176.means
apply(cc176.means, 1, mean)
## ')


###################################################
### code chunk number 5: dsgn.tex:193-198
###################################################
## hhpdf("cc176-2.pdf", col=col3x2) ## col is not an argument for grDevices:::pdf
interaction2wt(y.adj ~ current + n.treats, data=cc176,
               main.cex=1.6,
               scales=list(x=list(cex=.7), y=list(cex=.9, alternating=FALSE)))
## hhdev.off()


###################################################
### code chunk number 6: dsgn.tex:233-252
###################################################
tmp <-
sapply(split(cc176$y.adj, cc176$current),
       function(x)
         c(min=min(x),
           "m-sd"=mean(x)-sd(x),
           mean=mean(x),
           "m+sd"=mean(x)+sd(x),
           max=max(x)))

t(tmp)[4:1,]
## hhpdf("cc176-microplot.pdf", height=1.7, width=3, col=col3x2) ## col is not an argument for grDevices:::pdf
cc176.bwplot <-
bwplot(unpositioned(current) ~ y.adj, data=cc176,
       panel=panel.bwplot.intermediate.hh,
       xlab=NULL)
cc176.bwplot$par.settings$axis.line$col <- 0 ## $
update(cc176.bwplot, par.settings=list(clip=list(panel=FALSE)), scales=list(x=list(cex=.7))) +
   layer(panel.axis("bottom", line.col="black", text.col=0, outside=TRUE))
## hhdev.off()


###################################################
### code chunk number 7: dsgn.tex:304-321
###################################################
## hhpdf("cc176-4.pdf", col=col3x2) ## col is not an argument for grDevices:::pdf
useOuterStrips(
xyplot(y.adj ~ as.position(minutes) | n.treats + current, data=cc176,
       panel=panel.bwplot.superpose, groups=minutes,
       horizontal=FALSE, xlim=c(.5, 5.5),
       strip=function(...) strip.default(strip.names=c(TRUE, TRUE), ...),
       par.strip.text=list(cex=1.4),
       scales=list(
         y=list(alternating=2),
         x=list(alternating=1)),
       xlab=list("minutes", cex=1.2),
       ylab.right=list("response wt.d adjusted for the covariate wt.n", cex=1.2),
       ylab=list("current", cex=1.2),
       xlab.top=list("n.treats, number of treatments", cex=1.2),
       between=list(x=1, y=1))
)
## hhdev.off()


###################################################
### code chunk number 8: dsgn.tex:356-413
###################################################
## hhpdf("cc176-5.pdf", height=14, width=14, col=col3x2) ## col is not an argument for grDevices:::pdf
##
cc176.5b <- ## identical lines
useOuterStrips(
xyplot(wt.d ~ wt.n | n.treats*current, data=cc176,
       group=n.treats, pch=19, ##pch=levels(cc176$n.treats),
       panel=function(x, y, ...) {
         panel.superpose(x, y, ...)
         panel.abline(lm(cc176$wt.d ~ cc176$wt.n), lty=1, lwd=2)
       },
       cex=1.4,
       par.strip.text=list(cex=1),
       strip=function(...) strip.default(strip.names=c(TRUE, TRUE), ...),
       scales=list(
         cex=1,
         x=list(alternating=1, at=seq(75, 200, 25)),
         y=list(alternating=2)),
       xlab=list("wt.n, weight of untreated other side", cex=1.2),
       ylab.right=list("wt.d, weight of treated muscle", cex=1.2),
       ylab=list("current", cex=1.2),
       xlab.top=list("n.treats, number of treatments", cex=1.2),
       between=list(x=1, y=1),
       main="b. Identical lines: wt.d ~ wt.n"
       )
)
##
cc176.5a <- update(cc176.5b, ## horizontal lines by group
                   panel=
                   function(x, y, ...) {
                     panel.superpose(x, y, ...)
                     panel.abline(a=mean(y), b=0, lty=1, lwd=2)
                   },
                   main="a. Separate horizontal lines: wt.d ~ n.c")
##
abline.args <- coef(aov(wt.d ~ -1 + wt.n + interaction(n.treats, current), data=cc176))
cc176.5cd <- update(cc176.5b, ## parallel lines by group
                   panel=
                   function(x, y, ...) {
                     panel.superpose(x, y, ...)
                     panel.abline(a=abline.args[panel.number()+1],
                                  b=abline.args[1], lty=1, lwd=2)
                   },
                   main="c,d. Parallel lines: wt.d ~ n.c + wt.n")
##
cc176.5e <- update(cc176.5b, ## separate lines
                   panel=
                   function(x, y, ...) {
                     panel.superpose(x, y, ...)
                     panel.abline(lm(y ~ x), lty=1, lwd=2)
                   },
                   main="e. Separate lines: wt.d ~ wt.n * n.c")
##
print(cc176.5a,  more=TRUE,  position=c(.00, .52,  .48, 1.0)) ## split=c(1,1,2,2))
print(cc176.5b,  more=TRUE,  position=c(.52, .52, 1.00, 1.0)) ## split=c(2,1,2,2))
print(cc176.5cd, more=TRUE,  position=c(.00, .00,  .48,  .48)) ## split=c(1,2,2,2))
print(cc176.5e,  more=FALSE, position=c(.52, .00, 1.00,  .48)) ## split=c(2,2,2,2))
## hhdev.off()


###################################################
### code chunk number 9: dsgn.tex:492-503
###################################################
## hhcapture("cc176-6.Rout", '
  cc176t <- cc176
  for (i in names(cc176t))
    if (is.factor(cc176t[[i]]))
      contrasts(cc176t[[i]]) <-
        contr.treatment(length(levels(cc176t[[i]])))
  sapply(cc176t, class)
  cc176t.aov <- aov(wt.d ~ rep + wt.n + n.treats + wt.n*current,
                    data=cc176t)
  summary(cc176t.aov)
## ')


###################################################
### code chunk number 10: dsgn.tex:517-546
###################################################
##   hhpdf("cc176-7.pdf", width=8)
  cc176.mmc <- mmc(cc176t.aov, focus="current")
  print(cc176.mmc)
  mmcplot(cc176.mmc, xlim=c(-21, 18), style="both")
##   hhdev.off()
  ##
##   hhpdf("cc176-7mmc.pdf", height=5, width=9)
  cc176mmc <- mmcplot(cc176.mmc, xlim=c(-21, 18), axis.right=2.2)
  cc176mmc
##   hhdev.off()
  ##
##   hhpdf("cc176-7contr.pdf")
  current.lmat <- cbind("cc-gf"=c(-1,-1, 1, 1),
                        "25-60"=c( 0, 0,-1, 1),
                        "g-f"  =c( 1,-1, 0, 0))
  dimnames(current.lmat)[[1]] <- levels(cc176$current)
  cc176.mmc <- mmc(cc176t.aov, focus="current", focus.lmat=current.lmat)
  print(cc176.mmc)
  mmcplot(cc176.mmc, xlim=c(-21, 18), type="lmat", style="both")
##   hhdev.off()
  ##
##   hhpdf("cc176-7mmccontr.pdf", height=5, width=9)
  cc176mmccontr <- mmcplot(cc176.mmc, xlim=c(-21, 18), type="lmat", axis.right=2.2)
  cc176mmccontr
##   hhdev.off()
##   hhpdf("cc176-7mmcandcontr.pdf", height=8.5, width=8)
  update(between=list(y=1), scales=list(alternating=1),
    c("Orthogonal Contrasts"=cc176mmccontr, MMC=cc176mmc, layout=c(1,2)))
##   hhdev.off()


###################################################
### code chunk number 11: dsgn.tex:709-726
###################################################
  data(tires)
  ## simpler
  bwplot(wear ~ car + position + brand, outer=TRUE, data=tires,
         between=list(x=1), layout=c(3,1))

  ## control of colors
##   hhpdf("tiresbw.pdf", width=6, height=2.5, col=col3x2) ## col is not an argument for grDevices:::pdf
  tpgcol <- trellis.par.get()$superpose.symbol$col
  tmp <-
    c(car=bwplot(wear ~ car, data=tires,
                 panel=panel.bwplot.superpose, groups=car, col=tpgcol[1]),
      position=bwplot(wear ~ position, data=tires,
                      panel=panel.bwplot.superpose, groups=position, col=tpgcol[2]),
      tires=bwplot(wear ~ brand, data=tires,
                   panel=panel.bwplot.superpose, groups=brand, col=tpgcol[3]))
   update(tmp, between=list(x=1), layout=c(3,1), par.strip.text=list(cex=1.2))
##   hhdev.off()


###################################################
### code chunk number 12: dsgn.tex:792-800
###################################################
##   hhcapture("tires.Rout", '
  data(tires)
  tires.aov <- aov(wear ~ car + position + brand, data=tires)
  summary(tires.aov)
  tapply(tires$wear, tires$car, "mean")
  tapply(tires$wear, tires$position, "mean")
  tapply(tires$wear, tires$brand, "mean")
## ')


###################################################
### code chunk number 13: dsgn.tex:816-832
###################################################
old.stars <- options(show.signif.stars=FALSE)
##   hhcapture("tires2.Rout", '
  tires.mmc.brand <- mmc(tires.aov, linfct=mcp(brand="Tukey"))
  ## print(tires.mmc.brand)
  brand.lmat <- cbind("1-43" =c( 2, 0,-1,-1),
                      "4-3"  =c( 0, 0,-1, 1),
                      "143-2"=c( 1,-3, 1, 1))
  dimnames(brand.lmat)[[1]] <- levels(tires$brand)
  tires.mmc.brand <- mmc(tires.aov, linfct=mcp(brand="Tukey"),
                         focus.lmat=brand.lmat)
  print(tires.mmc.brand)
  contrasts(tires$brand) <- brand.lmat
  tires.aov <- aov(wear ~ car + position + brand, data=tires)
  summary(tires.aov, split=list(brand=list("1-43"=1, rest=2:3)))
## ')
options(old.stars)


###################################################
### code chunk number 14: dsgn.tex:850-862
###################################################
##   hhpdf("tiresmmc.pdf", width=6, height=4)
  tiresmmc <- mmcplot(tires.mmc.brand, ylim=c(10.4, 14.4))
  tiresmmc
##   hhdev.off()
##   hhpdf("tiresmmccontr.pdf", width=6, height=4)
  tiresmmccontr <- mmcplot(tires.mmc.brand, type="lmat", ylim=c(10.4, 14.4))
  tiresmmccontr
##   hhdev.off()
##   hhpdf("tiresmmcandcontr.pdf", width=8.5, height=8)
  update(between=list(y=1), scales=list(alternating=1),
    c("Orthogonal Contrasts"=tiresmmccontr, MMC=tiresmmc, layout=c(1,2)))
##   hhdev.off()


###################################################
### code chunk number 15: dsgn.tex:981-987
###################################################
data(filmcoat)
## display data in table dsgntwo.t.filmcoat
print(quote=FALSE,
      reshape2::acast(temprt ~ pressure, value.var="coat",
                      fun=paste, collapse=",", data=filmcoat)
      )


###################################################
### code chunk number 16: dsgn.tex:1010-1021
###################################################
##   hhpdf("filmcoat.pdf", width=3.5, height=2.8, col=likertColor(2)[2]) ## col is not an argument for grDevices:::pdf
  useOuterStrips(
    bwplot(coat ~ 1 | pressure*temprt, data=filmcoat, horizontal=FALSE,
           xlab=NULL, xlab.top="Pressure",
           ylab="Temperature", ylab.right="coat",
           scales=list(x=list(draw=FALSE), y=list(alternating=2)),
           as.table=TRUE, aspect=1,
           par.settings=list(box.dot=list(
              col=trellis.par.get()$superpose.symbol$col[1])))
  )
##   hhdev.off()


###################################################
### code chunk number 17: dsgn.tex:1046-1052
###################################################
##   hhcapture("filmcoat.ma.Rout", '
  reshape2::acast(filmcoat, temprt ~ pressure, mean,
                  value.var="coat", margins=TRUE)
  film.aov1 <- aov(coat ~ temprt*pressure, data=filmcoat)
  summary(film.aov1)
## ')


###################################################
### code chunk number 18: dsgn.tex:1067-1072
###################################################
## hhpdf("filmcoatIntr.pdf", width=7, height=6, col=col3x2) ## col is not an argument for grDevices:::pdf
interaction2wt(data=filmcoat, coat ~ temprt+pressure,
               simple=TRUE, simple.scale=list(temprt=.3, pressure=.3),
               xlim=c(.5, 3.5), between=list(x=.5, y=.5))
## hhdev.off()


###################################################
### code chunk number 19: dsgn.tex:1093-1122
###################################################
## hhpdf("mcout5.pdf", height=6, width=12)
ResidMS <- function(x) summary(x)[[1]]["Residuals","Mean Sq"]
ResidMSAvg <- ResidMS(film.aov1)
crit.val <- qtukey(.95, 3, 18, 3)/sqrt(2)
## separate ANOVA for each pressure
filmcoat.aov.3p <- sapply(levels(filmcoat$pressure),
                          function(i) aov(coat ~ temprt,
                                          data=filmcoat,
                                          subset=(pressure==i)),
                          simplify=FALSE)
print(lapply(filmcoat.aov.3p, anova))
filmcoat.mmc.3pc <-
  sapply(filmcoat.aov.3p, simplify = FALSE,
         function(x) mmc(x, calpha=crit.val * sqrt(ResidMSAvg/ResidMS(x))))
## filmcoat.mmc.3pc
mmc3pb <- sapply(filmcoat.mmc.3pc, mmcplot, style="both", simplify=FALSE,
                 axis.right=2, xlim=c(-14, 14), ylim=c(33, 46),
                 ylab.right=NULL, ylab=NULL,
                 col.iso="gray40", col.contr0='gray40')
## mmc3pb
mmc3pb[[1]]$condlevels[[1]][1] <- names(mmc3pb)[1]
mmc3pb[[2]]$condlevels[[1]][1] <- names(mmc3pb)[2]
mmc3pb[[3]]$condlevels[[1]][1] <- names(mmc3pb)[3]
old.digits <- options(digits=4)
print(mmc3pb[[1]], position=c(0/3-.00,0,0/3-.00+.26,1), more=TRUE)
print(mmc3pb[[2]], position=c(1/3-.02,0,1/3-.02+.26,1), more=TRUE)
print(mmc3pb[[3]], position=c(2/3-.04,0,2/3-.04+.26,1), more=FALSE)
options(old.digits)
## hhdev.off()


###################################################
### code chunk number 20: dsgn.tex:1140-1166
###################################################
## hhpdf("mcout6.pdf", height=6, width=12)
## separate ANOVA for each temperature
filmcoat.aov.3t <- sapply(levels(filmcoat$temprt),
                          function(i) aov(coat ~ pressure,
                                          data=filmcoat,
                                          subset=(temprt==i)),
                          simplify=FALSE)
print(lapply(filmcoat.aov.3t, anova))
filmcoat.mmc.3tc <-
  sapply(filmcoat.aov.3t, simplify = FALSE,
         function(x) mmc(x, calpha=crit.val * sqrt(ResidMSAvg/ResidMS(x))))
## filmcoat.mmc.3tc
mmc3tb <- sapply(filmcoat.mmc.3tc, mmcplot, style="both", simplify=FALSE,
                 axis.right=2, xlim=c(-14, 14), ylim=c(33, 46),
                 ylab.right=NULL, ylab=NULL,
                 col.iso="gray40", col.contr0='gray40')
## mmc3tb
mmc3tb[[1]]$condlevels[[1]][1] <- names(mmc3tb)[1]
mmc3tb[[2]]$condlevels[[1]][1] <- names(mmc3tb)[2]
mmc3tb[[3]]$condlevels[[1]][1] <- names(mmc3tb)[3]
old.digits <- options(digits=4)
print(mmc3tb[[1]], position=c(0/3-.00,0,0/3-.00+.26,1), more=TRUE)
print(mmc3tb[[2]], position=c(1/3-.02,0,1/3-.02+.26,1), more=TRUE)
print(mmc3tb[[3]], position=c(2/3-.04,0,2/3-.04+.26,1), more=FALSE)
options(old.digits)
## hhdev.off()


###################################################
### code chunk number 21: dsgn.tex:1310-1325
###################################################
## hhpdf("dsgnGunload.pdf", col=likertColor(2)[2]) ## col is not an argument for grDevices:::pdf
data(gunload)
useOuterStrips(
bwplot(rounds  ~ team | method*group, data=gunload,
       ylab=list("group", cex=1.4),
       ylab.right=list("rounds", cex=1.4),
       xlab=list("team %in% group", cex=1.4),
       xlab.top=list("method", cex=1.4),
       par.strip.text=list(cex=1.4),
       scales=list(x=list(cex=1), y=list(cex=1.2, alternating=2)),
       par.settings=list(box.dot=list(
          col=trellis.par.get()$superpose.symbol$col[1]))
       )
)
## hhdev.off()


###################################################
### code chunk number 22: dsgn.tex:1537-1545
###################################################
##   hhcapture("gunloada.Rout", '
gunload.aov <-
   aov(rounds ~ method*group + Error((team %in% group)/method),
       data=gunload)
## The R warning is ok.  There are no treatment terms inside the Error stratum.
summary(gunload.aov)
model.tables(gunload.aov, type="means")
## ')


###################################################
### code chunk number 23: dsgn.tex:1663-1681
###################################################
##   hhcapture("turkeyFactors.Rout", '
data(turkey)

turkey[c(1,7,13,19,25),]

turkey$trt.vs.control <-
   factor(rep(c("control","treatment"), c(6,24)))
contrasts(turkey$trt.vs.control) <- c(4,-1)

turkey$additive <- factor(rep(c("control","A","B"), c(6,12,12)),
                          levels=c("control","A","B"))
contrasts(turkey$additive) <- c(0,1,-1)

turkey$amount <- factor(rep(c(0,1,2,1,2), c(6,6,6,6,6)))
contrasts(turkey$amount) <- c(0,1,-1)

turkey[c(1,7,13,19,25),]
## ')


###################################################
### code chunk number 24: dsgn.tex:1694-1699
###################################################
##   hhcapture("turkeyAov2.Rout", '
turkey3.aov <- aov(wt.gain ~ trt.vs.control / (additive*amount),
                   data=turkey, x=TRUE)
summary(turkey3.aov)
## ')


###################################################
### code chunk number 25: dsgn.tex:1715-1722
###################################################
##   hhcapture("turkeyMeans.Rout", '
print(na.print="",
tapply(turkey$wt.gain,
       turkey[,c("additive","amount")],
       mean)
)
## ')


###################################################
### code chunk number 26: dsgn.tex:1736-1750
###################################################
## hhpdf("turkeyF2V.pdf", height=4, width=4, col=likertColor(2)[2]) ## col is not an argument for grDevices:::pdf
additive.rev <- ordered(turkey$additive, rev(levels(turkey$additive)))
useOuterStrips(
bwplot( wt.gain ~ rep(1, 30) | amount * additive.rev,
       data=turkey,
       horizontal=FALSE, aspect=1,
       scales=list(x=list(draw=FALSE), y=list(alternating=2)),
       xlab=NULL, xlab.top=list("amount", cex=1.4),
       ylab=list("additive", cex=1.4),
       ylab.right=list("wt.gain", cex=1.4),
       par.strip.text=list(cex=1.4),
       groups=rep(1,30), panel=panel.bwplot.superpose)
)
## hhdev.off()


###################################################
### code chunk number 27: dsgn.tex:1808-1809
###################################################
data(abc)


###################################################
### code chunk number 28: dsgn.tex:1856-1873
###################################################
##   hhcapture("abcrearrange1.Rout", '
data(abc)

abc
abc.oneway <- ## one-way
with(abc,
     matrix(y, 4, 3, dimnames=list(1:4, A=unique(A)))
     )
abc.oneway

abc.crossed <- ## crossed
with(abc,
     matrix(y, 3, 4, byrow=TRUE,
            dimnames=list(A=unique(A), B=unique(B)))
     )
abc.crossed
## ')


###################################################
### code chunk number 29: dsgn.tex:1887-1902
###################################################
##   hhcapture("abcrearrange2.Rout", '
abc.nested <- ## nested
with(abc,
     matrix(c(y[1:4],    rep(NA,8),
              rep(NA,4), y[5:8],    rep(NA,4),
              rep(NA,8),            y[9:12]),
              3, 12, byrow=TRUE,
              dimnames=list(A=unique(A), BwA=BwA))
            )
print(abc.nested, na.print="")

abc.double.indexed <- ## doubly-indexed
  abc[,"y",drop=FALSE]
abc.double.indexed
## ')


###################################################
### code chunk number 30: dsgn.tex:1922-1966
###################################################
##  hhcode("modelsAB.R", '
## one-way
abc.A.aov <- aov(y ~ A, data=abc)
anova(abc.A.aov)
coef(abc.A.aov)
contrasts(abc$A)
model.matrix(abc.A.aov)


## crossed: no interaction
abc.ApB.aov <- aov(y ~ A+B, data=abc)
anova(abc.ApB.aov)
coef(abc.ApB.aov)
contrasts(abc$A)
contrasts(abc$B)
model.matrix(abc.ApB.aov)


## crossed: with interaction
abc.AsB.aov <- aov(y ~ A*B, data=abc)
anova(abc.AsB.aov)
coef(abc.AsB.aov)
contrasts(abc$A)
contrasts(abc$B)
contrasts(abc$AB)
model.matrix(abc.AsB.aov)


## nested
abc.BwA.aov <- aov(y ~ A/B, data=abc)
anova(abc.BwA.aov)
coef(abc.BwA.aov)
contrasts(abc$A)
contrasts(interaction(abc$A, abc$B))
model.matrix(abc.BwA.aov)


## doubly-indexed
abc.AB.aov <- aov(y ~ AB, data=abc)
anova(abc.AB.aov)
coef(abc.AB.aov)
contrasts(abc$AB)
model.matrix(abc.AB.aov)
## ')


###################################################
### code chunk number 31: dsgn.tex:1984-1994
###################################################
##   hhcapture("contrasts-contr1.Rout", '
model.matrix(~A, data=abc,
   contrasts=
      list(A=contr.treatment))
## ')
##   hhcapture("contrasts-contr2.Rout", '
model.matrix(~A, data=abc,
   contrasts=
      list(A=contr.sum))
## ')


###################################################
### code chunk number 32: dsgn.tex:2034-2049
###################################################
##   hhcapture("contrasts-AB.Rout", '
old.width <- options(width=70)
mm <- model.matrix(~A*B, data=abc,
             contrasts=list(A=contr.sum, B=contr.sum))
mm[,]

print(AA <- mm["r.z", c("A1","A2")])
print(BBB <- mm["r.z", c("B1","B2","B3")])

outer(AA, BBB)
as.vector(outer(AA, BBB))
mm["r.z", c("A1:B1","A2:B1","A1:B2","A2:B2","A1:B3","A2:B3")]

options(old.width)
## ')


###################################################
### code chunk number 33: dsgn.tex:2090-2104
###################################################
##   hhcapture("turkeyAov3.Rout", '
match(dimnames(coef(summary.lm(turkey3.aov)))[[1]],
      dimnames(turkey3.aov$x)[[2]])
turkey[c(1,7,13,19,25),]
turkey3.coef <- summary.lm(turkey3.aov)$coef
turkey3.x <- turkey3.aov$x
term.names <-
    c("(Intercept)","trt.vs.control","additive","amount",
      "additive:amount")
dimnames(turkey3.coef)[[1]] <- term.names
dimnames(turkey3.x)[[2]][c(1,2,4,8,12)] <- term.names
zapsmall(turkey3.coef)
turkey3.x[c(1,7,13,19,25), c(1,2,4,8,12)]
## ')


###################################################
### code chunk number 34: dsgn.tex:2304-2341
###################################################
##   hhcapture("cloverT123sc.Rout", '
data(rhiz.clover)
## drop two observation to illustrate Type II and III sums of squares
## I am dropping the non-outlier observations in 3D0k5
cloverD <- rhiz.clover[-c(7,9,10),]

old.opt <- options(show.signif.stars=FALSE, digits=3)

cloverDsc.aov <-
  aov(Npg ~ strain * comb,
      data=cloverD,
      contrasts=
        list(strain=contr.sum,
             comb=contr.sum))

anova(cloverDsc.aov)[,c(2,1,4,5)]

car::Anova(cloverDsc.aov, type=2)

car::Anova(cloverDsc.aov, type=3)
## ')
##   hhcapture("cloverT123cs.Rout", '
cloverDcs.aov <-
  aov(Npg ~ comb * strain,
      data=cloverD,
      contrasts=
        list(strain=contr.sum,
             comb=contr.sum))

anova(cloverDcs.aov)[,c(2,1,4,5)]

car::Anova(cloverDcs.aov, type=2)

car::Anova(cloverDcs.aov, type=3)

options(old.opt)
## ')


###################################################
### code chunk number 35: dsgn.tex:2387-2404
###################################################
## hhpdf("cloverD.pdf", height=6, width=4)

print(position = c(0, .50, 1, 1.00), more = TRUE,  # top
dotplot(Npg ~ strain | comb, data=rhiz.clover,
        main=list("clover: Nitrogen per Gram --- full data", cex=.8),
        layout=c(2,1), col=likertColor(2)[2],
        scales=list(cex=.75, x=list(rot=90)))
)

print(position = c(0,  0, 1, .50), more = FALSE,  # bottom
dotplot(Npg ~ strain | comb, data=cloverD,
        main=list("clover: Nitrogen per Gram --- Observations 7,9,10 dropped", cex=.8),
        layout=c(2,1), col=likertColor(2)[2],
        scales=list(cex=.75, x=list(rot=90)))
)

## hhdev.off()


###################################################
### code chunk number 36: dsgn.tex:2421-2434
###################################################
##   hhcapture("fatLst.Rout", '

data(fat)
fat.lm <- lm(bodyfat ~ abdomin + biceps, data=fat)
## regression coefficients
coef(summary(fat.lm))
## sequential sums of squares (Type I)
anova(fat.lm)
## weighted squares of means (Type III)
car::Anova(fat.lm, type="III")
## model sum of squares
var(fat$bodyfat) * (nrow(fat)-1) - sum(fat.lm$residuals^2)
## ')


###################################################
### code chunk number 37: dsgn.tex:2587-2604
###################################################
## hhpdf("vulcanInteraction.pdf", height=6, width=6, col=col3x2) ## col is not an argument for grDevices:::pd)
data(vulcan)
levels(vulcan$raw)
vulcan$raw <- ordered(vulcan$raw, levels=c(4, 1, 2, 3))
levels(vulcan$raw)
position(vulcan$raw) <- c(1.1, 2.3, 3.5, 4.8)

levels(vulcan$filler)
vulcan$filler <- ordered(vulcan$filler, levels=c(1, 3, 5, 2, 4))
levels(vulcan$filler)

interaction2wt(wear ~ filler+raw+pretreat, data=vulcan,
               main.cex=1.6,
               par.strip.text=list(cex=.6),
               xlim=c(.3, 5.6))

## hhdev.off()


###################################################
### code chunk number 38: dsgn.tex:2617-2638
###################################################
## The two-factor 2-way interaction figure is not in the book.
## It is here for comparison with the simple effects figure.
## hhpdf("vulcan2factor.pdf", height=6, width=6, col=col3x2) ## col is not an argument for grDevices:::pdf
interaction2wt(wear ~ filler+raw, data=vulcan,
               main.cex=1.6,
               par.strip.text=list(cex=.8),
               simple.scale=list(filler=.18, raw=.15),
               box.ratio=.05,
               xlim=c(.3, 5.6))
## hhdev.off()

## hhpdf("vulcanSimple.pdf", height=6, width=6, col=col3x2) ## col is not an argument for grDevices:::pdf

interaction2wt(wear ~ filler+raw, data=vulcan, simple=TRUE,
               main.cex=1.6,
               par.strip.text=list(cex=.8),
               simple.scale=list(filler=.18, raw=.15),
               box.ratio=.05,
               xlim=c(.3, 5.6))

## hhdev.off()


###################################################
### code chunk number 39: dsgn.tex:2681-2700
###################################################
##  hhcode("tiresLatin.R", '
## R defaults to treatment contrasts for factors.
## We need an orthogonal set of factors for this exercise.
##
data(tires)
contrasts(tires$car)      <- contr.helmert(4)
contrasts(tires$position) <- contr.helmert(4)
contrasts(tires$brand)    <- contr.helmert(4)

tires.aov <- aov(wear ~ car + position + brand, data=tires, x=TRUE)
anova(tires.aov)
tires.rc.aov <- aov(wear ~ car * position, data=tires, x=TRUE)
anova(tires.rc.aov)

t(tires.aov$x[,8:10])
t(tires.rc.aov$x[,8:16])
tr1.lm <- lm(tires.aov$x[,8] ~ tires.rc.aov$x[,8:16])
anova(tr1.lm)
## ')


###################################################
### code chunk number 40: dsgn.tex:2809-2822
###################################################
##  hhcode("turkeyAov3Match.R", '
                    summary.lm(turkey3.aov)
               coef(summary.lm(turkey3.aov))
      dimnames(coef(summary.lm(turkey3.aov)))
      dimnames(coef(summary.lm(turkey3.aov)))[[1]]

               turkey3.aov$x
      dimnames(turkey3.aov$x)
      dimnames(turkey3.aov$x)[[2]]

match(dimnames(coef(summary.lm(turkey3.aov)))[[1]],
      dimnames(turkey3.aov$x)[[2]])
## ')


###################################################
### code chunk number 41: dsgn.tex:2884-2897
###################################################
## hhpdf("turkeyF2H.pdf", height=4, width=4, col=likertColor(2)[2]) ## col is not an argument for grDevices:::pdf
useOuterStrips(
bwplot( ~ wt.gain | amount * additive.rev,
       data=turkey,
       aspect=1,
       scales=list(x=list(alternating=1)),
       xlab=list(cex=1.4),
       xlab.top=list("amount", cex=1.4),
       ylab=list("additive", cex=1.4),
       par.strip.text=list(cex=1.4),
       groups=rep(1,30), panel=panel.bwplot.superpose)
)
## hhdev.off()



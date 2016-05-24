### R code from vignette source '~/WindowsC/HOME/rmh/hh.e2/hh2/regc.tex'

###################################################
### code chunk number 1: regc.tex:9-10
###################################################
library(HH)


###################################################
### code chunk number 2: regc.tex:13-18
###################################################
## the standard lattice color 2 is difficult for people with color deficient vision
data(col3x2)
## These colors look like a 3x2 color array when run through
## the vischeck simulator to see how they look for the three most
## common color vision deficiencies: Protanope, Deuteranope, Tritanope.


###################################################
### code chunk number 3: regc.tex:97-105
###################################################
## hhpdf("rent1.pdf", width=9, height=5.5)
data(rent)
splom( ~ rent[,c(1,2,6,3,4)] | rent$lime, pch=16,
      varname.cex=.8, col=likertColor(2)[2],
      axis.text.cex=.5, pscales=2,
      par.strip.text=list(cex=1.4),
      xlab=NULL)
## hhdev.off()


###################################################
### code chunk number 4: regc.tex:133-140
###################################################
## hhcapture("rent-lm3l.Rout", '
rent.lm3l <-
    lm(rnt.alf ~ rnt.till + cow.dens + prop.past + lime,
       data=rent)
summary(rent.lm3l)
anova(rent.lm3l)
## ')


###################################################
### code chunk number 5: regc.tex:171-176
###################################################
## hhpdf("rent2.pdf", width=8, height=4.5)
splom(~rent[,c(1,2,3)] | rent$lime, pch=16, col=likertColor(2)[2],
      varname.cex=1.2, axis.text.cex=.8, pscales=3,
      par.strip.text=list(cex=1.4), xlab=NULL)
## hhdev.off()


###################################################
### code chunk number 6: regc.tex:192-198
###################################################
## hhcapture("rent-lm4ln.Rout", '
rent.lm4ln <- lm(rnt.alf ~ rnt.till + cow.dens +
                 lime + cow.dens:lime, data=rent)
summary(rent.lm4ln)
anova(rent.lm4ln)
## ')


###################################################
### code chunk number 7: regc.tex:212-223
###################################################
## hhpdf("rent4lnres.pdf", width=7, height=5)
tmp <-
xyplot(resid(rent.lm4ln) ~ rnt.till + cow.dens | lime, groups=lime, data=rent,
       scales=list(alternating=FALSE), pch=19, xlab=NULL, col=likertColor(2))
tmp2 <- update(transpose(tmp),
               scales=list(x=list(relation="free",
                          limits=list(c(7,85), c(0,58)))),
               between=list(x=1),
               xlab=c("rnt.till","cow.dens"))
update(useOuterStrips(combineLimits(tmp2)), strip=FALSE)
## hhdev.off()


###################################################
### code chunk number 8: regc.tex:252-257
###################################################
## hhcapture("rent-lm12p.Rout", '
rent.lm12p <- lm(alf.till ~ lime * cow.dens + prop.past, data=rent)
summary(rent.lm12p)
anova(rent.lm12p)
## ')


###################################################
### code chunk number 9: regc.tex:281-286
###################################################
## hhcapture("rent-lm12.Rout", '
rent.lm12m <- aov(alf.till ~ lime * cow.dens, data=rent)
anova(rent.lm12m)
summary.lm(rent.lm12m)
## ')


###################################################
### code chunk number 10: regc.tex:301-305
###################################################
## hhpdf("rent-lm12m.pdf", width=7, height=3.75)
ancovaplot(alf.till ~ lime * cow.dens, data=rent, col=likertColor(2),
           scales=list(alternating=FALSE), between=list(x=c(0,1)))
## hhdev.off()


###################################################
### code chunk number 11: regc.tex:334-337
###################################################
## hhpdf("rent-plot-lm12m.pdf", height=7, width=9)
lmplot(rent.lm12m, col=likertColor(2))
## hhdev.off()


###################################################
### code chunk number 12: regc.tex:352-359
###################################################
## hhcapture("rent-diag-lm12m.Rout", '
rent.case12m <- case(rent.lm12m)
rent.case12m.trellis <-
   plot(rent.case12m, rent.lm12m, par.strip.text=list(cex=1.2),
        layout=c(3,3), main.cex=1.6, col=likertColor(2)[2], lwd=4)
rent.case12m.trellis ## display both graph and list of noteworthy cases
## ')


###################################################
### code chunk number 13: regc.tex:376-379
###################################################
## hhpdf("rent-diag-lm12m.pdf", width=13, height=9)
rent.case12m.trellis
## hhdev.off()


###################################################
### code chunk number 14: regc.tex:407-422
###################################################
## hhpdf("rent-text-lm12m.pdf", width=8, height=4)
xyplot(alf.till ~ cow.dens | lime, data=rent, RowNames=row.names(rent),
       col=likertColor(2)[2],
       panel=function(x, y, subscripts, RowNames, ...) {
         panel.xyplot(x, y, ...)
         subs <- match(c(19,33,60, 49), subscripts, 0)
         panel.text(x[subs], y[subs],
                    RowNames[subscripts][subs], adj=0, cex=1.5)
       },
       par.strip.text=list(cex=1.4),
       between=list(x=1),
       scales=list(alternating=FALSE, cex=1.2),
       xlab=list(cex=1.4), ylab=list(cex=1.4),
       pch=16, cex=.8)
## hhdev.off()


###################################################
### code chunk number 15: regc.tex:444-449
###################################################
## hhcapture("rent-lm12ms.Rout", '
rent.lm12ms.aov <- aov(alf.till ~ lime * cow.dens,
                       data=rent[-c(19, 33, 60, 49),])
anova(rent.lm12ms.aov)
## ')


###################################################
### code chunk number 16: regc.tex:465-471
###################################################
## hhpdf("rent-lm12ms.pdf", width=7, height=3.75)
ancovaplot(alf.till ~ lime * cow.dens, col=likertColor(2),
           data=rent[-c(19, 33, 60, 49),],
           ylim=range(rent$alf.till, .35, 2.1),
           scales=list(alternating=FALSE), between=list(x=c(0,1)))
## hhdev.off()


###################################################
### code chunk number 17: regc.tex:568-586
###################################################
## hhpdf("rent-residn.pdf", width=7, height=4.5)
## normal plot with straight line, identified outliers, ylim control
qqmath( ~ sort(resid(rent.lm12m)), location.points=c(1,65:67), pch=19,
       col=likertColor(2)[2],
       panel=function(x, location.points=location.points, ...) {
         id.points <- as.numeric(names(x))[location.points]
         panel.qqmathline(x)
         panel.qqmath(x, ...)
         panel.text(x=qnorm(ppoints(x))[location.points],
                    y=x[location.points],
                    names(x)[location.points],
                    adj=0, cex=1.4)
       }
       )
## this graph looks not normal
## also do the test
shapiro.test(resid(rent.lm12m))
## hhdev.off()


###################################################
### code chunk number 18: regc.tex:606-643
###################################################
y.normal <- data.frame(y=rnorm(67*6), group=rep(paste("Normal:", 1:6), each=67))

## hhpdf("norm-prob-plot.pdf", width=5.5, height=4)
qqmath(~ y | group, data=y.normal,
       panel = function(x, ...) {
         panel.qqmathline(x, ...)
         panel.qqmath(x, ...)
       },
       ylab=NULL, pch=19, col=likertColor(2)[2],
       scales=list(alternating=FALSE, tck=c(1,0)),
       main="six randomly generated normal plots",
       between=list(x=1.6, y=1.6),
       par.settings=list(
         layout.heights=list(main.key.padding=0),
         layout.widths=list(axis.key.padding=0, right.padding=0)))
## hhdev.off()

y.nonnormal <- reshape2::melt(data.frame("runif(67)"=runif(67),
                                         "rf(67, 4, 24)"=rf(67, 4, 24),
                                         "rbinom(67, 10, .2)"=rbinom(67, 10, .2),
                                         "rt(67, 2)"=rt(67, 2),
                                         "rchisq(67, 5)"=rchisq(67, 5),
                                         "rpois(67, 8)"=rpois(67, 8),
                                         check.names=FALSE),
                              id.vars=NULL,
                              value.name="y")

## hhpdf("nonnorm-prob-plot.pdf", width=5.5, height=4)
qqmath(~ y | variable, data=y.nonnormal,
       panel = function(x, ...) {
         panel.qqmathline(x, ...)
         panel.qqmath(x, ...)
       },
       ylab=NULL, pch=19, col=likertColor(2)[2],
       scales=list(relation="free"),
       main="six randomly generated nonnormal plots")
## hhdev.off()


###################################################
### code chunk number 19: regc.tex:698-704
###################################################
## hhpdf("rent-resid-plots.pdf", width=6, height=8)
print(A4.left=.019, panel.width=list(4,"cm"),
residual.plots.lattice(rent.lm12m, X=model.matrix(rent.lm12m)[,-1],
                       pch=19, col=likertColor(2)[2])
)
## hhdev.off()


###################################################
### code chunk number 20: regc.tex:867-870
###################################################
## hhpdf("rent-diag-lm12m%03d.pdf", width=7, height=4, onefile=FALSE)
update(rent.case12m.trellis, layout=c(1,1), main=NULL, scales=list(x=list(tck=c(0,1)))) ## all 9, one per file
## hhdev.off()


###################################################
### code chunk number 21: regc.tex:1213-1216
###################################################
## hhcapture("dfbetas.Rout", '
stats:::dfbetas.lm
## ')


###################################################
### code chunk number 22: regc.tex:1248-1251
###################################################
## hhpdf("rent-diag-lm12mNarrow%03d.pdf", width=4.5, height=4, onefile=FALSE)
update(rent.case12m.trellis, layout=c(1,1), main=NULL, scales=list(x=list(tck=c(0,1)))) ## all 9, one per file
## hhdev.off()


###################################################
### code chunk number 23: regc.tex:1290-1294
###################################################
## hhpdf("plot-lm-rent.pdf", width=7, height=7, lwd=2) ## lwd is not an argument to grDevices::pdf
par(mfrow=c(2,2))
plot(rent.lm12m, pch=16, col=likertColor(2)[2])
## hhdev.off()


###################################################
### code chunk number 24: regc.tex:1309-1312
###################################################
## hhpdf("plot-lm-rent-Cook.pdf", width=5, height=5, lwd=2) ## lwd is not an argument to grDevices::pdf
plot(rent.lm12m, 5, add.smooth=FALSE, col=likertColor(2)[2], pch=16)
## hhdev.off()



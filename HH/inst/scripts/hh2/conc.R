### R code from vignette source '~/WindowsC/HOME/rmh/hh.e2/hh2/conc.tex'

###################################################
### code chunk number 1: conc.tex:7-8
###################################################
library(HH)


###################################################
### code chunk number 2: conc.tex:11-16
###################################################
## the standard lattice color 2 is difficult for people with color deficient vision
data(col3x2)
## These colors look like a 3x2 color array when run through
## the vischeck simulator to see how they look for the three most
## common color vision deficiencies: Protanope, Deuteranope, Tritanope.


###################################################
### code chunk number 3: conc.tex:151-176
###################################################
BB <- matrix(c(2,4,4,5), 2, 2,
             dimnames=list(First=c("R.","W."), Second=c(".R",".W")))
require(vcd)
## hhpdf("RW.pdf", height=4.5, width=4.5)
mosaic(BB,
       rot_labels=c(0,0,0,0),  ## zero is horizontal
       rot_varnames=c(0,0,0,0),
       gp=gpar(fill=c("red","#ff8888","#ff8888","white"), col="black"),
       spacing=spacing_highlighting(.5)
       )
grid.text("A",                   x=.92, y=.343)
grid.text(expression(bar(A)),    x=.92, y=.727)
grid.text("B",                   x=.62, y=.070)
grid.text(expression(bar(B)),    x=.26, y=.070)
grid.text(expression(scriptstyle(RR~""==""~bar(A)*intersect(bar(B)))), x=.255, y=.752)
grid.text(expression(scriptstyle(RW~""==""~bar(A)*intersect(B))),      x=.624, y=.752)
grid.text(expression(scriptstyle(WR~""==""~A*intersect(bar(B)))),      x=.295, y=.368)
grid.text(expression(scriptstyle(WW~""==""~A*intersect(B))),           x=.662, y=.368)
grid.text(expression(scriptstyle((4/10)(3/9)~""==""~2/15)),            x=.255, y=.702)
grid.text(expression(scriptstyle((4/10)(6/9)~""==""~4/15)),            x=.624, y=.702)
grid.text(expression(scriptstyle((6/10)(4/9)~""==""~4/15)),            x=.295, y=.318)
grid.text(expression(scriptstyle((6/10)(5/9)~""==""~5/15)),            x=.662, y=.318)
## hhdev.off()
## The x and y values in the grid.text statements probably will not align correctly on
## any graphics device or size other than the one listed in the hhpdf statement.


###################################################
### code chunk number 4: conc.tex:275-280
###################################################
## hhpdf("discuniv.pdf", width=3, height=2)
discuniv <- data.frame(y=c(.25, .50, .25), x=factor(0:2))
barchart(y ~ x, data=discuniv, ylab="P(X = x)", xlab="x", origin=0, col="#ff4444",
         scales=list(y=list(at=c(0, .25, .50))))
## hhdev.off()


###################################################
### code chunk number 5: conc.tex:316-343
###################################################
require(reshape2)
pmf <- matrix(c(.10, .05, .20, .10, .30, .25), 2, 3,
              dimnames=list(x=1:2, y=0:2))
## hhpdf("discbiv1.pdf", width=2.5, height=2)
barchart(value ~ rep("", 6) | y*x,
         data=melt(pmf, as.is=TRUE),
         origin=0, box.ratio=100,
         scales=list(x=list(limits=c(.65, 1.35), alternating=0),
                     y=list(relation="free", limits=c(0,.35),
                            at=.15, labels=list("1", "","","2","",""))),
         strip=FALSE,
         par.settings=list(axis.line=list(col=0), clip=list(panel=FALSE)),
         as.table=TRUE, col="#ff4444",
         xlab.top="y", ylab=list("x", rot=0)) +
    layer(panel.axis(side="top", at=1, rot=0, tck=0, outside=TRUE,
          labels=c("0","1","2","","","")[panel.number()]))
## hhdev.off()
## hhpdf("discbiv2.pdf", width=3, height=2)
mosaic(t(pmf), split_vertical=c(TRUE,FALSE),
       highlighting=2, highlighting_fill=c("#ff8888","red"),
       rot_labels=c(0,0,0,0), rot_varnames=c(0,0,0,0), keep_aspect_ratio=FALSE)
## hhdev.off()
## hhpdf("discbiv3.pdf", width=3, height=2)
mosaic(pmf,
       highlighting=2, highlighting_fill=c("#ffBBBB","#ff8888","red"),
       rot_labels=c(0,0,0,0), rot_varnames=c(0,0,0,0), keep_aspect_ratio=FALSE)
## hhdev.off()


###################################################
### code chunk number 6: conc.tex:453-471
###################################################
## hhpdf("bimodal-shade.pdf", height=4, width=7, col=likertColor(2)[2], lwd=3) ## col and lwd are not arguments for grDevices:::pdf
xx <- seq(-3, 6, .025)
dd <- (dnorm(xx, mean=0, sd=1) + dnorm(xx, mean=2.5, sd=1.1))/2
pr <- (pnorm(c(2,4), mean=0, sd=1) + pnorm(c(2,4), mean=2.5, sd=1.1))/2
xyplot(dd ~ xx, type="l",
       panel=function(x, y, ...) {
         grid.polygon(x=x[c(201,201:281,281)],
                      y=c(0,y[201:281],0),
                      gp=gpar(fill="gray60", col="gray60"),
                      default.units="native")
          panel.xyplot(x, y, ...)
         panel.abline(h=0, lty=1, col="gray70")
      },
       ylab=list("f(x)", rot=0),
       xlab=list("x"),
       main=list(paste("Prob(2 < X < 4) =", round(diff(pr),3)))
       )
## hhdev.off()


###################################################
### code chunk number 7: conc.tex:634-651
###################################################
pp <- ppoints(101)

zz.norm <- qnorm(pp, s=2.2)
dd.norm <- dnorm(zz.norm, s=2.2)
N <- xyplot(dd.norm ~ zz.norm, type="l") + layer(panel.abline(h=0, col="gray70"))

xx.chisq4 <- c(0,qchisq(pp, 4))
dd.chisq4 <- dchisq(xx.chisq4, 4)
C4 <- xyplot(dd.chisq4 ~ xx.chisq4, type="l") + layer(panel.abline(h=0, col="gray70"))
C4R <- xyplot(rev(dd.chisq4) ~ -rev(xx.chisq4), type="l") + layer(panel.abline(h=0, col="gray70"))

## hhpdf("skewdens2.pdf", height=3.5, width=8, col=likertColor(2)[2], lwd=3) ## col and lwd are not arguments for grDevices:::pdf
update(c("negatively skewed"=C4R,
         symmetric=N, "positively skewed"=C4,
         layout=c(3,1), y.same=TRUE),
       between=list(x=1), xlab="x", ylab="density")
## hhdev.off()


###################################################
### code chunk number 8: conc.tex:704-713
###################################################
## hhcapture("tv-freq.Rout", '
data(tv)
tmp <- as.matrix(table(cut(tv$male.life.exp, breaks=seq(49.5,79.5,5))))
dimnames(tmp) <-
    list("Male Life Expectancy"=
              c("50--54","55--59","60--64","65--69","70--74","75--79"),
         " "="Frequency")
tmp
## ')


###################################################
### code chunk number 9: conc.tex:760-769
###################################################
## hhpdf("tv-hist.pdf", height=3.5, width=7)
tmp <-
histogram( ~ male.life.exp, data = tv,
          breaks=seq(49.5, 79.5, 5), type="count", col="gray60")
update(tmp, par.settings=list(clip=list(panel=FALSE),
                              layout.widths=list(axis.right=1.5, right.padding=1.5)),
       ylab.right="Proportion") +
           layer(panel.axis("right", at=seq(0,10,2), labels=seq(0,10,2)/40, outside=TRUE))
## hhdev.off()


###################################################
### code chunk number 10: conc.tex:806-809
###################################################
## hhcapture("conc-stem-male.Rout", '
stem(tv$male.life.exp)
## ')


###################################################
### code chunk number 11: conc.tex:906-909
###################################################
## hhcapture("conc-quartiles-male.Rout", '
quantile(tv$male.life.exp)
## ')


###################################################
### code chunk number 12: conc.tex:920-923
###################################################
## hhpdf("tv-bw.pdf", height=2.5, width=7, col=likertColor(2)[2]) ## col is not an argument for grDevices:::pdf
bwplot( ~ male.life.exp, data = tv)
## hhdev.off()


###################################################
### code chunk number 13: conc.tex:940-974
###################################################
pp <- ppoints(100)
pp <- c(0, pp[1:25], .25, pp[26:50], .50, pp[51:75], .75, pp[76:100])
q.pp <- qf(pp, df1=3, df2=36)
dd <- df(q.pp, df1=3, df2=36)

## hhpdf("quartiles.pdf", height=4.5, width=7, col=likertColor(2)[2]) ## col is not an argument for grDevices:::pdf
     xyplot(dd ~ q.pp, type="l",
            par.settings = list(clip = list(panel = "off")),
ylim=c(-.1, .78),
            panel=function(x, y, ...) {
              grid.polygon(x=x[c(1,1:27,27)],
                           y=c(0,y[1:27],0),
                           gp=gpar(fill="gray60", col="gray60"), default.units="native")
              grid.polygon(x=x[c(53,53:79,79)],
                           y=c(0,y[53:79],0),
                           gp=gpar(fill="gray60", col="gray60"), default.units="native")
              panel.axis("bottom", at=x[c(27,53,79)],
                         tck = .5, # line.col = "transparent",
                         labels=c("Q1","Med","Q3"),
                         half=FALSE,
                         rot=0, text.cex=1)
              panel.axis("bottom", at=x[c(27,53,79)],
                         tck = 2.5,
                         labels=round(x[c(27,53,79)], 2),
                         outside=TRUE,
                         rot=0, text.cex=1)
              panel.abline(h=0, col="gray70")
              panel.xyplot(x, y, ...)
            },
            ylab=list("density", cex=1.6),
            xlab=list("quantile", cex=1.6),
            main=list("Quartiles of F(3,36)", cex=1.5)
            )
## hhdev.off()


###################################################
### code chunk number 14: conc.tex:992-1011
###################################################
## hhpdf("skew.pdf", height=3, width=7, col=likertColor(2)[2]) ## col is not an argument for grDevices:::pdf
sym <- rnorm(100)

neg <- sym[sym<0]
pos <- sym[sym>0]

neg.skew <- c(-neg^2, pos^.5)

pos.skew <- c(-(-neg)^.5, pos^2)

skew.levels <- c("negatively skewed", "symmetric", "positively skewed")
skew.df <- data.frame(y=c(neg.skew, sym, pos.skew),
                      dist=ordered(rep(skew.levels, c(100,100,100)),
                        levels=skew.levels))

bwplot(dist ~ y, data=skew.df, xlab="",
 par.settings=list(plot.symbol=list(pch=19),
                   box.dot=list(pch=19, col=trellis.par.get()$plot.symbol$col)))
## hhdev.off()


###################################################
### code chunk number 15: conc.tex:1083-1099
###################################################
## hhpdf("corr-eps.pdf", height=2.5, width=8)
x <- rnorm(100)
e <- rnorm(100)
r <- c(-1, -.9, -.5, 0, .5, .9, 1)

corr.data <- data.frame(e, x, x %o% r + e %o% (1-r^2)^.5)
names(corr.data)[3:9] <- r

corr.data.melt <- reshape2::melt(corr.data[,-1], id="x", variable.name="correlation", value.name="y")
maxabs <- c(-1,1) * 3.9
xyplot(y ~ x | correlation, data=corr.data.melt, layout=c(7, 1), aspect="iso",
       col=likertColor(2)[2],
       xlim=maxabs, ylim=maxabs, ylab=list(rot=0),
       scales=list(at=c(-2,0,2), alternating=FALSE),
       strip=strip.custom(var.name=expression(rho), strip.names=c(TRUE,TRUE), sep=" = "))
## hhdev.off()


###################################################
### code chunk number 16: conc.tex:1144-1151
###################################################
## hhpdf("bivnorm8.pdf", height=8, width=8)
bv8 <- bivariateNormal(.7)  ## all views on one page
bv8
## hhdev.off()
## hhpdf("bivnorm1125.pdf", height=4, width=4)
update(bv8[3], layout=c(1,1))
## hhdev.off()


###################################################
### code chunk number 17: conc.tex:1259-1273
###################################################
## hhpdf("normApproxBin.pdf", height=4.5, width=6.5)
BB <- HH:::dstrplotDiscrete("binom", ok0=TRUE, args=list(size=15, prob=.4),
      x.tick.number=15) +
      layer(HH:::panel.dstrplotDiscreteFill(..., X=6))

NT2 <- NTplot(xbar=.4*15, n=1, sd=sqrt(15*.4*(1-.4)), type="confidence",
                     alpha.left=.025, alpha.right=0.025, float=FALSE, cex.top.axis=1)

update(BB, ylim=c(-.03, .24), box.width=.99) +
as.layer(NT2, under=TRUE) +
 layer({panel.abline(v=c(5.5, 6.5), lty=5, col="gray30", lwd=.5)
        panel.axis("bottom", at=c(5.5, 6.5), tck=.7, half=FALSE,
                   text.cex=.6, line.col="transparent")}, under=TRUE)
## hhdev.off()


###################################################
### code chunk number 18: conc.tex:1275-1281
###################################################
## hhcapture("normApproxBin.Rout", '
pbinom(size=15, prob=.4, q=6)
pnorm(q=6.5, mean=15*.4, sd=sqrt(15*.4*(1-.4)))
dbinom(size=15, prob=.4, x=6)
diff(pnorm(q=c(5.5, 6.5), mean=15*.4, sd=sqrt(15*.4*(1-.4))))
## ')


###################################################
### code chunk number 19: conc.tex:1354-1383
###################################################
## hhpdf("normalPhi.pdf", height=6.5, width=6.5)
AA <- HH:::dstrplotContinuous("norm", args=list(m=0, s=1), X=1.645, par.settings=list(clip=list(panel=FALSE))) +
        layer({
          panel.axis("left", at=.1031, rot=0, outside=TRUE, tck=5)
          panel.abline(h=.1031, lty=2, col="gray50")
          HH:::panel.dstrplotContinuousFill(..., X=1.645)
        })
qq <- seq(-26, 26, 1)/10
BB <- xyplot(pnorm(qq) ~ qq, type="l", col="black", xlab="z", ylim=c(0,1),
             scales=list(y=list(at=seq(0,1,.2), rot=0)),
             par.settings=list(clip=list(panel=FALSE))) +
             layer({
               panel.axis("left", at=.95, rot=0, outside=TRUE, tck=7)
               panel.axis("bottom", at=1.645, outside=TRUE, tck=3)
               ## panel.segments(1.645, c(.95, 0), 1.645, c(1, .95),
               ##                col=c("skyblue1", "skyblue3"), lwd=8, lend=1)
               ## lend=1 not honored by panel.segments
               panel.rect(1.645-.03, c(.95, 0), 1.645+.03, c(1, .95),
                          col=c("skyblue1", "skyblue3"),
                          border=c("skyblue1", "skyblue3"))
               panel.abline(h=.95, lty=2, col="gray50")
             })
BBAA <- resizePanels(c(BB=BB, AA=AA, layout=c(1,2), x.same=TRUE), h=c(.35, .65))
update(BBAA, ylab=NULL, xlim=c(-2.58, 2.58),
      ylim=list(c(.06, .94), c(.024, .38)),
      strip.left=strip.custom(factor.levels=c(expression("probability CDF"~Phi(z)), expression("density"~phi(z)))),
      par.strip.text=list(lines=1.5),
      strip=FALSE, between=list(y=2))
## hhdev.off()


###################################################
### code chunk number 20: conc.tex:1386-1391
###################################################
## hhcapture("normalPhi.Rout", '
  dnorm(1.645, m=0, s=1)
  pnorm(1.645, m=0, s=1)
  qnorm(0.95, m=0, s=1)
## ')


###################################################
### code chunk number 21: conc.tex:1454-1460
###################################################
## hhpdf("norm.pdf", height=5.5, width=7)
tmp.norm <-
NTplot(mean0=100, mean1=NA,  xbar=NA,  xlim=c(75, 125), sd=5,
       digits=6, zaxis=TRUE, cex.z=0.6, cex.prob=.9)
print(tmp.norm, cex.table=.8)
## hhdev.off()


###################################################
### code chunk number 22: conc.tex:1462-1468
###################################################
## These lines show how to access the xbar and z values and their p-values
## from within the R session.
## hhcapture("norm.Rout", '
attr(tmp.norm, "scales")
attr(tmp.norm, "prob")
## ')


###################################################
### code chunk number 23: conc.tex:1498-1533
###################################################
## hhpdf("tt.pdf", height=9, width=8)
normal  <- NTplot(mean0=100, mean1=NA,  xbar=NA,
                 xlim=c(75, 125), sd=5, digits=6, distribution.name="z",
                 zaxis=TRUE, cex.z=0.6, cex.prob=.9, xhalf.multiplier=1.2,
                 key.axis.padding=6)
`t[30]` <- NTplot(mean0=100, mean1=NA,  xbar=NA,
                  xlim=c(75, 125), sd=5, digits=6, distribution.name="t", df=30,
                  zaxis=TRUE, cex.z=0.6, cex.prob=.9, xhalf.multiplier=1.2)
`t[10]` <- NTplot(mean0=100, mean1=NA,  xbar=NA,
                  xlim=c(75, 125), sd=5, digits=6, distribution.name="t", df=10,
                  zaxis=TRUE, cex.z=0.6, cex.prob=.9, xhalf.multiplier=1.2)
`t[2]`  <- NTplot(mean0=100, mean1=NA,  xbar=NA,
                  xlim=c(75, 125), sd=5, digits=6, distribution.name="t", df=2,
                  zaxis=TRUE, cex.z=0.6, cex.prob=.9, xhalf.multiplier=1.2)

scales <-
cbind(
normal  = attr(normal , "scales")[,2],
`t[30]` = attr(`t[30]`, "scales")[,2],
`t[10]` = attr(`t[10]`, "scales")[,2],
`t[2]`  = attr(`t[2]` , "scales")[,2]
)
rownames(scales) <- c(colnames(attr(normal , "scales"))[2], "t")

NTTT <-
  update(
    c(normal=normal, "t[30]"=`t[30]`, "t[10]"=`t[10]`, "t[2]"=`t[2]`,
      layout=c(2,2), y.same=TRUE, x.same=FALSE),
    strip=FALSE, strip.left=TRUE, scales=list(y=list(alternating=1)),
    between=list(x=2, y=4), ylab="density",
    main=expression("normal and three t distributions, " ~ sigma[bar(x)]==5 ~ "," ~ n==1)) +
  layer(panel.abline(h=.08, col="gray60", lty=2, lwd=.5))
attr(NTTT, "scales") <- scales
print(NTTT, digits=6, cex.table=.9)
## hhdev.off()


###################################################
### code chunk number 24: conc.tex:1624-1653
###################################################
## hhpdf("normalSamplingDist.pdf", height=8, width=8)
N1  <- reshape2::melt(data.frame(matrix(rnorm( 10, mean=100, sd=5), nrow= 1, ncol=10)), id=NULL)
N4  <- reshape2::melt(data.frame(matrix(rnorm( 40, mean=100, sd=5), nrow= 4, ncol=10)), id=NULL)
N16 <- reshape2::melt(data.frame(matrix(rnorm(160, mean=100, sd=5), nrow=16, ncol=10)), id=NULL)
N64 <- reshape2::melt(data.frame(matrix(rnorm(640, mean=100, sd=5), nrow=64, ncol=10)), id=NULL)


blue4Transparent <- HH:::ColorWithAlpha("blue4", 150)
update(strip=FALSE, strip.left=TRUE, between=list(x=2, y=3), xlab=NULL,
       scale=list(x=list(limits=c(84,116), alternating=1), y=list(at=1:10, labels=1:10)),
       key=list(
         text=list(expression(x, bar(x))),
         points=list(pch=c(1, 17), col=c("DodgerBlue1", blue4Transparent), cex=c(1.2, 1.3)),
         columns=2, space="bottom"),
c("n=1"=dotplot(variable ~ value, data=N1,  pch=1, cex=1.2, col="DodgerBlue1") +
    xyplot(variable ~ value, data=aggregate(value ~ variable, data=N1,  FUN=mean),
           cex=1.3, pch=17, col=blue4Transparent),
  "n=4"=dotplot(variable ~ value, data=N4,  pch=1, cex=1.2, col="DodgerBlue1") +
    xyplot(variable ~ value, data=aggregate(value ~ variable, data=N4,  FUN=mean),
           cex=1.3, pch=17, col=blue4Transparent),
  "n=16"=dotplot(variable ~ value, data=N16, pch=1, cex=1.2, col="DodgerBlue1") +
    xyplot(variable ~ value, data=aggregate(value ~ variable, data=N16, FUN=mean),
           cex=1.3, pch=17, col=blue4Transparent),
  "n=64"=dotplot(variable ~ value, data=N64, pch=1, cex=1.2, col="DodgerBlue1") +
    xyplot(variable ~ value, data=aggregate(value ~ variable, data=N64, FUN=mean),
           cex=1.3, pch=17, col=blue4Transparent),
  layout=c(2, 2), x.same=TRUE, y.same=TRUE) +
    layer(panel.abline(v=c(92,98,102,108), col="gray60")))
## hhdev.off()


###################################################
### code chunk number 25: conc.tex:1681-1721
###################################################
## hhpdf("normalCLT.pdf", height=9, width=8)
n1  <- NTplot(mean0=100, mean1=NA,  xbar=NA,
              xlim=c(85, 115), sd=5, n=1,  digits=6,
              zaxis=TRUE, cex.z=0.6, cex.prob=.9, ylim=c(0, .65),
              xhalf.multiplier=1.2, key.axis.padding=6)
n4  <- NTplot(mean0=100, mean1=NA,  xbar=NA,
              xlim=c(85, 115), sd=5, n=4,  digits=6,
              zaxis=TRUE, cex.z=0.6, cex.prob=.9, ylim=c(0, .65),
              xhalf.multiplier=1.2)
n16 <- NTplot(mean0=100, mean1=NA,  xbar=NA,
              xlim=c(85, 115), sd=5, n=16, digits=6,
              zaxis=list(at=seq(85, 115, 5), labels=seq(-12, 12, 4)),
              cex.z=0.6, cex.prob=.9, ylim=c(0, .65),
              xhalf.multiplier=1.2)
n64 <- NTplot(mean0=100, mean1=NA,  xbar=NA,
              xlim=c(85, 115), sd=5, n=64, digits=6,
              zaxis=list(at=seq(85, 115, 5), labels=seq(-24, 24, 8)),
              cex.z=0.6, cex.prob=.9, ylim=c(0, .65),
              xhalf.multiplier=1.2)

scales <-
cbind(
"n==1"  = c(attr(n1 , "scales")[,2], attr(n1 , "table")["bar(x)", "sigma"]),
"n==4"  = c(attr(n4 , "scales")[,2], attr(n4 , "table")["bar(x)", "sigma"]),
"n==16" = c(attr(n16, "scales")[,2], attr(n16, "table")["bar(x)", "sigma"]),
"n==64" = c(attr(n64, "scales")[,2], attr(n64, "table")["bar(x)", "sigma"])
)
rownames(scales)[1] <- c(colnames(attr(n1 , "scales"))[2])
rownames(scales)[3] <- "sigma[bar(x)]"

NNNN <-
  update(
    c("n=1"=n1, "n=4"=n4, "n=16"=n16, "n=64"=n64,
      layout=c(2,2), y.same=TRUE, x.same=FALSE),
    strip=FALSE, strip.left=TRUE, scales=list(y=list(alternating=1)),
    between=list(x=2, y=4), ylab="density",
    main=expression("Normal with increasing sample size n, " ~ sigma[bar(x)]==5/sqrt(n)))
attr(NNNN, "scales") <- scales
print(NNNN, digits=6, cex.table=.9)
## hhdev.off()


###################################################
### code chunk number 26: conc.tex:1861-1871
###################################################
## hhpdf("normalconf.pdf", height=6, width=8)
tmp <-
NTplot(xbar=8.5, sd=2, df=24, n=25,
       xlim=c(7,10), ylim=c(0,0.96),
       alpha.right=0.025, alpha.left=0.025,
       digits=5,
       distribution.name="t", type="confidence",
       z1axis=TRUE, zaxis=TRUE, cex.z=0.7, cex.prob=1)
print(tmp, cex.table=.8)
## hhdev.off()


###################################################
### code chunk number 27: conc.tex:2051-2061
###################################################
## hhpdf("NormConf.pdf", height=6, width=8)
tmp <-
NTplot(xbar=8.5, sd=2, n=25,
       xlim=c(7,10), ylim=c(0,0.96),
       alpha.right=0.025, alpha.left=0.025,
       digits=5,
       distribution.name="z", type="confidence",
       z1axis=TRUE, zaxis=TRUE, cex.z=0.7, cex.prob=1)
print(tmp, cex.table=.8)
## hhdev.off()


###################################################
### code chunk number 28: conc.tex:2095-2105
###################################################
## hhpdf("NormConfRight.pdf", height=6, width=8)
tmp <-
NTplot(xbar=8.5, sd=2, n=25,
       xlim=c(7,10), ylim=c(0,0.96),
       alpha.right=0, alpha.left=0.05,
       digits=5,
       distribution.name="z", type="confidence",
       z1axis=TRUE, zaxis=TRUE, cex.z=0.7, cex.prob=1)
print(tmp, cex.table=.8)
## hhdev.off()


###################################################
### code chunk number 29: conc.tex:2243-2281
###################################################
tmp32  <- NTplot(mean0=8, mean1=8.411, sd=2, n=32,
                 xlim=c(7.5, 8.9), ylim=c(0, 2.2), cex.top.axis=1.8, cex.prob=1.5,
                 cex.z=0.7, prob.labels=FALSE,
                 digits.axis=5, digits.float=3, digits.left=3,
                 power=TRUE, beta=TRUE, col.power="#FF808B",
                 xhalf.multiplier=.8, yhalf.multiplier=.6)
## tmp32

tmp64  <- NTplot(mean0=8, mean1=8.411, sd=2, n=64,
                 xlim=c(7.5, 8.9), ylim=c(0, 2.2), cex.top.axis=1.8, cex.prob=1.5,
                 cex.z=0.7, prob.labels=FALSE,
                 digits.axis=5, digits.float=3, digits.left=3,
                 power=TRUE, beta=TRUE, col.power="#FF808B",
                 xhalf.multiplier=.8, yhalf.multiplier=.6)
## tmp64

tmp128 <- NTplot(mean0=8, mean1=8.411, sd=2, n=128,
                 xlim=c(7.5, 8.9), ylim=c(0, 2.2), cex.top.axis=1.8, cex.prob=1.5,
                 cex.z=0.7, prob.labels=FALSE,
                 digits.axis=5, digits.float=3, digits.left=3,
                 power=TRUE, beta=TRUE, col.power="#FF808B", ## darker than pink
                 xhalf.multiplier=.8, yhalf.multiplier=.6)
## tmp128

## hhpdf("powern.pdf", height=11, width=19) ## density, power, beta by sample size

print(position=c(.010, 0, .334, .95), more=TRUE,
      update(tmp32, ylab=NULL, main=list('n = 32\n', cex=2),
             strip.left=TRUE, par.strip.text=list(cex=1.5)),
      tables=FALSE)
print(position=c(.357, 0, .667, .95), more=TRUE,
      update(tmp64, ylab=NULL, main=list('n = 64\n', cex=2), strip.left=FALSE),
      tables=FALSE)
print(position=c(.690, 0, 1.000, .95), more=FALSE,
      update(tmp128, ylab=NULL, main=list('n = 128\n', cex=2), strip.left=FALSE),
      tables=FALSE)

## hhdev.off()


###################################################
### code chunk number 30: conc.tex:2497-2505
###################################################
## hhpdf("bottlefill.pdf", height=5, width=8)
tmp <-
NTplot(mean0=32, xbar=31.94, sd=.3, n=100,
       xlim=c(31.88,32.12), ylim=c(0,13),
       alpha.right=0.005, alpha.left=0.005,
       zaxis=TRUE, cex.z=0.7)
print(tmp, cex.table=.8)
## hhdev.off()


###################################################
### code chunk number 31: conc.tex:2582-2590
###################################################
## hhpdf("bottlefillonetail.pdf", height=5, width=8)
tmp <-
NTplot(mean0=32, xbar=31.94, sd=.3, n=100,
       xlim=c(31.88,32.12), ylim=c(0,13),
       alpha.right=0, alpha.left=0.01,
       zaxis=TRUE, cex.z=0.7)
print(tmp, cex.table=.8)
## hhdev.off()


###################################################
### code chunk number 32: conc.tex:2704-2741
###################################################
tmp8     <- NTplot(mean0=8, mean1=8, sd=2, n=64,
                   xlim=c(7.3, 9.5), cex.top.axis=1.8, cex.prob=1.5,
                   cex.z=0.7, prob.labels=FALSE,
                   digits.axis=5, digits.float=3, digits.left=3,
                   power=TRUE, beta=TRUE, col.power="#FF808B",
                   xhalf.multiplier=.8, yhalf.multiplier=.6)
## tmp8

tmp8411  <- NTplot(mean0=8, mean1=8.411, sd=2, n=64,
                   xlim=c(7.3, 9.5), cex.top.axis=1.8, cex.prob=1.5,
                   cex.z=0.7, prob.labels=FALSE,
                   digits.axis=5, digits.float=3, digits.left=3,
                   power=TRUE, beta=TRUE, col.power="#FF808B",
                   xhalf.multiplier=.8, yhalf.multiplier=.6)
## tmp8411

tmp87314 <- NTplot(mean0=8, mean1=8.7314, sd=2, n=64, xbar=NA,
                   xlim=c(7.3, 9.5), cex.top.axis=1.8, cex.prob=1.5,
                   cex.z=0.7, prob.labels=FALSE,
                   digits.axis=5, digits.float=3, digits.left=3,
                   power=TRUE, beta=TRUE, col.power="#FF808B", ## darker than pink
                   xhalf.multiplier=.8, yhalf.multiplier=.6)
## tmp87314

## hhpdf("powerbeta64.pdf", height=11, width=19) ## density, power, beta: by mean1

print(position=c(.010, 0, .334, .95), more=TRUE,
      update(tmp8, ylab=NULL, main=NULL, strip.left=TRUE, par.strip.text=list(cex=1.5)),
      tables=FALSE)
print(position=c(.357, 0, .667, .95), more=TRUE,
      update(tmp8411, ylab=NULL, main=NULL, strip.left=FALSE),
      tables=FALSE)
print(position=c(.690, 0, 1.000, .95), more=FALSE,
      update(tmp87314, ylab=NULL, main=NULL, strip.left=FALSE),
      tables=FALSE)

## hhdev.off()


###################################################
### code chunk number 33: conc.tex:2893-2913
###################################################
tmpn <- NTplot(mean0=8, mean1=9.4, sd=2, n=12, cex.main=1.4,
                   xlim=c(7, 10.5), ylim=c(0, .7), cex.top.axis=1.5, cex.prob=1.5,
                   cex.z=0.7, prob.labels=FALSE,
                   digits.axis=5, digits.float=3, digits.left=3,
                   power=TRUE, beta=TRUE, col.power="#FF808B", ## darker than pink
                   xhalf.multiplier=.8, yhalf.multiplier=.6, hh=c(7, 1.5, 1.5))
tmpt <- NTplot(mean0=8, mean1=9.4, sd=2, n=12, cex.main=1.4, distribution.name="t", df=11,
                   xlim=c(7, 10.5), ylim=c(0, .7), cex.top.axis=1.5, cex.prob=1.5,
                   cex.z=0.7, prob.labels=FALSE,
                   digits.axis=5, digits.float=3, digits.left=3,
                   power=TRUE, beta=TRUE, col.power="#FF808B", ## darker than pink
                   xhalf.multiplier=.8, yhalf.multiplier=.6, hh=c(7, 1.5, 1.5))
## hhpdf("ncpt.pdf", height=10, width=12) ## density, power, beta: by mean1
print(position=c(.010, 0, .520, .95), more=TRUE,
      update(tmpn, strip.left=TRUE, par.strip.text=list(cex=1.5)),
      tables=FALSE)
print(position=c(.515, 0, 1.000, .95), more=FALSE,
      update(tmpt, strip.left=FALSE),
      tables=FALSE)
## hhdev.off()


###################################################
### code chunk number 34: conc.tex:2933-2939
###################################################
## hhcode("power.t.test.R", '
PowerT <- power.t.test(n=12, sd=2, delta=1.4,
                       type="one.sample",
                       alternative="one.sided")
NTplot(PowerT, beta=TRUE, power=TRUE)
## ')



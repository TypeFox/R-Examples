### R code from vignette source '~/WindowsC/HOME/rmh/hh.e2/hh2/dstr.tex'

###################################################
### code chunk number 1: dstr.tex:8-9
###################################################
library(HH)


###################################################
### code chunk number 2: dstr.tex:82-86
###################################################
##   hhpdf("dstr-beta.pdf", height=3, width=6.5)
  HH:::dstrplotContinuous("beta", args=list(shape1=85.5, shape2=15.5), ok0=TRUE, ok1=TRUE, X=.85) +
    layer(HH:::panel.dstrplotContinuousFill(..., X=.85))
##   hhdev.off()


###################################################
### code chunk number 3: dstr.tex:88-93
###################################################
## hhcapture("dstr-beta.Rout", '
  dbeta(.85, shape1=85.5, shape2=15.5)
  pbeta(.85, shape1=85.5, shape2=15.5)
  qbeta(0.5131489, shape1=85.5, shape2=15.5)
## ')


###################################################
### code chunk number 4: dstr.tex:109-113
###################################################
##   hhpdf("dstr-cauchy.pdf", height=3, width=6.5)
  HH:::dstrplotContinuous("cauchy", xlim=c(-21, 21), X=1.96) +
    layer(HH:::panel.dstrplotContinuousFill(..., X=1.96))
##   hhdev.off()


###################################################
### code chunk number 5: dstr.tex:115-120
###################################################
## hhcapture("dstr-cauchy.Rout", '
  dcauchy(1.96)
  pcauchy(1.96)
  qcauchy(0.8498286)
## ')


###################################################
### code chunk number 6: dstr.tex:131-135
###################################################
##   hhpdf("dstr-chisq.pdf", height=3, width=6.5)
  HH:::dstrplotContinuous("chisq", ok0=TRUE, args=list(df=10), X=18.31) +
    layer(HH:::panel.dstrplotContinuousFill(..., X=18.31))
##   hhdev.off()


###################################################
### code chunk number 7: dstr.tex:137-142
###################################################
## hhcapture("dstr-chisq.Rout", '
  dchisq(18.31, df=10)
  pchisq(18.31, df=10)
  qchisq(0.9500458, df=10)
## ')


###################################################
### code chunk number 8: dstr.tex:173-177
###################################################
##   hhpdf("dstr-exp.pdf", height=3, width=6.5)
  HH:::dstrplotContinuous("exp", ok0=TRUE, args=list(rate=.6), X=1) +
    layer(HH:::panel.dstrplotContinuousFill(..., X=1))
##   hhdev.off()


###################################################
### code chunk number 9: dstr.tex:179-184
###################################################
## hhcapture("dstr-exp.Rout", '
  dexp(1, rate=.6)
  pexp(1, rate=.6)
  qexp(0.4511884, rate=.6)
## ')


###################################################
### code chunk number 10: dstr.tex:212-216
###################################################
##   hhpdf("dstr-f.pdf", height=3, width=6.5)
  HH:::dstrplotContinuous("f", ok0=TRUE, args=list(df1=4, df2=20), X=3) +
    layer(HH:::panel.dstrplotContinuousFill(..., X=3))
##   hhdev.off()


###################################################
### code chunk number 11: dstr.tex:218-223
###################################################
## hhcapture("dstr-f.Rout", '
  df(3, df1=4, df2=20)
  pf(3, df1=4, df2=20)
  qf(0.956799, df1=4, df2=20)
## ')


###################################################
### code chunk number 12: dstr.tex:249-253
###################################################
##   hhpdf("dstr-gamma.pdf", height=3, width=6.5)
  HH:::dstrplotContinuous("gamma", ok0=TRUE, args=list(shape=3), X=6) +
    layer(HH:::panel.dstrplotContinuousFill(..., X=6))
##   hhdev.off()


###################################################
### code chunk number 13: dstr.tex:255-260
###################################################
## hhcapture("dstr-gamma.Rout", '
  dgamma(6, shape=3)
  pgamma(6, shape=3)
  qgamma(0.9380312, shape=3)
## ')


###################################################
### code chunk number 14: dstr.tex:305-309
###################################################
##   hhpdf("dstr-lnorm.pdf", height=3, width=6.5)
  HH:::dstrplotContinuous("lnorm", ok0=TRUE, X=5) +
    layer(HH:::panel.dstrplotContinuousFill(..., X=5))
##   hhdev.off()


###################################################
### code chunk number 15: dstr.tex:311-316
###################################################
## hhcapture("dstr-lnorm.Rout", '
  dlnorm(5)
  plnorm(5)
  qlnorm(0.9462397)
## ')


###################################################
### code chunk number 16: dstr.tex:341-345
###################################################
##   hhpdf("dstr-logis.pdf", height=3, width=6.5)
  HH:::dstrplotContinuous("logis", X=2) +
    layer(HH:::panel.dstrplotContinuousFill(..., X=2))
##   hhdev.off()


###################################################
### code chunk number 17: dstr.tex:347-352
###################################################
## hhcapture("dstr-logis.Rout", '
  dlogis(2)
  plogis(2)
  qlogis(0.8807971)
## ')


###################################################
### code chunk number 18: dstr.tex:381-385
###################################################
##   hhpdf("dstr-norm.pdf", height=3, width=6.5)
  HH:::dstrplotContinuous("norm", args=list(m=0, s=1), X=1.645) +
    layer(HH:::panel.dstrplotContinuousFill(..., X=1.645))
##   hhdev.off()


###################################################
### code chunk number 19: dstr.tex:387-392
###################################################
## hhcapture("dstr-norm.Rout", '
  dnorm(1.645, m=0, s=1)
  pnorm(1.645, m=0, s=1)
  qnorm(0.95, m=0, s=1)
## ')


###################################################
### code chunk number 20: dstr.tex:433-444
###################################################
##   hhpdf("dstr-tukey.pdf", height=3, width=6.5)
  qX <- 4.199
  pX <- ptukey(qX, nmeans = 4, df=12)
  pp <- sort(c(0, ppoints(101), pX))
  qq <- qtukey(pp, nmeans = 4, df=12)
  xyplot(pp ~ qq, type="l", xlab="x", ylab=NULL, main="ptukey(x, nmeans=4, df=12)") +
    layer({panel.abline(h=.95, v=4.199, lty=3); panel.abline(h=c(0, 1), col="gray70")})
## ## there is no dtukey function in R (as of March 2015)
##  HH:::dstrplotContinuous("tukey", args=list(nmeans=4, df=12), X=4.199) +
##    layer(HH:::panel.dstrplotContinuousFill(..., X=4.199))
##   hhdev.off()


###################################################
### code chunk number 21: dstr.tex:446-450
###################################################
## hhcapture("dstr-tukey.Rout", '
  ptukey(4.199, nmeans=4, df=12)
  qtukey(0.95, nmeans=4, df=12)
## ')


###################################################
### code chunk number 22: dstr.tex:479-483
###################################################
##   hhpdf("dstr-t.pdf", height=3, width=6.5)
  HH:::dstrplotContinuous("t", args=list(df=4), X=2) +
    layer(HH:::panel.dstrplotContinuousFill(..., X=2))
##   hhdev.off()


###################################################
### code chunk number 23: dstr.tex:485-490
###################################################
## hhcapture("dstr-t.Rout", '
  dt(2, df=4)
  pt(2, df=4)
  qt(0.9419417, df=4)
## ')


###################################################
### code chunk number 24: dstr.tex:510-514
###################################################
##   hhpdf("dstr-unif.pdf", height=3, width=6.5)
  HH:::dstrplotContinuous("unif", ok0=TRUE, ok1=TRUE, ylim=c(-.05, 1.15), X=.7) +
    layer(HH:::panel.dstrplotContinuousFill(..., X=.7))
##   hhdev.off()


###################################################
### code chunk number 25: dstr.tex:516-521
###################################################
## hhcapture("dstr-unif.Rout", '
  dunif(.7)
  punif(.7)
  qunif(.7)
## ')


###################################################
### code chunk number 26: dstr.tex:532-536
###################################################
##   hhpdf("dstr-weibull.pdf", height=3, width=6.5)
  HH:::dstrplotContinuous("weibull", ok0=TRUE, args=list(shape=4), X=1.3) +
    layer(HH:::panel.dstrplotContinuousFill(..., X=1.3))
##   hhdev.off()


###################################################
### code chunk number 27: dstr.tex:538-543
###################################################
## hhcapture("dstr-weibull.Rout", '
  dweibull(1.3, shape=4)
  pweibull(1.3, shape=4)
  qweibull(0.9425075, shape=4)
## ')


###################################################
### code chunk number 28: dstr.tex:680-687
###################################################
##   hhpdf("dstr-chisq-ncp.pdf", height=3, width=6.5)
  HH:::dstrplotContinuous("chisq", ok0=TRUE, args=list(df=10, ncp=4), ylim=c(-.01, .105),
                          X=18.31,
                          key=HH:::key.dstrplotContinuous()) +
    layer(HH:::panel.dstrplotContinuousFill(..., X=18.31)) +
    HH:::dstrplotContinuous("chisq", ok0=TRUE, args=list(df=10), lty=3, lwd=4, col="deepskyblue3")
##   hhdev.off()


###################################################
### code chunk number 29: dstr.tex:689-694
###################################################
## hhcapture("dstr-chisq-ncp.Rout", '
  dchisq(18.31, df=10, ncp=4)
  pchisq(18.31, df=10, ncp=4)
  qchisq(0.7852264, df=10, ncp=4)
## ')


###################################################
### code chunk number 30: dstr.tex:703-709
###################################################
##   hhpdf("dstr-t-ncp.pdf", height=3, width=6.5)
  HH:::dstrplotContinuous("t", args=list(df=4, ncp=2), X=2, xlim=c(-4, 10), ylim=c(-.02, .4),
                          key=HH:::key.dstrplotContinuous()) +
    layer(HH:::panel.dstrplotContinuousFill(..., X=2)) +
    HH:::dstrplotContinuous("t", args=list(df=4), lty=3, lwd=4, col="deepskyblue3")
##   hhdev.off()


###################################################
### code chunk number 31: dstr.tex:711-716
###################################################
## hhcapture("dstr-t-ncp.Rout", '
  dt(2, df=4, ncp=2)
  pt(2, df=4, ncp=2)
  qt(0.455672, df=4, ncp=2)
## ')


###################################################
### code chunk number 32: dstr.tex:726-732
###################################################
##   hhpdf("dstr-f-ncp.pdf", height=3, width=6.5)
  HH:::dstrplotContinuous("f", ok0=TRUE, args=list(df1=4, df2=20, ncp=2), X=3, ylim=c(-.03, .75),
                          key=HH:::key.dstrplotContinuous()) +
    layer(HH:::panel.dstrplotContinuousFill(..., X=3)) +
    HH:::dstrplotContinuous("f", ok0=TRUE, args=list(df1=4, df2=20), lty=3, lwd=4, col="deepskyblue3")
##   hhdev.off()


###################################################
### code chunk number 33: dstr.tex:734-739
###################################################
## hhcapture("dstr-f-ncp.Rout", '
  df(3, df1=4, df2=20, ncp=2)
  pf(3, df1=4, df2=20, ncp=2)
  qf(0.8710256, df1=4, df2=20, ncp=2)
## ')


###################################################
### code chunk number 34: dstr.tex:796-807
###################################################
## hhcapture("dstr-DiscreteRounding.Rout", '
## this is printing precision, not internal representation
old.digits <- options(digits=7)
ddiscunif(1:6, 6)
pdiscunif(1:6, 6)
qdiscunif(pdiscunif(1:6, 6), 6)
round(pdiscunif(1:6, 6), 4) ## rounded to four decimal digits
## inverse after rounding to four decimal digits
qdiscunif(round(pdiscunif(1:6, 6), 4), 6)
options(old.digits)
## ')


###################################################
### code chunk number 35: dstr.tex:858-862
###################################################
##   hhpdf("dstr-discunif.pdf", height=3, width=6.5)
  HH:::dstrplotDiscrete("discunif", args=list(size=12), x.tick.number=12) +
    layer(HH:::panel.dstrplotDiscreteFill(..., X=6))
##   hhdev.off()


###################################################
### code chunk number 36: dstr.tex:864-869
###################################################
## hhcapture("dstr-discunif.Rout", '
  ddiscunif(6, size=12)
  pdiscunif(6, size=12)
  qdiscunif(.5, size=12)
## ')


###################################################
### code chunk number 37: dstr.tex:895-899
###################################################
##   hhpdf("dstr-binom.pdf", height=3, width=6.5)
  HH:::dstrplotDiscrete("binom", ok0=TRUE, args=list(size=15, prob=.4), x.tick.number=15) +
    layer(HH:::panel.dstrplotDiscreteFill(..., X=6))
##   hhdev.off()


###################################################
### code chunk number 38: dstr.tex:901-911
###################################################
## hhcapture("dstr-binom.Rout", '
## probability of exactly 6 Heads
dbinom(6, size=15, prob=.4)
## probability of 6 or fewer Heads
## extra precision is needed
print(pbinom(6, size=15, prob=.4), digits=17)
## q, the number for which the probability of seeing
## q or fewer Heads is 0.60981315570892769
qbinom(0.60981315570892769, size=15, prob=.4)
## ')


###################################################
### code chunk number 39: dstr.tex:924-928
###################################################
##   hhpdf("dstr-geom.pdf", height=3, width=6.5)
  HH:::dstrplotDiscrete("geom", ok0=TRUE, args=list(prob=.5), size=9, x.tick.number=9) +
    layer(HH:::panel.dstrplotDiscreteFill(..., X=3))
##   hhdev.off()


###################################################
### code chunk number 40: dstr.tex:930-935
###################################################
## hhcapture("dstr-geom.Rout", '
  dgeom(3, prob=.5)
  pgeom(3, prob=.5)
  qgeom(0.9375, prob=.5)
## ')


###################################################
### code chunk number 41: dstr.tex:952-956
###################################################
##   hhpdf("dstr-hyper.pdf", height=3, width=6.5)
  HH:::dstrplotDiscrete("hyper", ok0=TRUE, args=list(m=5, n=6, k=7), size=7) +
    layer(HH:::panel.dstrplotDiscreteFill(..., X=3))
##   hhdev.off()


###################################################
### code chunk number 42: dstr.tex:958-963
###################################################
## hhcapture("dstr-hyper.Rout", '
  dhyper(3, m=5, n=6, k=7)
  print(phyper(3, m=5, n=6, k=7), digits=17)
  qhyper(0.65151515151515149, m=5, n=6, k=7)
## ')


###################################################
### code chunk number 43: dstr.tex:978-982
###################################################
##   hhpdf("dstr-nbinom.pdf", height=3, width=6.5)
  HH:::dstrplotDiscrete("nbinom", ok0=TRUE, args=list(size=8.5, prob=.4), size=30) +
    layer(HH:::panel.dstrplotDiscreteFill(..., X=17))
##   hhdev.off()


###################################################
### code chunk number 44: dstr.tex:984-989
###################################################
## hhcapture("dstr-nbinom.Rout", '
  dnbinom(17, size=8.5, prob=.4)
  print(pnbinom(17, size=8.5, prob=.4), digits=17)
  qnbinom(0.81209497223034977, size=8.5, prob=.4)
## ')


###################################################
### code chunk number 45: dstr.tex:1010-1014
###################################################
##   hhpdf("dstr-pois.pdf", height=3, width=6.5)
  HH:::dstrplotDiscrete("pois", ok0=TRUE, args=list(lambda=4), size=12) +
    layer(HH:::panel.dstrplotDiscreteFill(..., X=7))
##   hhdev.off()


###################################################
### code chunk number 46: dstr.tex:1016-1021
###################################################
## hhcapture("dstr-pois.Rout", '
  dpois(7, lambda=4)
  print(ppois(7, lambda=4), digits=17)
  qpois(0.94886638420715264, lambda=4)
## ')


###################################################
### code chunk number 47: dstr.tex:1040-1044
###################################################
##   hhpdf("dstr-signrank.pdf", height=3, width=6.5)
  HH:::dstrplotDiscrete("signrank", ok0=TRUE, args=list(n=5), size=5*(5+1)/2, x.tick.number=15) +
    layer(HH:::panel.dstrplotDiscreteFill(..., X=11))
##   hhdev.off()


###################################################
### code chunk number 48: dstr.tex:1046-1051
###################################################
## hhcapture("dstr-signrank.Rout", '
  dsignrank(11, n=5)
  psignrank(11, n=5)
  qsignrank(0.84375, n=5)
## ')


###################################################
### code chunk number 49: dstr.tex:1060-1064
###################################################
##   hhpdf("dstr-wilcox.pdf", height=3, width=6.5)
  HH:::dstrplotDiscrete("wilcox", args=list(m=4, n=12), size=4*12, x.tick.number=10) +
    layer(HH:::panel.dstrplotDiscreteFill(..., X=35))
##   hhdev.off()


###################################################
### code chunk number 50: dstr.tex:1066-1071
###################################################
## hhcapture("dstr-wilcox.Rout", '
  dwilcox(35, m=4, n=12)
  print(pwilcox(35, m=4, n=12), digits=17)
  qwilcox(0.9148351648351648, m=4, n=12)
## ')


###################################################
### code chunk number 51: dstr.tex:1084-1125
###################################################
old.width=options(width=70)
## hhcapture("dstr-multinom.Rout", '
## This example is based on ?dmultinom in R
## all possible outcomes of Multinom(N = 3, K = 3)
X <- t(as.matrix(expand.grid(0:3, 0:3)))
X <- X[, colSums(X) <= 3]
X <- rbind(X, 3:3 - colSums(X))
dimnames(X) <- list(letters[1:3], apply(X, 2, paste, collapse=""))
Y <- round(apply(X, 2, function(x) dmultinom(x, prob = c(1,2,5))), 3)
rbind(X, Y)
## ')
options(old.width)

xy <- rbind(x=-c(0.5, 0, 1),
            R=c(1.0, 0, 0))

xRBY <- rbind(c(3,0) + (xy %*% X), B=X[2,], Y)
xRBY

multinomialPlot <-
xyplot((R+.10) ~ (x+.3), data=data.frame(t(xRBY)),
       panel=panel.text, labels=xRBY["Y",], cex=.8,
       xlim=c(-.25, 3.55), ylim=c(-.5, 3.5),
       main="dmultinom(x, prob = c(1,2,5))") +
xyplot(R-.10 ~ x, data=data.frame(t(xRBY)),
       panel=panel.text, labels=colnames(xRBY)) +
xyplot(R ~ x, data=data.frame(t(xRBY)),
       panel=function(x, y, height, ...)
         panel.rect(x-.02, y, x+.02, y+height, col="blue", border="blue", ...),
        height=2.5*xRBY["Y",]) +
layer(under=TRUE, {
   panel.segments(c(0, .5, 1, 1.5), c(0, 1, 2, 3), c(3, 2.5, 2, 1.5), c(0, 1, 2, 3), col="gray30", lty=3, lwd=1.2) ## horizontal
   panel.segments(c(0, 1, 2, 3),    c(0, 0, 0, 0), c(1.5, 2, 2.5, 3), c(3, 2, 1, 0), col="gray30", lty=3, lwd=1.2) ## SW-NE
   panel.segments(c(0, .5, 1, 1.5), c(0, 1, 2, 3), c(0, 1, 2, 3),     c(0, 0, 0, 0), col="gray30", lty=3, lwd=1.2) ## NW-SE
})
multinomialPlot

## hhpdf("dstr-multinom.pdf", height=4, width=7)
update(multinomialPlot,
       xlab=NULL, ylab=NULL, scales=list(col="transparent", col.line="transparent"))
## hhdev.off()


###################################################
### code chunk number 52: dstr.tex:1153-1165
###################################################
##   hhpdf("dstr-mvnorm.pdf", height=7, width=7)
  z <- matrix(nrow=81, ncol=81)
  seqij <- seq(-2, 2, .05)
  for (i in 1:81) for (j in 1:81)
  z[i,j] <- dmvnorm(seqij[c(i, j)])
  zblue <- z
  zblue[] <- NA
  zblue[1:21, 1:21] <- z[1:21, 1:21]
  wireframe(z, row.values=seqij, column.values=seqij, zlim=c(0, max(z))) +
  wireframe(zblue, row.values=seqij, column.values=seqij, zlim=c(0, max(z)), col="blue", zlab=NULL)
  grid.text("Bivariate Normal: dmvnorm(c(-1, -1))", x=.5, y=.9, gp=gpar(fontface="bold", cex=1.2))
##   hhdev.off()


###################################################
### code chunk number 53: dstr.tex:1167-1172
###################################################
## hhcapture("dstr-mvnorm.Rout", '
dmvnorm(c(-1, -1))
pmvnorm(upper=c(-1, -1))[1]
qmvnorm(0.02517, mean=c(0,0))$quantile
## ')



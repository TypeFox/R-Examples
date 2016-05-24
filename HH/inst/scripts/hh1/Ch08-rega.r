library(HH)

## ##  hh/rega/code/

##   53 2002-08-01  readfat.s  ## read data
data(fat)

##  111 2004-04-26  rega.f1.s  ## splom of all data Figure 8.1
## follows rega/code/readfat.s

splom( ~ fat, main="Fat data", pch=16)
## export.eps(hh("rega/figure/f1.eps"))

##  101 2004-03-29  ls.s       ## regression analysis Table 8.1
## least-squares fit
fat.lm <- lm(bodyfat ~ abdomin, data=fat)
summary(fat.lm, corr=FALSE)
anova(fat.lm)

##  209 2004-04-26  rega.f6.s  ## diagnostic plots Figure 8.5
old.par <- par(pch=16, mar=c(5,6,4,2)+.1)
par(mfrow=c(2,3))
plot(fat.lm, cex=.8, which=1:6)
mtext("diagnostics for lm(bodyfat ~ abdomin, data=fat)", outer=TRUE, cex=1.4)
par(old.par)
## export.eps(hh("rega/figure/f6.eps"))

##  182 2004-02-15  fat-ci.s   ## confidence and prediction intervals Figure 8.4
## follows readfat.s. ls.s

tmp <- ci.plot(fat.lm, xlab=list(cex=1.4), ylab=list(cex=1.4), main.cex=1.4)
print(tmp, position=c(0,0,1,.8))
## export.eps(hh("rega/figure/fat-ci.eps"))


## ## Figure 8.2 and relatives

## 2058 2004-03-29  residuals.s     ## Figure 8.2
## follows from
##   rega/code/readfat.s

## Explain linear regression and residuals, HH section 8.2.1
## This is done as a mockup in regular graphics

old.par <- par(mfrow=c(2,3), oma=c(3,2,4,9), mar=c(4,5,1,1)+.1, pch=16)
frame()

par(mfg=c(1,1, 2,3))
plot(fat[,"abdomin"], fat[,"bodyfat"],
     ylab="", xlab="",
     xlim=c(70,195),ylim=c(0,50))
axis(3,labels=FALSE)
axis(4,labels=FALSE)
mtext("abdomen", line=3, cex=1.2, side=1)
if.R(r=text("bodyfat", x=10, y=40, cex=1.8, srt=90, xpd=NA),
     s=text("bodyfat", x=20, y=40, cex=1.2, srt=90))


## draw the residuals
par(mfg=c(1,2, 2,3))
regr1.plot(fat[,"abdomin"], fat[,"bodyfat"],
           xaxt="n", yaxt="n",
     ylab="", xlab="",
           resid.plot="line", main="",
           xlim=c(70,195), ylim=c(0,50))
axis(1,labels=FALSE)
axis(2,labels=FALSE)
axis(3,labels=FALSE)
axis(4,labels=FALSE)
mtext("least-squares fit:", line=4, cex=1.2)
mtext("y = -28.56 + .505 x", line=1.8, cex=1.2)

par(mfg=c(1,3, 2,3))
regr1.plot(fat[,"abdomin"], fat[,"bodyfat"],
           coef=c(20, .1),
           xaxt="n", yaxt="n",
     ylab="", xlab="",
           resid.plot="line", main="",
           xlim=c(70,195), ylim=c(0,50))
axis(1,labels=FALSE)
axis(2,labels=FALSE)
axis(3,labels=FALSE)
axis(4)
mtext("too shallow:", line=4, cex=1.2)
mtext("y = 20 + .1x", line=1.8, cex=1.2)
if.R(r=text("residuals", x=300, y=40, cex=1.8, xpd=NA),
     s=text("residuals", x=270, y=40, cex=1.2))

## square the residuals
par(mfg=c(2,2, 2,3))
regr1.plot(fat[,"abdomin"], fat[,"bodyfat"],
           xaxt="n", yaxt="n",
     ylab="", xlab="",
           resid.plot="square", main="",
           xlim=c(70,195), ylim=c(0,50))
axis(1)
axis(2)
axis(3,labels=FALSE)
axis(4,labels=FALSE)
if.R(r=text("bodyfat", x=10, y=40, cex=1.8, srt=90, xpd=NA),
     s=text("bodyfat", x=20, y=40, cex=1.2, srt=90))
mtext("abdomen", line=3, cex=1.2, side=1)

par(mfg=c(2,3, 2,3))
regr1.plot(fat[,"abdomin"], fat[,"bodyfat"],
           coef=c(20, .1),
     ylab="", xlab="",
           xaxt="n", yaxt="n",
           resid.plot="square", main="",
           xlim=c(70,195), ylim=c(0,50))
axis(1)
axis(2,labels=FALSE)
axis(3,labels=FALSE)
axis(4)
if.R(r=text("squared\nresiduals", x=300, y=40, cex=1.8, xpd=NA),
     s=text("squared\nresiduals", x=270, y=40, cex=1.2))
mtext("abdomen", line=3, cex=1.2, side=1)

par(old.par)
## export.eps(hh("rega/figure/resid2x2.eps"))

##  557 2004-03-29  rega.f5.s       ## Figure 8.3
old.par <- par(mfrow=c(1,2), mar=par("mar")+c(0,2,3,0))

regr1.plot(fat[,"abdomin"], fat[,"bodyfat"],
              main="Variance about\nmean of y",
           xlab="abdomin", ylab="bodyfat",
      coef=c(mean(fat[,"bodyfat"]), 0),
      resid.plot="line",
      ylim=c(0,45), cex=1.4, pch=16)

regr1.plot(fat[,"abdomin"], fat[,"bodyfat"],
              main="Variance about\nleast-squares line",
           xlab="abdomin", ylab="bodyfat",
      resid.plot="line",
      ylim=c(0,45), cex=1.4, pch=16)

par(old.par)
## export.eps(hh("rega/figure/f5.eps"))



## ## Figure 8.6
##  626 2004-04-26  diagnostic.s ## Figure 8.6
## Explanation of panel 5 of
##      plot(fat.lm)
## where
##      fat.lm <- lm(bodyfat ~ abdomin, data=fat)

## follows files
## rega/code/fat-ci.s
## rega/code/ls.s

old1.par <- par(mfrow=c(1,2), mar=c(6,6,10,2)+.1, mgp=c(4,1,0), xpd=NA)
old2.par <- par(pch=16)

old3.par <- par(cex=1.5)
##   89 2002-11-10  diag1.s      ## 8.6a
pp <- ppoints(101)
x <- qnorm(pp)
plot(pp ~ x, main="Cumulative Distribution of N(0,1)")
## export.eps(hh.file("rega/figure/diag1.eps"))
if.R(r=frame(),
     s={frame();frame()})
##  101 2002-12-17  diag2f.s     ## 8.6b fitted values
n <- nrow(fat)
f <- (1:n)/n
plot(f ~ sort(predict(fat.lm)))
title("empirical cdf of\nfitted values")

##   67 2002-12-17  diag2r.s     ## 8.6b residuals
plot(f ~ sort(resid(fat.lm)))
title("empirical cdf of\nresiduals")
## export.eps(hh.file("rega/figure/diag2.eps"))


##  360 2002-12-17  diag3.s      ## 8.6c
bodyfat.range <- range(predict(fat.lm) - mean(predict(fat.lm)),
                       resid(fat.lm)                          )
plot(f ~ sort(predict(fat.lm) - mean(predict(fat.lm))),
     xlim=bodyfat.range,
     main="empirical cdf of\ncentered fitted values")
plot(f ~ sort(resid(fat.lm)),
     xlim=bodyfat.range,
     main="empirical cdf of\n residuals")
## export.eps(hh.file("rega/figure/diag3.eps"))

##  381 2002-12-17  diag4.s      ## 8.6d
bodyfat.range <- range(predict(fat.lm) - mean(predict(fat.lm)),
                       resid(fat.lm)                          )
plot(sort(predict(fat.lm) - mean(predict(fat.lm))) ~ f,
     ylim=bodyfat.range,
     main="transposed empirical cdf\nof centered fitted values")
plot(sort(resid(fat.lm)) ~ f,
     ylim=bodyfat.range,
     main="transposed empirical\ncdf of residuals")
## export.eps(hh.file("rega/figure/diag4.eps"))

par(old3.par)
par(old2.par)
par(old1.par)

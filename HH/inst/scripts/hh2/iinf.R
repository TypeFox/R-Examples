### R code from vignette source '~/WindowsC/HOME/rmh/hh.e2/hh2/iinf.tex'

###################################################
### code chunk number 1: iinf.tex:9-10
###################################################
library(HH)


###################################################
### code chunk number 2: iinf.tex:13-18
###################################################
## the standard lattice color 2 is difficult for people with color deficient vision
data(col3x2)
## These colors look like a 3x2 color array when run through
## the vischeck simulator to see how they look for the three most
## common color vision deficiencies: Protanope, Deuteranope, Tritanope.


###################################################
### code chunk number 3: iinf.tex:167-192
###################################################
## hhpdf("Hyz6.pdf", height=8.5, width=7)

NTfun <- function(...)
   NTplot(..., main=NULL, xlab=NULL, ylab=NULL, prob.labels=FALSE,
          xhalf.multiplier=1.6, yhalf.multiplier=2, cex.top.axis=.7,
          scales=list(cex=.5, y=list(at=1)),
          xlim=c(-3, 4))

Nhr  <- NTfun(mean0=0, xbar= 1.8, alpha.right=.05,  alpha.left=0)
Nhl  <- NTfun(mean0=0, xbar=-1.8, alpha.right=0,    alpha.left=.05)
Nhlr <- NTfun(mean0=0, xbar= 1.8, alpha.right=.025, alpha.left=.025)

Ncl  <- NTfun(mean0=0, xbar= 1.8, alpha.right=0,    alpha.left=.05,  type="confidence")
Ncr  <- NTfun(mean0=0, xbar=-1.8, alpha.right=.05,  alpha.left=0,    type="confidence")
Nclr <- NTfun(mean0=0, xbar= 1.8, alpha.right=.025, alpha.left=.025, type="confidence")

print(Nhr,  tablesOnPlot=FALSE, position=c(.00, .630, .55, 1.000), more=TRUE)
print(Nhl,  tablesOnPlot=FALSE, position=c(.00, .315, .55,  .685), more=TRUE)
print(Nhlr, tablesOnPlot=FALSE, position=c(.00, .000, .55,  .370), more=TRUE)

print(Ncl,  tablesOnPlot=FALSE, position=c(.45, .630, 1.0, 1.000), more=TRUE)
print(Ncr,  tablesOnPlot=FALSE, position=c(.45, .315, 1.0,  .685), more=TRUE)
print(Nclr, tablesOnPlot=FALSE, position=c(.45, .000, 1.0,  .370), more=FALSE)

## hhdev.off()


###################################################
### code chunk number 4: iinf.tex:427-452
###################################################
## hhpdf("Hyp6.pdf", height=8.5, width=7)

NTfun <- function(...)
   NTplot(..., main=NULL, xlab=NULL, ylab=NULL, prob.labels=FALSE,
          xhalf.multiplier=1.6, yhalf.multiplier=2, cex.top.axis=.7,
          scales=list(cex=.5),
          xlim=c(0, 1), ylim=c(0, 4), distribution.name="binomial")

Nphr  <- NTfun(p0=.5, p.hat= .75, n=20, alpha.right=.05,  alpha.left=0)
Nphl  <- NTfun(p0=.5, p.hat= .25, n=20, alpha.right=0,    alpha.left=.05)
Nphlr <- NTfun(p0=.5, p.hat= .75, n=20, alpha.right=.025, alpha.left=.025)

Npcl  <- NTfun(p0=.5, p.hat= .75, n=20, alpha.right=0,    alpha.left=.05,  type="confidence")
Npcr  <- NTfun(p0=.5, p.hat= .25, n=20, alpha.right=.05,  alpha.left=0,    type="confidence")
Npclr <- NTfun(p0=.5, p.hat= .75, n=20, alpha.right=.025, alpha.left=.025, type="confidence")

print(Nphr,  tablesOnPlot=FALSE, position=c(.00, .630, .55, 1.000), more=TRUE)
print(Nphl,  tablesOnPlot=FALSE, position=c(.00, .315, .55,  .685), more=TRUE)
print(Nphlr, tablesOnPlot=FALSE, position=c(.00, .000, .55,  .370), more=TRUE)

print(Npcl,  tablesOnPlot=FALSE, position=c(.45, .630, 1.0, 1.000), more=TRUE)
print(Npcr,  tablesOnPlot=FALSE, position=c(.45, .315, 1.0,  .685), more=TRUE)
print(Npclr, tablesOnPlot=FALSE, position=c(.45, .000, 1.0,  .370), more=FALSE)

## hhdev.off()


###################################################
### code chunk number 5: iinf.tex:568-572
###################################################
data(vocab)
## hhcapture("vocab-stem.Rout", '
stem(vocab$score, scale=2)
## ')


###################################################
### code chunk number 6: iinf.tex:586-593
###################################################
## hhcapture("vocab-t.Rout", '
vocab.t <- t.test(vocab$score, mu=10)
vocab.t
## ')
## hhpdf("vocab.pdf", height=6, width=7)
NTplot(vocab.t, type="confidence", cex.prob=1)
## hhdev.off()


###################################################
### code chunk number 7: iinf.tex:656-677
###################################################
## hhpdf("chisqCI.pdf", height=4, width=7)
s2 <- 15
df <- 12
qch <- qchisq(c(.025, .975), df=12)
old.omd <- par(omd=c(.05,.88, .10,1))
chisq.setup(df=12, ylab="Density")
chisq.curve(df=12, col='blue', alpha=c(.025, .025), axis.name=expression(chi^2))
(df*s2)/rev(qch)
qs <- c(.1, .5, 1, 2, 3)
axis(1, at=qch, labels=round(12/qch, 2), line=3.5, lwd=0, lwd.ticks=0)
axis(1, at=-5, line=3.5, labels=expression(nu ~ "/" ~ chi^2), xpd=TRUE, lwd=0, lwd.ticks=0)
axis(1, at=qch, labels=round(s2*12/qch, 2), line=5, lwd=0, lwd.ticks=0)
axis(1, at=-5, line=5, labels=expression(s^2 ~ nu ~ "/" ~ chi^2), xpd=TRUE, lwd=0, lwd.ticks=0)
## axis(1, at=12/qs, labels=qs, line=5, lwd=0, lwd.ticks=0)
## axis(1, at=-5, line=5, labels=expression(nu ~ "/" ~ chi^2), xpd=TRUE, lwd=0, lwd.ticks=0)
## chisq.observed(15, df=12, axis.name=expression(chi^2)) ## not part of this example
axis(3, at=13, line=-.5, labels=expression(s^2==15), xpd=TRUE, lwd=0, lwd.ticks=0)
par(old.omd)
rev(s2*12/qch)
sqrt(rev((s2*12/qch)))
## hhdev.off()


###################################################
### code chunk number 8: iinf.tex:911-956
###################################################
## hhcapture("ttest2.Rout", '
data(cereals)
table(cereals[,c("mfr","type")])
C.KG <- cereals$type=="C" & cereals$mfr %in% c("K","G")
cerealsC <- cereals[C.KG, c("mfr", "carbo") ]
cerealsC$mfr <- factor(cerealsC$mfr)
bwplot(carbo ~ mfr, data=cerealsC) +
dotplot(carbo ~ mfr, data=cerealsC)
t.t <- t.test(carbo ~ mfr, data=cerealsC, var.equal=TRUE)
t.t
## ')
## hhcapture("ttest2B.Rout", '
mm <- tapply(cerealsC$carbo, cerealsC$mfr, mean)
vv <- tapply(cerealsC$carbo, cerealsC$mfr, var)
ll <- tapply(cerealsC$carbo, cerealsC$mfr, length)
s2p <- ((ll-1) %*% vv) / sum(ll-1)
tt <- -diff(mm) / (sqrt(s2p) * sqrt(sum(1/ll)))
tt
## ')
NTplot(t.t, zaxis=TRUE) ## This is sufficient for working on the screen.
## hhpdf("ttest2.pdf", height=6, width=12)
## The rest of the arguments adjust fonts and placement and colors for an extra-wide window
## and prevent overprinting of tick labels.

printbook.colors <- c(  ## based on "original" colors
  col.alpha             = "blue",
  col.notalpha          = "lightcyan",  ## this value is nonstandard
  col.beta              = "red",
  col.power             = "pink",
  col.pvalue            = "green",
  col.pvaluetranslucent = HH:::ColorWithAlpha("springgreen"),  ## this value is nonstandard
  col.critical          = "gray50",
  col.border            = HH:::ColorWithAlpha("black"),
  col.text              = "black",
  col.conf              = "lightgreen"
)

NT <- NTplot(t.t, zaxis=TRUE, cex.z=1, xlim=c(-2.7, 2.5),
             ntcolors=printbook.colors, ## this line uses nonstandard colors for the printed book
             cex.top.axis=1.1, digits=3,
             cex.main=1.5, cex.prob=1.2, key.axis.padding=5,
             xhalf.multiplier=.35, yhalf.multiplier=1.2)

print(update(NT, scales=list(cex=1.2)), cex.table=1.2)
## hhdev.off()


###################################################
### code chunk number 9: iinf.tex:1124-1204
###################################################
data(teachers)
teachers$"English-Greek" <- teachers$English - teachers$Greek
teachers$sentence <- ordered(1:32, order(teachers$"English-Greek"))
teachers$smaller <- factor(sign(teachers$"English-Greek"),
                           labels=c("English\nsmaller error count",
                           "Greek\nsmaller error count"))
teachers[c(1,2,32),]

teachers.melt <- reshape2:::melt(cbind(teachers, sentence=1:32),
                                 id.vars=c("sentence","English-Greek", "smaller"),
                                 variable.name="language",
                                 value.name="errors")
teachers.melt[c(1,2,32,33,63,64),]

## hhpdf("teachers-dot.pdf", height=6.9, width=9)
trellis.par.set(  ## must be inside device
   superpose.symbol=
       list(cex=rep(1.8, 7),
            col=c(likertColor(2,colorFunctionOption ="default"),3:7)),
   plot.symbol=list(cex=rep(1.1, 7)))

tpg <- trellis.par.get("superpose.symbol")
tpg$pch[1:3] <- c(69,71,17) ## E G triangle

AAAA <-
resizePanels(h=table(teachers$smaller),
dotplot(rep(teachers$sentence, 2) ~ errors | smaller,
        groups=language, data=teachers.melt,
        scales=list(y=list(relation="free")), pch=tpg$pch,
        par.strip.text=list(cex=1.1, lines=2),
        ylab="sentence",
        layout=c(1,2), between=list(y=1),
        key=list(text=
                   list(c("English for\nGreek Speakers",
                          "Greek for\nEnglish Speakers")),
          points=Rows(tpg, 1:2), padding.text=3,
          space="right",  border=TRUE)
        ))
## AAAA

BBBB <-
resizePanels(h=table(teachers$smaller),
             dotplot(sentence ~ `English-Greek` | smaller,
                     data=teachers, ref=0, pch=tpg$pch[3], col="gray30",
                     scales=list(y=list(relation="free")),
                     par.strip.text=list(cex=1.1, lines=2),
                     main="English - Greek",
                     layout=c(1,2), between=list(y=1)
                     )
             )
## BBBB + layer(panel.abline(v=0, col="gray60"))

CCCC <-
resizePanels(h=table(teachers$smaller),
             dotplot(sentence ~ sqrt(`English-Greek` + 17) | smaller,
                     data=teachers, ref=0, pch=tpg$pch[3], col="gray30",
                     par.strip.text=list(cex=1.1, lines=2),
                     scales=list(y=list(relation="free")),
                     main=expression(sqrt("English - Greek" + 17)),
                     layout=c(1,2), between=list(y=1))
                     )
## CCCC + layer(panel.abline(v=sqrt(17), col="gray60"))

AAAABBBBCCCC <-
resizePanels(h=table(teachers$smaller), w=c(30, 25, 25),
combineLimits(
update(cbind(errors=AAAA, "English-Greek"=BBBB, "sqrt(17 +\nEnglish-Greek)"=CCCC),
       xlab=NULL, ylab=list("sentence", cex=1.2),
       par.strip.text=list(cex=1.1, lines=1.2),
       scales=list(
         relation="free",
         x=list(limits=list(c(9, 42), c(-20, 20), c(0,6))))
       )))

xlab.top <- c(AAAABBBBCCCC$condlevels[[1]][1:2], expression(sqrt("English-Greek + 17")~" "))
AAAABBBBCCCC$condlevels[[1]] <- c(" "," "," ")
AAAABBBBCCCC + layer(panel.abline(v=c(0,0,sqrt(17))[current.column()], col="gray60"))
grid.text(xlab.top, x=c(.25, .475, .686), y=.92)

## hhdev.off()


###################################################
### code chunk number 10: iinf.tex:1251-1259
###################################################
## hhcapture("teachersLeft.Rout", '
stem(teachers$"English-Greek")
t.test(teachers$"English-Greek")
## ')
## hhcapture("teachersRight.Rout", '
stem(sqrt(teachers$"English-Greek" + 17), scale=.5)
t.test(sqrt(teachers$"English-Greek" + 17), mu=sqrt(17))
## ')


###################################################
### code chunk number 11: iinf.tex:1367-1392
###################################################
## hhpdf("normalOneTailFig.pdf", height=7.2, width=9)
mean0 <- 1
delta <- 2
mean1 <- mean0 + delta
sd <- 3
alpha <- .05
beta <- .20
power <- 1-beta
n <- sd^2 * (qnorm(1-alpha) + qnorm(1-beta))^2 / delta^2
n
tmp <-
NTplot(mean0=mean0, mean1=mean1, sd=sd, n=n,
        key.axis.padding=6, cex.main=1.5,
        cex.top.axis=1.8, cex.prob=1.5,
        prob.labels=FALSE,
        digits.axis=5, digits.float=3, digits.left=3,
        xhalf.multiplier=.8, yhalf.multiplier=.6,
        zaxis=TRUE, z1axis=TRUE, cex.z=1)
print(tmp +
        layer({panel.segments(1, -.080, 3.000, -.080, col="gray60",    lwd=6) ## xbar scale
               panel.segments(1, -.098, 2.323, -.098, col="lightblue", lwd=6) ## z scale
               panel.segments(3, -.122, 2.323, -.122, col="#FF808B",   lwd=6) ## z1 scale
              }, under=TRUE),
      cex.table=1.2, digits=5)
## hhdev.off()


###################################################
### code chunk number 12: iinf.tex:1470-1483
###################################################
## hhcapture("sampleSizeA.Rout", '
## one sided
alpha <- .05
power <- .80
beta <- 1-power
delta <- 1
sd <- 2

## Approximation using formula assuming normal is appropriate
sd^2*(qnorm(1-alpha) + qnorm(1-beta))^2 / delta^2
## [1] 24.73
## n is slightly smaller with the normal assumption.
## ')


###################################################
### code chunk number 13: iinf.tex:1499-1508
###################################################
## hhcapture("sampleSizeB.Rout", '
## solve using power.t.test
PTT <-
power.t.test(delta=delta, sd=sd, sig.level=alpha, power=power,
             type="one.sample", alternative="one.sided")
PTT
NTplot(PTT, zaxis=TRUE)  ## static plot
## NTplot(PTT, zaxis=TRUE, shiny=TRUE)  ## dynamic plot
## ')


###################################################
### code chunk number 14: iinf.tex:1525-1551
###################################################
## hhcapture("sampleSizeC.Rout", '
## solve manually with t distribution.  Use ncp for alternative.
n0 <- 30 ## pick an n0 for starting value
t.critical <- qt(1-alpha, df=n0-1)
t.critical
## [1] 1.699

## a series of n values
nn <- 23:30
names(nn) <- nn
nn

## find the power for a series of n values for the specified critical value
pt(t.critical, df=nn-1, ncp=delta/(sd/sqrt(nn)), lower=FALSE)
##     23     24     25     26     27     28     29     30
## 0.7568 0.7722 0.7868 0.8006 0.8136 0.8258 0.8374 0.8483

## recalculate critical value with new n=26
t.critical <- qt(1-alpha, df=26-1)
t.critical
## find the power for a series of n values for the new critical value
pt(t.critical, df=nn-1, ncp=delta/(sd/sqrt(nn)), lower=FALSE)
##     23     24     25     26     27     28     29     30
## 0.7540 0.7695 0.7842 0.7981 0.8112 0.8235 0.8352 0.8461
## conclude n between 26 and 27
## ')


###################################################
### code chunk number 15: iinf.tex:1609-1635
###################################################
tmp26 <- NTplot(NTplot(PTT), main=NA,
                key.axis.padding=6, cex.main=1.5,
                cex.top.axis=1.8, cex.prob=1.5,
                prob.labels=FALSE,
                digits.axis=5, digits.float=3, digits.left=3,
                xhalf.multiplier=.8, yhalf.multiplier=.6,
                zaxis=TRUE, cex.z=1)
## tmp26

tmp30 <- NTplot(NTplot(PTT), n=30, df=29, main=NA,
                key.axis.padding=6, cex.main=1.5,
                cex.top.axis=1.8, cex.prob=1.5,
                prob.labels=FALSE,
                digits.axis=5, digits.float=3, digits.left=3,
                xhalf.multiplier=.8, yhalf.multiplier=.6,
                zaxis=TRUE, cex.z=1)
## tmp30

## hhpdf("sampleSizeFig26.pdf", height=8, width=6) ## density, power, beta by sample size
print(update(tmp26, ylab=NULL), cex.table=1.2)
## hhdev.off()

## hhpdf("sampleSizeFig30.pdf", height=8, width=6) ## density, power, beta by sample size
print(update(tmp30, ylab=NULL), cex.table=1.2)
## hhdev.off()



###################################################
### code chunk number 16: iinf.tex:1801-1807
###################################################
## hhcapture("fairDice.Rout", '
dice <- sample(rep(1:6, c(3,7,5,8,1,6)))
dice
table(dice)
chisq.test(table(dice))
## ')


###################################################
### code chunk number 17: iinf.tex:1818-1825
###################################################
## hhpdf("fairDice.pdf", height=4, width=7)
old.omd <- par(omd=c(.05, .88, .05, 1))
chisq.setup(df=5, ylab="Density")
chisq.curve(df=5, col='blue', axis.name=expression(chi^2))
chisq.observed(6.8, df=5, axis.name=expression(chi^2))
par(old.omd)
## hhdev.off()


###################################################
### code chunk number 18: iinf.tex:1875-1896
###################################################
## hhpdf("family.pdf", height=3, width=7)
Observed <- c(13, 18, 20, 18, 6, 5)
names(Observed) <- 0:5
## binomial proportion p=.4 is specified
Expected <- dbinom(0:5, size=5, p=.4)*80
names(Expected) <- 0:5
update(
c(
Observed = xyplot(Observed ~ 0:5,
                  horizontal=FALSE, origin=0, col="#727EB5",
                  panel=panel.barchart, ylim=c(-1, 28)),
"Binomial p=.4" = xyplot(Expected ~ 0:5,
                  horizontal=FALSE, origin=0, col="#727EB5",
                  panel=panel.barchart, ylim=c(-1, 28))
),
between=list(x=1), scales=list(alternating=FALSE), ylab=NULL)
## hhdev.off()

var(rep(0:5, times=Observed))

var(rep(0:5, times=Expected))


###################################################
### code chunk number 19: iinf.tex:1914-1934
###################################################
## hhcapture("ex.chisq.Rout", '
Observed <- c(13, 18, 20, 18, 6, 5)
names(Observed) <- 0:5
## binomial proportion p=.4 is specified
Expected <- dbinom(0:5, size=5, p=.4)*80
names(Expected) <- 0:5
chisq.test(Observed, p=Expected, rescale.p=TRUE)


## binomial proportion p is calculated from the observations
p <- sum(Observed * (0:5)/5)/sum(Observed)
p
Expected <- dbinom(0:5, size=5, p=p)*80
names(Expected) <- 0:5
WrongDF <- chisq.test(Observed, p=Expected, rescale.p=TRUE)
WrongDF
c(WrongDF$statistic, WrongDF$parameter)
## correct df and p-value
pchisq(WrongDF$statistic, df=WrongDF$parameter - 1, lower=FALSE)
## ')


###################################################
### code chunk number 20: iinf.tex:2088-2109
###################################################
## hhpdf("iinf-f-qqnorm.pdf", height=7.2, width=9)

qqnorm.dat <- data.frame(
   normal=rnorm(100,0,1),
   uniform=runif(100,-3,3))
qqnorm.dat <- within(qqnorm.dat, {
  `heavy-tailed` <- normal/uniform
  `thin-tailed` <- normal*exp(-abs(normal))
  `positively skewed` <- normal^2
  `negatively skewed` <- -`positively skewed`
  quantiles <- qnorm(ppoints(length(normal)))
})
qqnorm.dat <- data.frame(lapply(qqnorm.dat, sort), check.names=FALSE) ## sort on quantile is redundant
head(qqnorm.dat)

xyplot(normal + uniform + `negatively skewed` + `positively skewed` +
      `thin-tailed` + `heavy-tailed` ~ quantiles, data=qqnorm.dat, outer=TRUE,
      ylab="Observed values", xlab="Quantiles of Standard Normal",
      col=likertColor(2)[2],
      scales=list(rot=0, alternating=FALSE, y=list(relation="free")), between=list(x=1, y=1), layout=c(2,3))
## hhdev.off()


###################################################
### code chunk number 21: iinf.tex:2120-2124
###################################################
## hhpdf("iinf-f-dist.pdf", height=7.2, width=9)
tmp <- lapply(qqnorm.dat, histogram, nint=15, col=likertColor(2)[2])
update(do.call(latticeExtra:::c.trellis, c(tmp[-3], list(layout=c(2,3)))), xlab="Observed Values")
## hhdev.off()


###################################################
### code chunk number 22: iinf.tex:2142-2163
###################################################
## hhpdf("iinf-f1.pdf", height=4, width=7)
points <- ppoints(100)
qqt.data <- data.frame(
  ppoints=points,
  q3=qt(points, 3),
  q5=qt(points, 5),
  q7=qt(points, 7),
  qn=qnorm(points))
qqt.data <- within(qqt.data, {
  d3 <- dt(q3, 3)
  d5 <- dt(q5, 5)
  d7 <- dt(q7, 7)
  dn <- dnorm(qn)
})
head(qqt.data)
xyplot(q3+q5+q7+qn ~ qn, data=qqt.data, type="l", aspect="iso", ylab="distribution",
       par.settings=list(superpose.line=list(lty=c(5,4,2,1))),
       auto.key=list(space="right", border=TRUE,
                     text=c("t, 3 df", "t, 5 df", "t, 7 df", "normal"),
                     lines=TRUE, points=FALSE))
## hhdev.off()


###################################################
### code chunk number 23: iinf.tex:2184-2196
###################################################
## hhpdf("iinf-f2.pdf", height=4, width=7)
qqt.melt <- cbind(reshape2::melt(qqt.data[,2:5], value.name="quantile",
                                 variable.name="distribution"),
                  reshape2::melt(qqt.data[,9:6], value.name="density"))
levels(qqt.melt$distribution) <- c("t, 3 df", "t, 5 df", "t, 7 df", "normal")
head(qqt.melt)

xyplot(density ~ quantile, group=distribution, data=qqt.melt, type="l",
       par.settings=list(superpose.line=list(lty=c(5,4,2,1))),
       auto.key=list(space="right", border=TRUE, lines=TRUE, points=FALSE)) +
          layer(panel.abline(a=0, b=0, col="gray60"))
## hhdev.off()


###################################################
### code chunk number 24: iinf.tex:2259-2266
###################################################
## hhcapture("iinf-f3.Rout", '
rt5 <- rt(300, df=5)
rnn <- rnorm(300)

ks.test(rt5, function(x)pt(x, df=2))
ks.test(rnn, function(x)pt(x, df=2))
## ')


###################################################
### code chunk number 25: iinf.tex:2280-2309
###################################################
plotks1 <- function(x, pfunction, ..., xlim=range(x), ylim.dev=c(-.10, .10)) {
  sx <- sort(x)
  xset <- seq(xlim[1], xlim[2], length=101)
  xref <- xyplot(pfunction(xset) ~ xset, type="l", col="black", xlim=xlim, ylim=c(-.05, 1.05))
  empcdf <- (1:length(sx))/length(sx)
  resizePanels(h=c(1, 2.5),
               update(strip=FALSE, strip.left=TRUE, ylab=NULL, scales=list(rot=0),
                      xlim=xlim, ylim=list(ylim.dev, c(-.05, 1.05)),
               c(deviations=
                   xyplot(empcdf-pfunction(sx) ~ sx, type="h",
                          col=likertColor(2)[2], ...),
                 "cumulative distribution with deviations"=
                   xref + layer(panel.segments(sx, pfunction(sx), sx, empcdf,
                                               col=likertColor(2)[2]),
                           data=list(sx=sx, pfunction=pfunction, empcdf=empcdf)),
                 layout=c(1,2), x.same=TRUE, y.same=FALSE)))
}

KSt5 <- plotks1(rt5, function(x) pt(x, df=2), xlim=range(rt5, rnn),
                xlab="random selection from t with 5 df",
                main="compare 't with 5 df' to 't with 2 df'")
KSnn <- plotks1(rnn, function(x) pt(x, df=2), xlim=range(rt5, rnn),
                xlab="random selection from normal",
                main="compare 'normal' to 't with 2 df'")

## hhpdf("iinf-f3.pdf", height=6, width=7)
print(KSt5, position=c(0,0,.5,1), more=TRUE)
print(KSnn, position=c(.5,0,1,1), more=FALSE)
## hhdev.off()


###################################################
### code chunk number 26: iinf.tex:2338-2341
###################################################
## hhcapture("iinf-f32.Rout", '
ks.test(rt5, rnn)
## ')


###################################################
### code chunk number 27: iinf.tex:2354-2368
###################################################
plotks2 <- function(x, y, ..., xlim=range(x,y), ylim.dev=c(-.10, .10)) {
  sx <- sort(x)
  sy <- sort(y)
  empcdfx <- (1:length(sx))/length(sx)
  empcdfy <- (1:length(sy))/length(sy)

  xyplot(empcdfx ~ sx, xlim=xlim, ylim=c(-.05, 1.05), type="l", col=likertColor(2)[1], ...) +
     xyplot(empcdfy ~ sy, xlim=xlim, ylim=c(-.05, 1.05), type="l", col=likertColor(2)[2], ...)
}

## hhpdf("iinf-f32.pdf", height=3, width=7)
plotks2(rt5, rnn, xlab="random selections from t with 5 df (red) and normal (blue)",
        ylab="empirical CDF", main="Two-Sample Kolmogorov-Smirnov Test")
## hhdev.off()



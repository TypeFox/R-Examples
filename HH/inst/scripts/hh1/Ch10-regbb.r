  ## \cref{htwt.s}
  ## \cref{fabricwear.s}
  ## \cref{hotdog.s}
  ## \cref{hotdog-mmc2.s}
  ## \cref[splus.library]{ancova.s}.
  ## \cref{hotdog-ancova.s}
  ## -rwx------+   1 rmh None    53436 2004-05-25 01:03 regbb/regbb.tex

library(HH)

## htwt.s
data(htwt)
htwt <- htwt ## local copy

## detect two missing values
length(htwt$ht)
for (i in tapply(htwt$ht, htwt$sex, c))
  if.R(s=
       stem(i, scale=-1, nl=2)
       ,r=
       stem(i)
       )

## assign values to the missing observations
match(NA, htwt$ht, 0)
htwt[4,]
htwt[4,"ht"] <- round(1.65 * 39.37)  ## this student answered in meters
match(NA, htwt$sex, 0)
## (I had the class list and was able to make this corrections)
htwt[27,]
htwt[27,"sex"] <- "m"

## redo the stem and leaf with the completed data
for (i in tapply(htwt$ht, htwt$sex, c))
  if.R(s=
       stem(i, scale=-1, nl=2)
       ,r=
       stem(i)
       )

## scatterplot matrix of completed data
splom(~ htwt[,c("lbs","months","sex","ht")])
## same splom with control of font sizes
if.R(r=
     splom(~ htwt[,c("lbs","months","sex","ht")],
           axis.text.cex=.7)
     ,s=
     splom(~ htwt[,c("lbs","months","sex","ht")],
           superpanel=panel.pairs.hh,
           subpanel.scales=list(cex=.7))
     )
## export.eps(hh("regbb/figure/htwt.splom.eps"))

tpg <- trellis.par.get("superpose.symbol")
tpg <- lapply(tpg, function(x) x[1:2])
print(position=c(.2,.1,.8,.9),
xyplot(lbs ~ ht, data=htwt, panel=panel.superpose, groups=sex,
       aspect=1,
       key=list(
         space="right",
         border=TRUE,
         text=list(labels=levels(htwt$sex)),
         points=tpg))
)
## export.eps(hh("regbb/figure/htwt.xy.eps"))


## one-way analysis of variance
htwt.aov <- aov(ht ~ sex, data=htwt)
summary(htwt.aov)
model.tables(htwt.aov, type="means")


## dummy variable
htwt$female <- as.numeric(htwt$sex == "f")
htwt.lm <- lm(ht ~ female, data=htwt)
summary(htwt.lm, corr=FALSE)
anova(htwt.lm)

## dummy variable
htwt$treat <- (htwt$sex == "f") - (htwt$sex == "m")
htwtb.lm <- lm(ht ~ treat, data=htwt)
summary(htwtb.lm, corr=FALSE)
anova(htwtb.lm)


## fabricwear.s
## \cite{Ott:1993} datafile Ch11_9.DAT

data(fabricwear)

if.R(s=
     print(position=c(0,.4, 1,1),
           t(bwplot(speed ~ wear, data=fabricwear,
                    scales=list(cex=1.4),
                    xlab=list(cex=1.4), ylab=list(cex=1.4)))
           )
     ,r=
     bwplot(wear ~ speed, data=fabricwear,
            scales=list(cex=1.4), xlab=list(cex=1.4), ylab=list(cex=1.4))
     )
## export.eps(hh("regbb/figure/fabricwear.eps"))

fabricwear.aov <- aov(wear ~ speed, data=fabricwear, x=TRUE)
summary(fabricwear.aov)
model.tables(fabricwear.aov, "mean")

tmp.c <- contrasts(fabricwear$speed)
dimnames(tmp.c)[[1]] <- levels(fabricwear$speed)
tmp.c
zapsmall(crossprod(tmp.c))
tmp.min <- apply(abs(tmp.c), 2, min)
sweep(tmp.c, 2, tmp.min, "/")

tmp.tr <- data.frame(polynomial=as.vector(tmp.c),
                     speed=rep(as.numeric(dimnames(tmp.c)[[1]]),5),
                     power=rep(ordered(dimnames(tmp.c)[[2]],
                       levels=dimnames(tmp.c)[[2]]), c(6,6,6,6,6)))
print(position=c(.25,0,.75,1),
xyplot(polynomial ~ speed | power, data= tmp.tr, type="b", layout=c(1,5),
       between=list(y=1),
       ylab="normalized orthogonal polynomials",
       panel=function(...) {
         panel.xyplot(...)
         panel.abline(h=0, lty=2)
       },
       strip=function(...) strip.default(..., strip.names = c(TRUE, TRUE)),
       scales=list(cex=1, x=list(alternating=1), y=list(alternating=1)),
       par.strip.text=list(cex=1.4))
)
## export.eps(hh("regbb/figure/orthpoly.eps"))

summary(fabricwear.aov,
        split=list(speed=list(speed.L=1, speed.Q=2, speed.C=3, rest=4:5)))

dimnames(fabricwear.aov$x)[[2]]
summary.lm(fabricwear.aov, corr=FALSE)


## hotdog.s
## http://lib.stat.cmu.edu/DASL/Stories/Hotdogs.html

## Datafile Name: Hot dogs
## Datafile Subjects: Food
## Story Names: Hot dogs

## Reference: Moore, David S., and George P. McCabe (1989). Introduction
## to the Practice of Statistics. Original source: Consumer Reports, June
## 1986, pp. 366-367.

## Description: Results of a laboratory analysis of calories and
## sodium content of major hot dog brands. Researchers for Consumer
## Reports analyzed three types of hot dog: beef, poultry, and meat
## (mostly pork and beef, but up to 15% poultry meat).

## Number of cases: 54

## Variable Names:
##     Type: Type of hotdog (beef, meat, or poultry)
## Calories: Calories per hot dog
##   Sodium: Milligrams of sodium per hot dog

data(hotdog)

## oneway analysis of variance
if.R(r=bwplot(Sodium ~ Type, data=hotdog),
     s=t(bwplot(Type ~ Sodium, data=hotdog)))
## export.eps(hh("regbb/figure/hotdog1.eps"))

## horizontal lines: zero slope and separate intercepts
  hotdog.aov <- ancova(Sodium ~ Type, data=hotdog, x=Calories,
                       par.strip.text=list(cex=1.2), ylim=c(140,700))
  print(position=c(0,0, 1,.6),
        attr(hotdog.aov,"trellis"))
## export.eps(hh("regbb/figure/hotdog.f0.eps"))

summary(hotdog.aov)
model.tables(hotdog.aov, type="means")

## multiple comparisons of ANOVA
if.R(r=
     {
       hotdog.mca <- glht(hotdog.aov, linfct = mcp(Type="Tukey"))
       print(hotdog.mca)
       old.omd <- par(omd=c(.1,1, 0,1))
       plot(hotdog.mca)
       par(old.omd)
     }, s={
       hotdog.mca <- multicomp(hotdog.aov, focus="Type")
       print(hotdog.mca)
       hotdog.mca <- multicomp.reverse(hotdog.mca)  ## positive differences
       print(hotdog.mca)
       plot(hotdog.mca)
     })
## export.eps(hh("regbb/figure/hotdog2.eps"))



## regression
## same line: common intercept and common slope
     hC.aov <- ancova(Sodium ~ Calories, groups=Type, data=hotdog,
                      par.strip.text=list(cex=1.2), ylim=c(140,700))
     print(position=c(0,0, 1,.6),
           attr(hC.aov,"trellis"))
## export.eps(hh("regbb/figure/hotdog.f3.eps"))
summary.aov(hC.aov)

## analysis of covariance
## analysis with a concomitant explanatory variable
## parallel lines: separate intercepts and common slope
hCT.aov <- ancova(Sodium ~ Calories + Type, data=hotdog,
                  par.strip.text=list(cex=1.2), ylim=c(140,700))
print(position=c(0,0, 1,.6),
      attr(hCT.aov,"trellis"))
## export.eps(hh("regbb/figure/hotdog.f1.eps"))
summary(hCT.aov)

## multiple comparisons of ANCOVA
if.R(r=
     {
       hotdog.mca <- glht(hCT.aov, linfct=mcp(Type="Tukey"))
       print(hotdog.mca)
       old.omd <- par(omd=c(.1,1, 0,1))
       plot(hotdog.mca)
       par(old.omd)
     }, s={
       hCT.mca <- multicomp(hCT.aov, focus="Type")
       print(hCT.mca)
       hCT.mca <- multicomp.reverse(hCT.mca)  ## positive differences
       print(hCT.mca)
       plot(hCT.mca)
     })
## export.eps(hh("regbb/figure/hotdog3.eps"))

## same graph as above, anova table shows a different sequence
## hTC.aov <- ancova(Sodium ~ Type + Calories, data=hotdog,
##                   par.strip.text=list(cex=1.2))
## summary(hTC.aov)


## analysis of response Sodium adjusted for the covariable Calories
## parallel lines from hCT.aov become horizontal with the same
## vertical differences between ?intercepts.
predict.lm(hCT.aov, type="terms")
hotdog <- cbind(hotdog, Sodium.Calories=hotdog$Sodium -
                predict.lm(hCT.aov, type="terms", terms="Calories")[,1])
hSCT.aov <- ancova(Sodium.Calories ~ Type, x=Calories, data=hotdog,
                  par.strip.text=list(cex=1.2), ylim=c(140,700))
print(position=c(0,0, 1,.6),
      attr(hSCT.aov,"trellis"))
## export.eps(hh("regbb/figure/hotdog.f4.eps"))
summary(hSCT.aov)

model.tables(hSCT.aov, type="means")


## interaction: separate intercepts and slopes
hCTi.aov <- ancova(Sodium ~ Calories * Type, data=hotdog,
                   par.strip.text=list(cex=1.2), ylim=c(140,700))
print(position=c(0,0, 1,.6),
      attr(hCTi.aov,"trellis"))
## export.eps(hh("regbb/figure/hotdog.f2.eps"))
summary(hCTi.aov)

## hTCi.aov <- ancova(Sodium ~ Type * Calories, data=hotdog,
##                   par.strip.text=list(cex=1.2), ylim=c(140,700))
## print(position=c(0,0, 1,.6),
##       attr(hCTi.aov,"trellis"))
## summary.aov(hTCi.aov)



coef(summary.lm(hCT.aov))
coef(summary.lm(hCTi.aov))

summary(hCT.aov)
summary(hCTi.aov)

hCTx.lm <- lm(Sodium ~ Calories + Type, data=hotdog, x=TRUE)

names(hCT.aov)
names(hCTx.lm)
hCTx.lm$x
contrasts(hotdog$Type)
contr.helmert(3)

Type1.h <- hCTx.lm$x[,3]
Type2.h <- hCTx.lm$x[,4]

hCT.h.lm <-
  lm(Sodium ~ Calories + Type1.h + Type2.h, data=hotdog)
summary(hCT.h.lm, corr=FALSE)
summary.aov(hCT.h.lm)
coef(summary(hCT.h.lm))
coef(hCT.h.lm)

summary(hCT.aov)
coef(summary.lm(hCT.aov))

hCTix.lm <- lm(Sodium ~ Calories + Type + Calories:Type,
                    data=hotdog, x=TRUE)
hCTix.lm$x



##Other sets of contrasts

contr.helmert(3)
contr.sum(3)
contr.poly(3)
contr.treatment(3)

contrasts(hotdog$Type)
contrasts(hotdog$Type) <- contr.sum(3)
contrasts(hotdog$Type)

hCTsx.lm <- lm(Sodium ~ Calories + Type, data=hotdog, x=TRUE)
coef(summary(hCTsx.lm))
coef(summary.lm(hCT.aov))
summary.aov(hCTsx.lm)
hCTsx.lm$x

## Conclusion. The regression coefficients with dummy variables depend
## on the choice of dummy variables.  The anova table does not depend
## on the choice of dummy variables.  The graphs do not depend on the
## choice of dummy variables.


## I decided to look at  Sodium ~ Calories


## We could also look at Calories ~ Sodium
## parallel lines
hST.aov <- ancova(Calories ~ Sodium + Type, data=hotdog,
                  par.strip.text=list(cex=1.2))
print(position=c(0,0, 1,.6),
      attr(hST.aov,"trellis"))
## export.eps(hh("regbb/figure/hotdog.f5.eps"))
summary(hST.aov)

if.R(r={
## hotdog-mmc2.r
## follows hotdog.s

old.par <- par(mar=c(5,4,4,5.5)+.1)

hCT.mmc <- mmc(hCT.aov)
print(hCT.mmc)
plot(hCT.mmc)
## export.eps(hh("regbb/figure/hotdog-mmc-mca.eps"))

##                                  B  M  P
hCT.lmat <- cbind("BeefMeat-Poul"=c(1, 1,-2),
                  "Beef-Meat"    =c(1,-1, 0))
dimnames(hCT.lmat)[[1]] <-  levels(hotdog$Type)
hCT.lmat

hCT.mmc <- mmc(hCT.aov, linfct=mcp(Type="Tukey"), focus.lmat=hCT.lmat)
print(hCT.mmc)
plot(hCT.mmc)
## export.eps(hh("regbb/figure/hotdog-mmc-lmat.eps"))

par(old.par)
}, s={
## hotdog-mmc2.s
## follows hotdog.s

## find out which rows of lmat we need
tmp <- multicomp(hCT.aov)$lmat
zapsmall(tmp)

old.par <- par(mar=c(5,4,4,5.5)+.1)

multicomp.mmc(hCT.aov, focus="Type")
## export.eps(hh("regbb/figure/hotdog-mmc-mca.eps"))

##                                 (I) C  B  M  P
hCT.lmat <- cbind("BeefMeat-Poul"=c(0, 0, 1, 1,-2),
                  "Beef-Meat"    =c(0, 0, 1,-1, 0))
dimnames(hCT.lmat)[[1]] <- dimnames(tmp)[[1]]
hCT.lmat

multicomp.mmc(hCT.aov, focus="Type", lmat=hCT.lmat)
## export.eps(hh("regbb/figure/hotdog-mmc-lmat.eps"))

par(old.par)
})


## splus.library/ancova.s
## functions in the HH package


## {
## This section collects the ancova calls from above into HH Table 10.14
## hotdog-ancova.s

## y ~ x                     ## constant line across all groups
ancova(Sodium ~ Calories,     data=hotdog, groups=Type)

## y ~ a                     ## different horizontal line in each group
ancova(Sodium ~            Type, data=hotdog, x=Calories)

## y ~ x + a  or  y ~ a + x  ## constant slope, different intercepts
ancova(Sodium ~ Calories + Type, data=hotdog)
ancova(Sodium ~ Type + Calories, data=hotdog)

## y ~ x * a  or  y ~ a * x  ## different slopes, and different intercepts
ancova(Sodium ~ Calories * Type, data=hotdog)
ancova(Sodium ~ Type * Calories, data=hotdog)
## }


{
## This section collects all ancova graphs from above into a single multi-graph plot.
## based on jsm2004.hotdog.r
## based on jsm2004.hotdog.s
## based on presentations/handbook/code/hotdog.csc.s

removeLegendAxes <-
  if.R(r=
       function(x) {
         x$ylab <- NULL
         x$xlab <- NULL
         x$legend <- NULL
         x$x.scales$alternating <- 0
         x$y.scales$alternating <- 0
         x$x.scales$tck <- c(0,0)
         x$y.scales$tck <- c(0,0)
         x$par.strip.text$cex <- .7
         x
       }, s=
       function(x) {
         x$ylab <- NULL
         x$xlab <- NULL
         x$key <- NULL
         x$scales$alternating <- 0
         x$par.strip.text$cex <- .8
         x
       }
       )

## 2 x 3, with empty spots
print(position=c(0.03, 1/3,  .53, 2/3), more=TRUE,
removeLegendAxes(attr(hC.aov,"trellis"))
)
print(position=c(.50, 0, 1.00, 1/3), more=TRUE,
removeLegendAxes(attr(hotdog.aov,"trellis"))
)
print(position=c(.50, 1/3, 1.00, 2/3), more=TRUE,
removeLegendAxes(attr(hCT.aov,"trellis"))
)
print(position=c(.50, 2/3, 1.00, 1), more=TRUE,
removeLegendAxes(attr(hCTi.aov,"trellis"))
)

if.R(r={
print(position=c(.17, 0, .42, .10), more=TRUE,
      xyplot(0 ~ .5, panel=function(...){}, xlab=expression("constant intercept" ~~ alpha),
             ylab="", scales=list(draw=FALSE),
             par.settings = list(axis.line = list(col = "transparent")))
)
print(position=c(5/8, 0, 7/8, .10), more=TRUE,
      xyplot(0 ~ .5, panel=function(...){}, xlab=expression("variable intercept" ~~ alpha),
             ylab="", scales=list(draw=FALSE),
             par.settings = list(axis.line = list(col = "transparent")))
)

print(position=c(0, .10, .10, .20), more=TRUE,
      xyplot(0 ~ .5, panel=function(...){}, ylab=expression("zero slope" ~~ beta==0),
             xlab="", scales=list(draw=FALSE),
             par.settings = list(axis.line = list(col = "transparent")))
)
print(position=c(0, .45, .10, .55), more=TRUE,
      xyplot(0 ~ .5, panel=function(...){}, ylab=expression("constant slope" ~~ beta),
             xlab="", scales=list(draw=FALSE),
             par.settings = list(axis.line = list(col = "transparent")))
)
print(position=c(0, .80, .10, .90), more=FALSE,
      xyplot(0 ~ .5, panel=function(...){}, ylab=expression("variable slope" ~~ beta),
             xlab="", scales=list(draw=FALSE),
             par.settings = list(axis.line = list(col = "transparent")))
)
}, s={
mtext(side=1, outer=TRUE, at=.30, text="constant intercept alpha", line=-1)
mtext(side=1, outer=TRUE, at=.65, text="variable intercept alpha", line=-1)

mtext(side=2, outer=TRUE, at=.45, text="zero slope beta=0",   line=-1, srt=90)
mtext(side=2, outer=TRUE, at=.62, text="constant slope beta", line=-1, srt=90)
mtext(side=2, outer=TRUE, at=.79, text="variable slope beta", line=-1, srt=90)
par(mfrow=c(1,1))
}
## export.eps(hh("regbb/figure/hotdog.jsm.eps"))
## export.eps(hh("presentations/handbook/figure/hotdog.csc.eps"))
)
}

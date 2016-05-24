## still need to verify all the position=pp examples

## splus.library/logit.s
## logit-plot.s
## spaceshuttle.s
## spac.glma.s
## spac.glmb.s
## spac.glmc.s
## budworm.s
## lymph.s
## lymph2.s
## work/logit-da.s
## work/logit-ga.s
## icu.s
##  -rwx------+   1 rmh None    67555 2004-05-25 01:03 logi/logi.tex


## splus.library/logit.s
## included in HH package



## logit-plot.s
p <- seq(0,1,length=51)
y <- logit(p)

print(position=c(0,.1,.4,.9), more=TRUE,
      xyplot(y ~ p, type="l", main=list("y = logit(p)",cex=1.5),
             xlab=list(cex=1.3), ylab=list(cex=1.3),
             scales=list(cex=1.3, x=list(at=c(0,.25,.5,.75,1),
                                    labels=c("0.0","","0.5","","1.0"))))
      )

print(position=c(.4,.2,1,.75), more=FALSE,
      xyplot(p ~ y, type="l", main=list("p = antilogit(y)",cex=1.5),
             xlab=list(cex=1.3), ylab=list(cex=1.3),
             scales=list(cex=1.3, y=list(at=c(0,.25,.5,.75,1),
                                    labels=c("0.0","","0.5","","1.0"))))
      )

## export.eps(hh("logi/figure/logit-plot.eps"))



## spaceshuttle.s
data(spacshu)

## We used
##     trellis.device.hh.bw(orientation="portrait")
## to create the graphs marked "portrait".
## We used
##     trellis.device.hh.bw()
## to create the graphs marked "landscape".

Jitter.factor <- if.R(s=1, r=.08)

Fig.a <- 
print(split=c(1,3, 1,3), more=TRUE,
      xyplot(jitter(damage, factor=Jitter.factor*8) ~ tempF,
             ylab="damage",
             data=spacshu,
             scales=list(y=list(at=c(0,1))),
             cex=.5, pch=16,
             main="a. observed",
             panel=function(...) {
               panel.xyplot(...)
               panel.abline(h=0:1, lty=2)
             }))

temp.range.left <- seq(50,80,5)
temp.range.right <- seq(55,85,5)
temp.range.both <- seq(50,85,5)
damage.freq <- table(cut(spacshu$tempF, temp.range.both), spacshu$damage)
damage.prop <- damage.freq[,"1"] / (damage.freq[,"0"]+damage.freq[,"1"])

print(split=c(1,2, 1,3), more=FALSE,
      xyplot(jitter(damage, factor=Jitter.factor*8) ~ tempF,
             ylab="damage",
             data=spacshu,
             scales=list(y=list(at=c(0,1))),
             cex=.5, pch=16,
             main="b. observed and sectioned proportions",
             panel=function(...) {
               panel.xyplot(...)
               if.R(r=segments <- panel.segments,
                    s={})
               segments(temp.range.left, damage.prop,
                        temp.range.right, damage.prop)
               panel.abline(h=0:1, lty=2)
             }))
## export.eps(hh("logi/figure/spaceshuttle.ab.eps")) ## "portrait"


print(split=c(1,3, 1,3), more=TRUE,
      xyplot(jitter(damage, factor=Jitter.factor*8) ~ tempF,
             ylab="damage",
             data=spacshu,
             scales=list(y=list(at=c(0,1))),
             cex=.5, pch=16,
             main="c. observed and sectioned proportions\nappropriate temperature scale",
             xlim=c(30,85), ylim=c(-.1,1.5),
             panel=function(...) {
               panel.xyplot(...)
               if.R(r=segments <- panel.segments,
                    s={})
               segments(temp.range.left, damage.prop,
                        temp.range.right, damage.prop)
               panel.abline(h=0:1, lty=2)
             }))



spacshu.bin.glm <- glm(damage ~ tempF, data=spacshu, family=binomial)
spacshu.bin.glm
anova(spacshu.bin.glm, test="Chi")
coef(summary(spacshu.bin.glm))

## prediction on response scale, in this case (0,1).
## leading to Figure logi/figure/spaceshuttle.cde.eps Panel d
spacshu.pred <- 
  predict(spacshu.bin.glm, data.frame(tempF=30:85), type="response",
          se.fit=TRUE, ci.fit=TRUE, pi.fit=TRUE)

print(split=c(1,2, 1,3), more=TRUE,
xyplot(jitter(damage, factor=Jitter.factor*4) ~ tempF, data=spacshu,
       ylab="proportion damaged",
       main="d. glm logit fit with pi, estimating p(damage in one ring)",
       xlim=c(30,85), ylim=c(-.1,1.5),
       scales=list(y=list(at=c(0,1))),
       cex=.5, pch=16,
       panel=function(...) {
	 panel.xyplot(...)
         if.R(r=lines <- panel.lines,
              s={})
	 lines(x=30:85, y=spacshu.pred$fit,
	       lty=trellis.par.get("superpose.line")$lty[1],
	       col=trellis.par.get("superpose.line")$col[1])
	 lines(x=30:85, y=spacshu.pred$pi.fit[,"lower"],
	       lty=trellis.par.get("superpose.line")$lty[4],
	       col=trellis.par.get("superpose.line")$col[4])
	 lines(x=30:85, y=spacshu.pred$pi.fit[,"upper"],
	       lty=trellis.par.get("superpose.line")$lty[4],
	       col=trellis.par.get("superpose.line")$col[4])
         panel.abline(h=0:1, lty=2)
       }
)
)

spacshu6 <- data.matrix(spacshu)
dimnames(spacshu6) <- NULL
dim(spacshu6) <- c(6,23,2)
spacshu6 <- data.frame(damage=apply(spacshu6[,,2],2,sum), tempF=spacshu6[1,,1])


print(split=c(1,1, 1,3), more=FALSE,
xyplot(jitter(damage, factor=Jitter.factor*3) ~ tempF, data=spacshu6,
       ylab="number of damaged rings",
       main="e. glm logit fit with pi, estimating number of damaged rings",
       xlim=c(30,85), ylim=c(-1,8),
       scales=list(y=list(at=c(0,6))),
       cex=.5, pch=4,
       panel=function(...) {
	 panel.xyplot(...)
         if.R(r=lines <- panel.lines,
              s={})
	 lines(x=30:85, y=6*spacshu.pred$fit,
	       lty=trellis.par.get("superpose.line")$lty[1],
	       col=trellis.par.get("superpose.line")$col[1])
	 lines(x=30:85, y=6*spacshu.pred$pi.fit[,"lower"],
	       lty=trellis.par.get("superpose.line")$lty[4],
	       col=trellis.par.get("superpose.line")$col[4])
	 lines(x=30:85, y=6*spacshu.pred$pi.fit[,"upper"],
	       lty=trellis.par.get("superpose.line")$lty[4],
	       col=trellis.par.get("superpose.line")$col[4])
         panel.abline(h=c(0,6), lty=2)
       }
)
)
## export.eps(hh("logi/figure/spaceshuttle.cde.eps")) ## "portrait"





## all on one graph

tpgsl <- trellis.par.get("superpose.line")
if.R(r=tpgsl$alpha <- rep(tpgsl$alpha, 4),
     s={})
tpgsl <- lapply(tpgsl, function(x) x[1:4])
tpgsl$col <- c(0,0,tpgsl$col[c(1,4)])
tpgsl$lty <- c(0,0,tpgsl$lty[c(1,4)])
tpgsl$lwd[] <- 2

tmp <-
xyplot(jitter(damage, factor=Jitter.factor*6) ~ tempF, data=spacshu,
       ylab="proportion damaged",
       main="glm logit fit with pi",
       xlim=c(30,85), ylim=c(-.1,1.3),
       scales=list(y=list(at=c(0,1))),
       cex=.6,
       key=list(
         ## y=-.15, ## portrait
         y=-.2, ## landscape
         text=list(c("individual ring","average","predicted","PI")),
         points=list(pch=c(16,4,32,32)),
         lines=tpgsl,
         space="bottom",
         border=TRUE),
       panel=function(x, y, ...) {
	 panel.xyplot(x, y, ..., pch=16)
         panel.xyplot(y=spacshu6$damage/6, x=spacshu6$tempF, pch=4)
         if.R(
              r={
                cpl <- current.panel.limits()
                pushViewport(viewport(xscale = cpl$xlim,
                                      yscale = cpl$ylim,
                                      clip = "off"))
                ## put anything you want unclipped inside this:
                panel.axis(side="right", at=(0:6)/6, label=0:6, outside=TRUE)
                panel.text(x=diff(cpl$xlim)*.05+cpl$xlim[2],
                           y=mean(cpl$ylim),
                           "number damaged", srt=90)
                ## end of unclipped part
                upViewport()
              },s={
                axis(4, at=(0:6)/6, label=0:6)
                mtext("number damaged", side=4, line=2)
              })
         if.R(r=lines <- panel.lines,
              s={})
	 lines(x=30:85, y=spacshu.pred$fit,
	       lty=trellis.par.get("superpose.line")$lty[1],
	       col=trellis.par.get("superpose.line")$col[1])
	 lines(x=30:85, y=spacshu.pred$pi.fit[,"lower"],
	       lty=trellis.par.get("superpose.line")$lty[4],
	       col=trellis.par.get("superpose.line")$col[4])
	 lines(x=30:85, y=spacshu.pred$pi.fit[,"upper"],
	       lty=trellis.par.get("superpose.line")$lty[4],
	       col=trellis.par.get("superpose.line")$col[4])
         panel.abline(h=0:1, lty=2)
       }
)
print(position=c(0,0,1,.7), tmp)

## export.eps(hh("logi/figure/spaceshuttle.all.eps")) ## "portrait"



print(split=c(1,3, 1,3), more=FALSE,
      xyplot(jitter(damage, factor=Jitter.factor*8) ~ tempF,
             ylab="damage",
             xlim=if.R(s=Fig.a$xlim, r=Fig.a$x.limits),
             ylim=if.R(s=Fig.a$ylim, r=Fig.a$y.limits),
             data=spacshu[spacshu[,"damage"]==1,],
             scales=list(y=list(at=c(0,1))),
             cex=.5, pch=16,
             main="f. observed damaged O-rings",
             panel=function(...) {
               panel.xyplot(...)
               panel.abline(h=0:1, lty=2)
             }))
## export.eps(hh("logi/figure/spaceshuttle.damaged.eps")) ## "portrait"


## transformations
tmp <- cbind.data.frame(tempF=30:85,
                        p=spacshu.pred$fit,
                        odds=spacshu.pred$fit/(1-spacshu.pred$fit),
                        logit.p=logit(spacshu.pred$fit),
                        which="fit", stringsAsFactors=FALSE)
p.obs <- spacshu6[,1]/6
tmp2 <- cbind.data.frame(spacshu6[,2,drop=FALSE],
                         p=p.obs,
                         odds=p.obs/(1-p.obs),
                         logit.p=logit(p.obs),
                         which="data", stringsAsFactors=FALSE)
tmp2$logit.p[tmp2$p==0] <- min(tmp2$logit.p[tmp2$p!=0])-2.5  ## pull infinity in
tmp.both <- rbind(tmp2,tmp)
splom(~tmp.both[,1:4], type=c("p","l"), groups=tmp.both$which,
      distribute.type=TRUE,  ## R only, ignored by S-Plus
      panel=panel.superpose, lty=1, pch=16,
      superpanel=if.R(s=panel.pairs.hh, r=panel.pairs),   ## HH
      panel.cex=1.4, subpanel.scales=list(cex=1),
      cex=.9, pscales=3)
## export.eps(hh("logi/figure/spaceshuttle.logit-splom.eps")) ## "portrait"

print(position=c(.2,0, .8,1),
xysplom(p + odds + logit.p ~ tempF, data=tmp.both,
        type=c("p","l"), groups=rep(tmp.both$which,3),
        distribute.type=TRUE,  ## R only, ignored by S-Plus
        panel=panel.superpose, lty=1,
        adj.ylim=c(0,.2),  ## force p ~ tempF to have ylim=c(0,1)
        par.strip.text=list(cex=1.4),
        scales=list(cex=1, y=list(relation="free"), alternating=FALSE),
        between=list(y=1.5), xlab="", ylab="")
)
## export.eps(hh("logi/figure/spaceshuttle.logit-xysplom.eps")) ## landscape



## prediction on link scale, in this case $(-\infty,\infty)
## leading to Figure logi/figure/spaceshuttle.pred.logit.f.eps Panel g
spacshu.pred.link <-
  predict(spacshu.bin.glm, data.frame(tempF=30:85), type="link",
          se.fit=TRUE, ci.fit=TRUE, pi.fit=TRUE)

names(spacshu.pred)

cbind(tempF=30:85,
      round(as.data.frame(spacshu.pred.link[-(3:4)]), digits=2))[c(1:3,54:56),]
cbind(tempF=30:85,
      round(as.data.frame(spacshu.pred[-(3:4)]), digits=2))[c(1:3,54:56),]

## approximate logit(1) with 4.1, and logit(0) with -4.1
jld <- jitter(ifelse(spacshu$damage==1, 4.1, -4.1),8) ## pull infinity in

print(split=c(1,1, 1,3), more=FALSE,
xyplot(jld ~ tempF, data=spacshu,
       ylab="logit(proportion) damaged",
       main="g. glm logit fit with logit(pi), estimating logit(p(damage in one ring))",
       xlim=c(30,85), ylim=c(-8,4),
       cex=.5, pch=16,
       panel=function(...) {
	 panel.xyplot(...)
         if.R(r=lines <- panel.lines,
              s={})
	 lines(x=30:85, y=spacshu.pred.link$fit,
	       lty=trellis.par.get("superpose.line")$lty[1],
	       col=trellis.par.get("superpose.line")$col[1])
	 lines(x=30:85, y=spacshu.pred.link$pi.fit[,"lower"],
	       lty=trellis.par.get("superpose.line")$lty[4],
	       col=trellis.par.get("superpose.line")$col[4])
	 lines(x=30:85, y=spacshu.pred.link$pi.fit[,"upper"],
	       lty=trellis.par.get("superpose.line")$lty[4],
	       col=trellis.par.get("superpose.line")$col[4])
       }
       
)
)
## export.eps(hh("logi/figure/spaceshuttle.pred.logit.f.eps")) ## "portrait"




## spac.glma.s
spacshu.bin.glm <- glm(damage ~ tempF, data=spacshu,
                       family=binomial)




## spac.glmb.s
spacshu.bin.glm <- glm(damage ~ tempF, data=spacshu,
                       family=binomial(link=logit))




## spac.glmc.s
anova(spacshu.bin.glm, test="Chi")




## budworm.s
##                         Logistic Regression
##
##                Richard M. Heiberger and Burt Holland
##                          Temple University

### A. logistic regression
##
## Based on Venables and Ripley, MASS 4th Edition, Chapter 7
##
## The budworm data and the VR analysis is in file
##          MASS/scripts/ch07.ssc
## This HH file is based on the VR file.


library(MASS)  ## we use the VR function dose.p

options(contrasts = c("contr.treatment", "contr.poly"))

data(budworm)
SF <- cbind(numdead=budworm$numdead,
            numalive = 20 - budworm$numdead)


## model with interaction term for sex and logdose, from VR
budworm.lg <-
  glm(SF ~ sex*ldose,
      data=budworm,
      family = binomial)
summary(budworm.lg, cor = FALSE)
anova(budworm.lg, test="Chisq")

## VR plot with interaction term using regular graphics
plot(c(1,32), c(0,1), type = "n",
     xlab = "dose", ylab = "prob",
     log = "x")
text(2^budworm$ldose, budworm$numdead/20,
     labels = as.character(budworm$sex))
ld <- seq(0, 5, 0.1)
lines(2^ld,
      predict(budworm.lg,
              data.frame(ldose = ld,
                         sex = factor(rep("M", length(ld)),
                           levels = levels(budworm$sex))),
              type = "response"), col = 3)
lines(2^ld,
      predict(budworm.lg,
              data.frame(ldose = ld,
                         sex = factor(rep("F", length(ld)),
                           levels = levels(budworm$sex))),
              type = "response"), lty = 2, col = 2)


## model with no interaction term
budworm.lg0 <- glm(SF ~ sex + ldose - 1,
                   data=budworm,
                   family = binomial)
summary(budworm.lg0, cor = FALSE)$coefficients
anova(budworm.lg0, test="Chisq")

## LD25 LD50 LD75
xp.M <- dose.p(budworm.lg0, cf = c(2,3), p = 1:3/4)
xp.M
xp.F <- dose.p(budworm.lg0, cf = c(1,3), p = 1:3/4)
xp.F


## HH plot using trellis graphics

## get same colors for both observed and fit
tpgss.original <- trellis.par.get("superpose.symbol")
tpgsl.original <- trellis.par.get("superpose.line")

tpgss <- lapply(tpgss.original,
                function(x) x[c(1:2,1:2)])
tpgss$pch[tpgss$pch==1] <- 16
tpgsl <- lapply(tpgsl.original,
                function(x) x[c(1:2,1:2)])
tpgsl$lty[tpgsl$lty==2] <- 3
if.R(r=
     tpgsl$alpha <- rep(tpgsl.original$alpha, 4)
     ,s={})

## the plot will use the colors that have been saved back to the device
trellis.par.set("superpose.symbol", tpgss)
trellis.par.set("superpose.line",   tpgsl)  


logit.p.hat <- predict.glm(budworm.lg0, type="link")
p.hat <- predict.glm(budworm.lg0, type="response")
odds.hat <- p.hat/(1-p.hat)

lhat <- cbind(budworm, p=budworm$numdead/20,
              p.hat=p.hat, odds.hat=odds.hat, logit.p.hat=logit.p.hat)
names(lhat)
lhat.srt <- lhat[order(lhat$sex, lhat$ldose),]
splom(~lhat.srt[-(2:3)] | lhat.srt$sex , type="b",
      par.strip.text=list(cex=1.2), main="observed p, fitted others")
splom(~lhat.srt[-(2:3)], groups=lhat.srt$sex , type="b",
      panel=panel.superpose, main="observed p, fitted others")

## finer interpolation for fit
lhat2 <- data.frame(ldose=rep(seq(-1, 6, 0.1),2))
lhat2$sex <- factor(rep(c("F","M"), c(length(lhat2$ld)/2, length(lhat2$ld)/2)),
                        levels = levels(budworm$sex))
lhat2$p <- NA
lhat2$p.hat=predict(budworm.lg0, lhat2, type = "response")
lhat2$odds.hat <- lhat2$p.hat/(1-lhat2$p.hat)
lhat2$logit.p.hat=predict(budworm.lg0, lhat2, type = "link")
lhat2[c(1,2,70,71,72,73,141,142),]
splom(~lhat2[-(2:3)] | lhat2$sex , type="l",
      par.strip.text=list(cex=1.2), main="fit only, no observed points")
splom(~lhat2[-(2:3)], groups=lhat2$sex, type="l",
      panel=panel.superpose, main="fit only, no observed points")


## put the observed points in all panels with fits
lhat.p <- lhat
lhat.p$p.hat <- lhat.p$p
lhat.p$odds.hat <- lhat.p$p.hat/(1-lhat.p$p.hat)
lhat.p$logit.p.hat <- log(lhat.p$odds.hat)
tmp <- lhat.p[-3]
tmp[abs(tmp)==Inf] <- NA
lhat.p[-3] <- tmp

tmp2 <- rbind(cbind(lhat.p[-2], observed=TRUE), cbind(lhat2, observed=FALSE))
tmp2$so <- interaction(tmp2$sex, tmp2$observed)
names(tmp2)

## observed and fit in same panels
tmp3 <- tmp2[-c(2,3,7,8)]
names(tmp3) <- c("log dose","prob","odds","logit")
splom(~tmp3, groups=tmp2$so,
      type=c("l","l","p","p"),
      panel=panel.superpose, main="too many colors, no key")



## subset the colors for the key
## reverse the order to match the appearance in the plot
tpgss <- lapply(tpgss,
                function(x) x[c(2:1)])
tpgsl <- lapply(tpgsl,
                function(x) x[c(2:1)])

## observed and fit with same colors in same panels with key
## bigger fonts
splom(~tmp3, groups=tmp2$so,
      type=c("l","l","p","p"),
      distribute.type=TRUE,  ## R only, ignored by S-Plus
      panel=panel.superpose,
      superpanel=if.R(s=panel.pairs.hh, r=panel.pairs),   ## HH
      subpanel.scales=list(cex=.7), ## HH
      panel.cex=1.4,                ## HH
      key=list(space="right", border=1,
        text=list(c("Male","Female")),
        points=tpgss,
        lines=tpgsl,
        title="sex   observed   fit",
        cex.title=1.2))
## export.eps(hh("logi/figure/budworm-7.eps"))
## export.eps(hh("logi/figure/budworm-7-color.eps"))

## the single panel that we actually publish
xyplot(p.hat ~ ldose, groups=so, data=tmp2,
       type=c("l","l","p","p"), cex=1,
       distribute.type=TRUE,  ## R only, ignored by S-Plus
       xlab=list(cex=1.4), ylab=list("prob",cex=1.4),
       scales=list(cex=1.2),
       panel=panel.superpose,
       main=list("observed points and fitted logistic regression",cex=1.4),
       key=list(space="right", border=1,
         text=list(c("Male","Female")),
         points=tpgss,
         lines=tpgsl,
         title="sex   observed   fit",
         cex.title=1.2))
## export.eps(hh("logi/figure/budworm-8.eps"))
## export.eps(hh("logi/figure/budworm-8-color.eps"))


## restore standard colors
trellis.par.set("superpose.symbol", tpgss.original)
trellis.par.set("superpose.line",   tpgsl.original)  

if.R(s=detach("MASS"),
     r=detach("package:MASS"))

## Design Issues
## 
## a. simple interface (relatively simple)
##
## b. put variables in sensible order
##
## c. show both the parallel lines of the logit scale and the
##    interpretable probability scale
##
## d. different scales for residuals, with an argument to select
##
## e. sort order matters when dots are connected, not necessarily otherwise
##
## f. statistics need standard errors. VR in dose.p provide them.
##
## g. lots of graphs, until we get it right, then use just the panel
##    that makes the point.
##
## h. choice of similar colors for similar objects
##
## i. construction of an interaction variable for use with superpose
##
## j. HH functions give finer control over parts of the plot




## lymph.s
data(lymph)
     
## larger character size in plots
tpgss <- Rows(trellis.par.get("superpose.symbol"), seq(levels(lymph$nodes)))
tpgss$cex[] <-1.4
statistic.name <- if.R(s="t value", r="z value")


## glm                
lymph1a.glm <- glm(nodes ~ X.ray + acid.ph, data=lymph,
                   family=binomial)

anova(lymph1a.glm, test="Chisq")

summary(lymph1a.glm)$coef




## 1                  
## a.eps.gz
print(position=c(0,.1, 1,1), more=FALSE,  ## top 90%
xyplot(age ~ acid.ph | stage * grade * X.ray, data=lymph,
       panel=panel.superpose, group=lymph$nodes,
       layout=c(4,2),
       main=list("age ~ acid.ph | stage * grade * X.ray", cex=1.6),
## above is necessary, below makes it prettier
       scales=list(cex=1, alternating=FALSE),
       xlab=list(cex=1.4), ylab=list(cex=1.4),
       between=list(y=1),
       par.strip.text=list(cex=1.4), cex=1.4,
       strip=function(...) strip.default(..., style=1, strip.names=c(TRUE,TRUE)),
       key=list(space="right",
         text=list(levels(lymph$nodes), adj=1),
         points=tpgss,
         border=1,
         title="nodes", cex.title=1.2,
         cex=1))
)
## export.eps(hh("logi/figure/a.eps"))



## 2
## b.eps.gz
print(position=c(0,.1, 1,1), more=FALSE,  ## top 90%
xyplot(nodes.j ~ acid.ph | stage * grade * X.ray, data=lymph,
       layout=c(4,2),
       main=list("nodes ~ acid.ph | stage * grade * X.ray", cex=1.6),
## above is necessary, below makes it prettier
       scales=list(cex=1, alternating=FALSE),
       xlab=list(cex=1.4), ylab=list(cex=1.4),
       between=list(y=1),
       par.strip.text=list(cex=1.4), cex=1.1,
       strip=function(...) strip.default(..., style=1, strip.names=c(TRUE,TRUE)))
)
## export.eps(hh("logi/figure/b.eps"))


## 3
## c.eps.gz
## ignore grade and stage, simplification for exposition
print(position=c(0,.5, 1,1), more=FALSE,  ## top 50%
xyplot(nodes.j ~ acid.ph | X.ray, data=lymph,
       main=list("nodes ~ acid.ph | X.ray", cex=1.6),
## above is necessary, below makes it prettier
       scales=list(cex=1, alternating=FALSE),
       xlab=list(cex=1.4), ylab=list(cex=1.4),
       par.strip.text=list(cex=1.4), cex=1.1,
       strip=function(...) strip.default(..., style=1, strip.names=c(TRUE,TRUE)))
)
## export.eps(hh("logi/figure/c.eps"))

## 4                  
attach(lymph)
apx1 <- seq(20,170,30)
apx2 <- seq(50,200,30)
apx12 <- seq(20,200,30)
  
## frequencies within acid.ph range for X.ray==0 and X.ray==1 and 
nodes.freq0 <- table(cut(acid.ph[X.ray==0], apx12), nodes[X.ray==0])
nodes.freq0
nodes.prop0 <- nodes.freq0[,"1"]/(nodes.freq0[,"0"]+nodes.freq0[,"1"])
nodes.prop0

nodes.freq1 <- table(cut(acid.ph[X.ray==1], apx12), nodes[X.ray==1])
nodes.freq1
nodes.prop1 <- nodes.freq1[,"1"]/(nodes.freq1[,"0"]+nodes.freq1[,"1"])
nodes.prop1
detach("lymph")

## cg.eps.gz
print(position=c(0,.5, 1,1), more=FALSE,  ## top 50%
xyplot(nodes.j ~ acid.ph | X.ray, data=lymph,
       panel=function(x,y,subscripts,...) {
         panel.xyplot(x,y,...)
         n <- if.R(s=get("n", frame=sys.parent()),
                   r=panel.number())
         if.R(s={}, r={segments <- lsegments})
         if(n==1) segments(apx1, nodes.prop0, apx2, nodes.prop0)
         if(n==2) segments(apx1, nodes.prop1, apx2, nodes.prop1)
       },
       main=list("nodes ~ acid.ph | X.ray", cex=1.6),
## above is necessary, below makes it prettier
       scales=list(cex=1, alternating=FALSE),
       xlab=list(cex=1.4), ylab=list(cex=1.4),
       par.strip.text=list(cex=1.4), cex=1.1,
       strip=function(...) strip.default(..., style=1, strip.names=c(TRUE,TRUE)))
)
## export.eps(hh("logi/figure/cg.eps"))


## 5                  
## logit with simplified model
lymph1.glm <- glm(nodes ~ acid.ph + X.ray, data=lymph, family=binomial)
anova(lymph1.glm, test="Chisq")
summary(lymph1.glm)$coef
summary(lymph1.glm)$coef[,statistic.name]^2  ## approx = ChiSquare values

logit.p.hat <- predict.glm(lymph1.glm, type="link")
p.hat <- predict.glm(lymph1.glm, type="response")
odds.hat <- p.hat/(1-p.hat)

lhat <- cbind(lymph, p.hat=p.hat, odds.hat=odds.hat, logit.p.hat=logit.p.hat)
names(lhat)
lhat.srt <- lhat[order(lhat$X.ray,lhat$acid.ph),]


## 7
## f.EPS.gz
print(position=c(0,.5, 1,1), more=FALSE,  ## top 50%
xyplot(nodes.j ~ acid.ph | X.ray, data=lhat.srt,
       panel=function(x,y,subscripts,...) {
         panel.xyplot(x,y, ..., pch=16)
         panel.xyplot(x=lhat.srt$acid.ph[subscripts],
                      y=lhat.srt$p.hat[subscripts], type="l")
       },
       main=list("observed and predicted probability(nodes)", cex=1.6),
## above is necessary, below makes it prettier
       scales=list(cex=1, alternating=FALSE),
       xlab=list(cex=1.4), ylab=list(cex=1.4),
       par.strip.text=list(cex=1.4), cex=1.1,
       strip=function(...) strip.default(..., style=1, strip.names=c(TRUE,TRUE)))
)
## export.eps(hh("logi/figure/f.eps"))


## d.EPS.gz
pp <- if.R(r=c(0,0,1,1),
           s=c(0,-.1, 1,1.1))
print(position=pp, more=FALSE,  ## top 120%
splom(~lhat[,c(1,2,3,8,9,10)],
      panel=panel.superpose, group=lhat$X.ray,
      main=list("p.hat", cex=1.6),
## above is necessary, below makes it prettier
      superpanel=if.R(s=panel.pairs.hh, r=panel.pairs),   ## HH
      panel.cex=1.4, subpanel.scales=list(cex=1),
      cex=1, pscales=3,
      key=list(space="right",
        text=list(levels(lhat$X.ray), adj=1),
        points=tpgss,
        border=1,
        title="X.ray", cex.title=1.2,
        cex=1))
)
## export.eps(hh("logi/figure/d.eps"))


## 6                  
## e.EPS.gz
pp <- if.R(r=c(0,0,1,1),
           s=c(-.1,-.1, 1.1,1.1))
print(position=c(-.1,-.1, 1.1,1.1), more=FALSE,  ## top 120%
splom(~lhat[,c(1,2,  8,9,10)] | lhat$X.ray,
      main=list("p-hat | X.ray", cex=1.6),
## above is necessary, below makes it prettier
      superpanel=if.R(s=panel.pairs.hh, r=panel.pairs),   ## HH
      panel.cex=1.3, subpanel.scales=list(cex=.8),
      pscales=3,
      par.strip.text=list(cex=1.4), cex=1,
      strip=function(...) strip.default(..., style=1, strip.names=c(TRUE,TRUE)))
)
## export.eps(hh("logi/figure/e.eps"))


## x                  
### need to discuss
##  unusual point at X.ray=0 and acid.ph=187
##  bring in rest of variables


## logit with full model
lymph2.glm <- glm(nodes ~ X.ray + stage + grade + age + acid.ph,
                  data=lymph, family=binomial)
anova(lymph2.glm, test="Chisq")
summary(lymph2.glm)$coef
summary(lymph2.glm)$coef[,statistic.name]^2  ## approx = ChiSquare values

logit.p.hat <- predict.glm(lymph2.glm, type="link")
p.hat <- predict.glm(lymph2.glm, type="response")
odds.hat <- p.hat/(1-p.hat)

lhat <- cbind(lymph, p.hat=p.hat, odds.hat=odds.hat, logit.p.hat=logit.p.hat)

## i.EPS.gz (not displayed in chapter)
pp <- if.R(r=c(0,0,1,1),
           s=c(0,-.1, 1,1.1))
print(position=pp, more=FALSE,  ## top 120%
splom(~lhat[,c(1,2,3,8,9,10)],
      panel=panel.superpose, group=lhat$X.ray,
      main=list("p.hat", cex=1.6),
## above is necessary, below makes it prettier
      superpanel=if.R(s=panel.pairs.hh, r=panel.pairs),   ## HH
      panel.cex=1.4, subpanel.scales=list(cex=1),
      cex=1, pscales=3,
      key=list(space="right",
        text=list(levels(lhat$X.ray), adj=1),
        points=tpgss,
        border=1,
        title="X.ray", cex.title=1.2,
        cex=1))
)
## export.eps(hh("logi/figure/i.eps"))


## h.EPS.gz (not displayed in chapter)
pp <- if.R(r=c(0,0,1,1),
           s=c(-.1,-.1, 1.1,1.1))
print(position=pp, more=FALSE,  ## top 120%
splom(~lhat[,c(1,2,  8,9,10)] |  stage * grade * X.ray, data=lhat,
      layout=c(4,2),
      main=list("p-hat | stage * grade *X.ray", cex=1.6),
## above is necessary, below makes it prettier
      pscales=3,
      par.strip.text=list(cex=1.4), cex=1,
      strip=function(...) strip.default(..., style=1, strip.names=c(TRUE,TRUE)))
)
## export.eps(hh("logi/figure/h.eps"))


lhat.srt <- lhat[order(lhat$X.ray,lhat$stage,lhat$grade,lhat$acid.ph),]
## g.EPS.gz (not displayed in chapter)
print(position=c(0,.1, 1,1), more=FALSE,  ## top 90%
xyplot(nodes.j ~ acid.ph |  stage * grade * X.ray, data=lhat.srt,
       layout=c(4,2),
       panel=function(x,y,subscripts,...) {
         panel.xyplot(x,y, ..., pch=16)
         panel.xyplot(x=lhat.srt$acid.ph[subscripts],
                      y=lhat.srt$p.hat[subscripts], type="l")},
       main=
       list("observed and predicted probability(nodes) (extra wiggles from age)",
            cex=1.2),
## above is necessary, below makes it prettier
       scales=list(cex=1, alternating=FALSE),
       xlab=list(cex=1.4), ylab=list(cex=1.4),
       between=list(y=1),
       par.strip.text=list(cex=1.4), cex=1.1,
       strip=function(...)strip.default(..., style=1, strip.names=c(TRUE,TRUE)))
)
## export.eps(hh("logi/figure/g.eps"))


## 8                  
## logit with full model except age
lymph3.glm <- glm(nodes ~ X.ray + stage + grade + acid.ph,
                  data=lymph, family=binomial)
anova(lymph3.glm, test="Chisq")
summary(lymph3.glm)$coef
summary(lymph3.glm)$coef[,statistic.name]^2  ## approx = ChiSquare values

logit.p.hat <- predict.glm(lymph3.glm, type="link")
p.hat <- predict.glm(lymph3.glm, type="response")
odds.hat <- p.hat/(1-p.hat)

lhat <- cbind(lymph, p.hat=p.hat, odds.hat=odds.hat, logit.p.hat=logit.p.hat)

lhat.srt <- lhat[order(lhat$X.ray,lhat$stage,lhat$grade,lhat$acid.ph),]
## j.EPS.gz
print(position=c(0,.1, 1,1), more=FALSE,  ## top 90%
xyplot(nodes.j ~ acid.ph |  stage * grade * X.ray, data=lhat.srt,
       layout=c(4,2),
       panel=function(x,y,subscripts,...) {
         panel.xyplot(x,y, ..., pch=16)
         panel.xyplot(x=lhat.srt$acid.ph[subscripts],
                      y=lhat.srt$p.hat[subscripts], type="l")},
       main=list("observed and predicted probability(nodes) (-age)", cex=1.6),
## above is necessary, below makes it prettier
       scales=list(cex=1, alternating=FALSE),
       xlab=list(cex=1.4), ylab=list(cex=1.4),
       between=list(y=1),
       par.strip.text=list(cex=1.4), cex=1.1,
       strip=function(...)strip.default(..., style=1, strip.names=c(TRUE,TRUE)))
)
## export.eps(hh("logi/figure/j.eps"))


## 9                  
## logit with alternate specification of simplified model
lymph1a.glm <- glm(nodes ~ X.ray + acid.ph, data=lymph, family=binomial)
anova(lymph1a.glm, test="Chisq")
summary(lymph1a.glm)$coef
summary(lymph1a.glm)$coef[,statistic.name]^2  ## approx = ChiSquare values

## logit with alternate specification of simplified model with nesting
lymph1b.glm <- glm(nodes ~ X.ray/acid.ph, data=lymph, family=binomial)
anova(lymph1b.glm, test="Chisq")
summary(lymph1b.glm)$coef
summary(lymph1b.glm)$coef[,statistic.name]^2  ## approx = ChiSquare values


## common and idiosyncratic scaling
## k.EPS.gz
print(split=c(1,1,1,2), more=TRUE,  ## bottom half
xyplot(logit.p.hat ~ acid.ph | X.ray, data=lhat.srt,
       scales=list(relation="free", cex=1),
       main=list("b. logit.phat ~ acid.ph | X.ray; relation='free'", cex=1.6),
## above is necessary, below makes it prettier
       xlab=list(cex=1.4), ylab=list(cex=1.4),
       par.strip.text=list(cex=1.4), cex=1.1,
       strip=function(...) strip.default(..., style=1, strip.names=c(TRUE,TRUE)))
)

print(split=c(1,2,1,2), more=FALSE, ## top half
xyplot(logit.p.hat ~ acid.ph | X.ray, data=lhat.srt,
       scales=list(relation="same", cex=1, alternating=FALSE),
       main=list("a. logit.phat ~ acid.ph | X.ray; relation='same'", cex=1.6),
## above is necessary, below makes it prettier
       xlab=list(cex=1.4), ylab=list(cex=1.4),
       par.strip.text=list(cex=1.4), cex=1.1,
       strip=function(...) strip.default(..., style=1, strip.names=c(TRUE,TRUE)))
)
## export.eps(hh("logi/figure/k.eps"))




## lymph2.s
## requires a prior run of first 21 lines of lymph.s

lymph3.glm <- glm(nodes ~ X.ray + stage + grade + acid.ph,  
                  data=lymph, family=binomial)
anova(lymph3.glm, test="Chisq")
summary(lymph3.glm)$coef
summary(lymph3.glm)$coef[,statistic.name]^2  ## approx = ChiSquare values

logit.p.hat <- predict.glm(lymph3.glm, type="link")
p.hat <- predict.glm(lymph3.glm, type="response")
odds.hat <- p.hat/(1-p.hat)

lhat <- cbind(lymph, p.hat=p.hat, odds.hat=odds.hat, logit.p.hat=logit.p.hat)
names(lhat)
lhat.srt <- lhat[order(lhat$X.ray,lhat$acid.ph),]

pp <- if.R(r=c(0,0,1,1),
           s=c(0,-.1, 1,1.1))
print(position=pp, more=FALSE,  ## top 120%
splom(~lhat[,c(1,2,3,8,9,10)],
      panel=panel.superpose, group=lhat$X.ray,
      main=list("p.hat", cex=1.6),
## above is necessary, below makes it prettier
      superpanel=if.R(s=panel.pairs.hh, r=panel.pairs),   ## HH
      panel.cex=1.4, subpanel.scales=list(cex=1),
      cex=1, pscales=3,
      key=list(space="right",
        text=list(as.character(unique(lhat$X.ray)), adj=1),
        points=tpgss,
        border=1,
        title="X.ray", cex.title=1.2,
        cex=1))
)
## export.eps(hh("logi/figure/m.eps"))

pp <- if.R(r=c(0,0,1,1),
           s=c(-.1,-.1, 1.1,1.1))
print(position=pp, more=FALSE,  ## top 120%
splom(~lhat[,c(1,2,  8,9,10)] | X.ray, data=lhat,
      main=list("p-hat | X.ray", cex=1.6),
## above is necessary, below makes it prettier
      superpanel=if.R(s=panel.pairs.hh, r=panel.pairs),   ## HH
      panel.cex=1.3, subpanel.scales=list(cex=.8),
      pscales=3,
      par.strip.text=list(cex=1.4), cex=1,
      strip=function(...) strip.default(..., style=1, strip.names=c(TRUE,TRUE)))
)
## export.eps(hh("logi/figure/n.eps"))

xyplot(nodes.j ~ acid.ph |  X.ray + stage + grade, data=lhat.srt,
       layout=c(4,2),
       panel=function(x,y,subscripts,...) {
         panel.xyplot(x,y, ..., pch=16)
         panel.xyplot(x=lhat.srt$acid.ph[subscripts],
                      y=lhat.srt$p.hat[subscripts], type="l")
       },
       main=list("observed and predicted probability(nodes)", cex=1.6),
## above is necessary, below makes it prettier
       scales=list(cex=1, alternating=FALSE),
       xlab=list(cex=1.4), ylab=list(cex=1.4),
       par.strip.text=list(cex=1.4), cex=1.1,
       strip=function(...) strip.default(..., style=1, strip.names=c(TRUE,TRUE)))
## export.eps(hh("logi/figure/p8.eps"))

xyplot(nodes.j ~ acid.ph |  X.ray + stage + grade, data=lhat.srt,
       layout=c(4,2),
       panel=function(x,y,subscripts,...) {
         panel.xyplot(x,y, ..., pch=16)
         cplx <- if.R(s=par()$usr[1:2],
                      r=current.panel.limits()$xlim)
         new.x <- seq(cplx[1], cplx[2], length=51)
         new.y <-
           predict.glm(lymph3.glm, type="response",
                       newdata=data.frame(
                         cbind(acid.ph=new.x,
                         lhat.srt[subscripts[rep(1,length(new.x))],c("X.ray","stage","grade")])))
         panel.xyplot(x=new.x, y=new.y, type="l")
       },
       main=list("observed and predicted probability(nodes)", cex=1.6),
## above is necessary, below makes it prettier
       scales=list(cex=1, alternating=FALSE),
       xlab=list(cex=1.4), ylab=list(cex=1.4),
       par.strip.text=list(cex=1.4), cex=1.1,
       strip=function(...) strip.default(..., style=1, strip.names=c(TRUE,TRUE)))
## export.eps(hh("logi/figure/p8a.eps"))




## work/logit-da.s
lymph1.glm <- glm(nodes ~ acid.ph + X.ray,
                  data=lymph, family=binomial)
anova(lymph1.glm, test="Chisq")
summary(lymph1.glm)$coef



## work/logit-ga.s
lymph2.glm <- glm(nodes ~ X.ray + stage + grade + age + acid.ph,
                  data=lymph, family=binomial)
anova(lymph2.glm, test="Chisq")
summary(lymph2.glm)$coef




## icu.s
### rmh Hosmer and Lemeshow ICU data
data(icu)

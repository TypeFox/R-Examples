### R code from vignette source '~/WindowsC/HOME/rmh/hh.e2/hh2/logi.tex'

###################################################
### code chunk number 1: logi.tex:12-13
###################################################
library(HH)


###################################################
### code chunk number 2: logi.tex:56-80
###################################################
## hhpdf("logit-plot.pdf", width=8.5, height=5.5)
p <- seq(0,1,length=51)
y <- logit(p)

print(position=c(0, .1, .4, .9),
      panel.width=list(2, "in"), panel.height=list(4, "in"),
      more=TRUE,
      xyplot(y ~ p, type="l", main=list("y = logit(p)", cex=1.5),
             col="darkblue",
             xlab=list(cex=1.3), ylab=list(cex=1.3, rot=0),
             scales=list(cex=1.3, x=list(at=c(0, .25, .5, .75, 1),
                                    labels=c("0.0","","0.5","","1.0"))))
      )

print(position=c(.4, .2, 1, .75),
      panel.width=list(4, "in"), panel.height=list(2, "in"),
      more=FALSE,
      xyplot(p ~ y, type="l", main=list("p = antilogit(y)", cex=1.5),
             col="darkblue",
             xlab=list(cex=1.3), ylab=list(cex=1.3, rot=0),
             scales=list(cex=1.3, y=list(at=c(0, .25, .5, .75, 1),
                                    labels=c("0.0","","0.5","","1.0"))))
      )
## hhdev.off()


###################################################
### code chunk number 3: logi.tex:195-209
###################################################
data(spacshu)
Jitter.factor <- .08
## hhpdf("spaceshuttle-a.pdf", width=6.5, height=3)
xyplot(jitter(damage, factor=Jitter.factor*8) ~ tempF,
       ylab="damage",
       data=spacshu,
       scales=list(y=list(at=c(0,1))),
       cex=.75, pch=16, col="darkblue",
       main="a. observed",
       panel=function(...) {
         panel.xyplot(...)
         panel.abline(h=0:1, lty=2)
       })
## hhdev.off()


###################################################
### code chunk number 4: logi.tex:212-232
###################################################
temp.range.left <- seq(50,80,5)
temp.range.right <- seq(55,85,5)
temp.range.both <- seq(50,85,5)
damage.freq <- table(cut(spacshu$tempF, temp.range.both), spacshu$damage)
damage.prop <- damage.freq[,"1"] / (damage.freq[,"0"]+damage.freq[,"1"])

## hhpdf("spaceshuttle-b.pdf", width=6.5, height=3)
xyplot(jitter(damage, factor=Jitter.factor*8) ~ tempF,
       ylab="damage",
       data=spacshu,
       scales=list(y=list(at=c(0,1))),
       cex=.75, pch=16, col="darkblue",
       main="b. observed and sectioned proportions",
       panel=function(...) {
         panel.xyplot(...)
         panel.segments(temp.range.left, damage.prop,
                        temp.range.right, damage.prop)
         panel.abline(h=0:1, lty=2)
       })
## hhdev.off()


###################################################
### code chunk number 5: logi.tex:246-261
###################################################
## hhpdf("spaceshuttle-c.pdf", width=6.5, height=3)
xyplot(jitter(damage, factor=Jitter.factor*8) ~ tempF,
       ylab="damage",
       data=spacshu,
       scales=list(y=list(at=c(0,1))),
       cex=.75, pch=16, col="darkblue",
       main="c. observed and sectioned proportions\nappropriate temperature scale",
       xlim=c(30, 85),
       panel=function(...) {
         panel.xyplot(...)
         panel.segments(temp.range.left, damage.prop,
                        temp.range.right, damage.prop)
         panel.abline(h=0:1, lty=2)
       })
## hhdev.off()


###################################################
### code chunk number 6: logi.tex:264-276
###################################################
## hhcapture("spaceshuttle-glm.Rout", '
spacshu.bin.glm <- glm(damage ~ tempF, data=spacshu, family=binomial)
spacshu.bin.glm
anova(spacshu.bin.glm, test="Chi")
coef(summary(spacshu.bin.glm))

## prediction on response scale, in this case (0,1).
## leading to Figure spaceshuttle-d.pdf, Panel d
spacshu.pred <-
  interval(spacshu.bin.glm, newdata=data.frame(tempF=30:85),
           type="response")
## ')


###################################################
### code chunk number 7: logi.tex:279-301
###################################################
## hhpdf("spaceshuttle-d.pdf", width=6.5, height=3)
xyplot(jitter(damage, factor=Jitter.factor*4) ~ tempF, data=spacshu,
       ylab="proportion damaged",
       main="d. glm logit fit with pi, estimating p(damage in one ring)",
       xlim=c(30, 85),
       scales=list(y=list(at=c(0,1))),
       cex=.75, pch=16, col="darkblue",
       panel=function(...) {
	 panel.xyplot(...)
	 panel.lines(x=30:85, y=spacshu.pred[,"fit"],
                     lty=1,
                     col="darkblue")
	 panel.lines(x=30:85, y=spacshu.pred[,"pi.low"],
                     lty=4,
                     col="mediumblue")
	 panel.lines(x=30:85, y=spacshu.pred[,"pi.hi"],
                     lty=4,
                     col="mediumblue")
         panel.abline(h=0:1, lty=2)
       }
)
## hhdev.off()


###################################################
### code chunk number 8: logi.tex:305-332
###################################################
spacshu6 <- data.matrix(spacshu)
dimnames(spacshu6) <- NULL
dim(spacshu6) <- c(6,23,2)
spacshu6 <- data.frame(damage=apply(spacshu6[,,2],2,sum), tempF=spacshu6[1,,1])

## hhpdf("spaceshuttle-e.pdf", width=6.5, height=3)
xyplot(jitter(damage, factor=Jitter.factor*3) ~ tempF, data=spacshu6,
       ylab="number of damaged rings",
       main="e. glm logit fit with pi, estimating number of damaged rings",
       xlim=c(30,85), ylim=c(-1,8),
       scales=list(y=list(at=c(0,6))),
       cex=.75, pch=4, col="darkblue",
       panel=function(...) {
	 panel.xyplot(...)
	 panel.lines(x=30:85, y=6*spacshu.pred[,"fit"],
                     lty=1,
                     col="darkblue")
	 panel.lines(x=30:85, y=6*spacshu.pred[,"pi.low"],
                     lty=4,
                     col="mediumblue")
	 panel.lines(x=30:85, y=6*spacshu.pred[,"pi.hi"],
                     lty=4,
                     col="mediumblue")
         panel.abline(h=c(0,6), lty=2)
       }
)
## hhdev.off()


###################################################
### code chunk number 9: logi.tex:407-419
###################################################
## hhpdf("spaceshuttle-f.pdf", width=6.5, height=3)
xyplot(jitter(damage, factor=Jitter.factor*8) ~ tempF,
       ylab="damage",
       data=spacshu[spacshu[,"damage"]==1,],
       scales=list(y=list(at=c(0,1))),
       cex=.75, pch=16, col="darkblue",
       main="f. observed damaged O-rings",
       panel=function(...) {
         panel.xyplot(...)
         panel.abline(h=0:1, lty=2)
       })
## hhdev.off()


###################################################
### code chunk number 10: logi.tex:492-497
###################################################
## hhcapture("three-scales.Rout", '
p.hat <- predict(spacshu.bin.glm, type="response")
odds.hat <- p.hat/(1-p.hat)
logit.p.hat <- log(odds.hat)
## ')


###################################################
### code chunk number 11: logi.tex:523-544
###################################################
## transformations
p <- spacshu.pred[,"fit"]
tmp <- cbind.data.frame(tempF=30:85,
                        p=p,
                        odds=p/(1-p),
                        logit.p=logit(p),
                        which="fit", stringsAsFactors=FALSE)
p.obs <- spacshu6[,1]/6
tmp2 <- cbind.data.frame(spacshu6[,2,drop=FALSE],
                         p=p.obs,
                         odds=p.obs/(1-p.obs),
                         logit.p=logit(p.obs),
                         which="data", stringsAsFactors=FALSE)
tmp2$logit.p[tmp2$p==0] <- min(tmp2$logit.p[tmp2$p!=0])-2.5  ## pull infinity in
tmp.both <- rbind(tmp2,tmp)

## This splom is not in the book.  The top three panels in the left row are in the next figure.
splom(~tmp.both[,1:4], type=c("p","l"), groups=tmp.both$which,
      distribute.type=TRUE,
      pch=19, panel.cex=1.4, col="darkblue",
      cex=.9, pscales=3, xlab=NULL)


###################################################
### code chunk number 12: logi.tex:548-558
###################################################
## hhpdf("spaceshuttle-logit-xysplom.pdf", width=4.5, height=6)
xyplot(p + odds + logit.p ~ tempF, data=tmp.both,
       type=c("p","l"), groups=rep(tmp.both$which,3), pch=19, col="darkblue",
       distribute.type=TRUE,
       layout=c(1, 3),
       ylim=list(c(0, 1), c(0, 5), c(-4, 2)),
       par.strip.text=list(cex=1.4),
       scales=list(cex=1, y=list(relation="free", tick.number=4), alternating=FALSE),
       between=list(y=1.5), xlab="Temperature Fahrenheit", ylab="")
## hhdev.off()


###################################################
### code chunk number 13: logi.tex:600-610
###################################################
## hhcapture("spaceshuttle-pred-link.Rout", '
## prediction on link scale, in this case (-Inf, Inf)
## leading to Figure spaceshuttle-g.pdf Panel g
spacshu.pred.link <-
  interval(spacshu.bin.glm, newdata=data.frame(tempF=30:85),
           type="link")

cbind(tempF=30:85, round(spacshu.pred.link, digits=2))[c(1:3,54:56),]
cbind(tempF=30:85, round(spacshu.pred, digits=2))[c(1:3,54:56),]
## ')


###################################################
### code chunk number 14: logi.tex:627-651
###################################################
## hhpdf("spaceshuttle-g.pdf", width=6.5, height=3)
## approximate logit(1) with 4.1, and logit(0) with -4.1
jld <- jitter(ifelse(spacshu$damage==1, 4.1, -4.1),8) ## pull infinity in

xyplot(jld ~ tempF, data=spacshu,
       ylab="logit(proportion) damaged",
       main="g. glm logit fit with logit(pi), estimating logit(p(damage in one ring))",
       xlim=c(30,85), ylim=c(-18, 12),
       cex=.75, pch=16, col="darkblue",
       panel=function(...) {
	 panel.xyplot(...)
	 panel.lines(x=30:85, y=spacshu.pred.link[,"fit"],
                     lty=1,
                     col="darkblue")
	 panel.lines(x=30:85, y=spacshu.pred.link[,"pi.low"],
                     lty=4,
                     col="mediumblue")
	 panel.lines(x=30:85, y=spacshu.pred.link[,"pi.hi"],
                     lty=4,
                     col="mediumblue")
       }

)
## hhdev.off()


###################################################
### code chunk number 15: logi.tex:878-891
###################################################
data(budworm)
col.bud <- likertColor(2, colorFunctionOption = "default")

## hhpdf("budworm-data.pdf", width=7, height=4.5)
xyplot(numdead ~ ldose, data=budworm, groups=sex,
       pch=c("F","M"), cex=1.5, col=col.bud, type="b", lty=2,
       xlab="log dose",
       ylab="number dead",
       ylab.right="proportion dead",
       par.settings=list(clip=list(panel=FALSE), layout.widths=list(axis.right=1.4))) +
  layer(panel.axis("right", at=seq(0,20,5), labels=format(seq(0,20,5)/20, 2), outside=TRUE)) +
  layer(panel.abline(h=c(0, 20), lty=3, col="grey60"))
## hhdev.off()


###################################################
### code chunk number 16: logi.tex:910-927
###################################################
## hhcapture("budworm-glm.Rout", '
SF <- cbind(numdead=budworm$numdead,
            numalive = 20 - budworm$numdead)

## model with interaction term for sex and logdose, from VR
budworm.lg <-
  glm(SF ~ sex*ldose,
      data=budworm,
      family = binomial)
anova(budworm.lg, test="Chisq")

## model with no interaction term
budworm.lg0 <- glm(SF ~ sex + ldose - 1,
                   data=budworm,
                   family = binomial)
anova(budworm.lg0, test="Chisq")
## ')


###################################################
### code chunk number 17: logi.tex:944-967
###################################################
## hhpdf("budworm-predict.pdf", width=7, height=4.5)

ldose=seq(-1, 6, 0.1)
data.predict <- data.frame(ldose=c(ldose, ldose),
                           sex=factor(rep(c("F","M"), each=length(ldose))))
data.predict$p.hat <- predict(budworm.lg0, data.predict, type = "response")

Ap <- xyplot(p.hat ~ ldose, groups=sex,
             data=data.predict, type="l",
             xlab="log dose", ylab="probability of death", col=col.bud,
             scales=list(rot=0, y=list(at=0:4/4)))

budworm$p.hat <- predict.glm(budworm.lg0, type="response")
budworm$p <- budworm$numdead / 20
Bp <- xyplot(p ~ ldose, groups=sex, type="b", lty=2,
             data=budworm, pch=c("F","M"), cex=1.5, col=col.bud)

update(Ap + Bp,
       main=list("observed points and\nfitted logistic regression", cex=1)) +
       layer(panel.abline(h=c(0, 1), lty=3, col="grey60"))


## hhdev.off()


###################################################
### code chunk number 18: logi.tex:993-1000
###################################################
## hhcapture("logit-LD.Rout", '
## LD25 LD50 LD75
xp.M <- MASS::dose.p(budworm.lg0, cf = c(2,3), p = 1:3/4)
xp.M
xp.F <- MASS::dose.p(budworm.lg0, cf = c(1,3), p = 1:3/4)
xp.F
## ')


###################################################
### code chunk number 19: logi.tex:1012-1022
###################################################
## hhpdf("budworm-LD.pdf", width=7, height=4.5)
update(Ap + Bp +
  layer(panel.abline(h=1:3/4, col="gray60", lty=3)) +
  layer(panel.segments(xp.M, 0,   xp.M, .15, col="gray60", lty=3)) +
  layer(panel.segments(xp.F, .83, xp.F, 1, col="gray60", lty=3)) +
  layer(panel.segments(xp.M, c(.15, .37, .62), xp.M, 1, col=col.bud[2], lty=2)) +
  layer(panel.segments(xp.F, 0, xp.F, c(.37, .62, .83), col=col.bud[1], lty=2)),
  main=list("observed points and\nfitted logistic regression", cex=1)) +
  layer(panel.abline(h=c(0, 1), lty=3, col="grey60"))
## hhdev.off()


###################################################
### code chunk number 20: logi.tex:1039-1055
###################################################
## hhpdf("budworm-AB.pdf", width=3, height=5.5)
A <- xyplot(p.hat + odds(p.hat) + logit(p.hat) ~ ldose, groups=sex,
            data=data.predict, type="l",
            xlab="log dose", ylab=NULL, col=col.bud,
            scales=list(relation="free", rot=0),
            layout=c(1,3))

B <- xyplot(p + odds(p) + logit(p) ~ ldose, groups=sex,
            data=budworm, pch=c("F","M"), cex=1.5, col=col.bud,
            scales=list(relation="free"),
            layout=c(1,3))

update(combineLimits.trellisvector(A + B),
       scales=list(y=list(relation="free")),
       between=list(y=1))
## hhdev.off()


###################################################
### code chunk number 21: logi.tex:1058-1073
###################################################
## hhpdf("budworm-ABprime.pdf", width=3, height=5.5)
budworm$oddsp  <-  odds(budworm$p)
budworm$logitp <- logit(budworm$p)
budworm$oddsp[6]  <- A$y.limits[[2]][2]
budworm$logitp[6:7] <- A$y.limits[[3]][2:1]

Bprime <- xyplot(p + oddsp + logitp ~ ldose, groups=sex,
            data=budworm, pch=c("F","M"), cex=1.5, col=col.bud,
            scales=list(relation="free"),
            layout=c(1,3))

update(combineLimits.trellisvector(A + Bprime),
       scales=list(y=list(relation="free")),
       between=list(y=1))
## hhdev.off()


###################################################
### code chunk number 22: logi.tex:1176-1199
###################################################
## hhpdf("logi-f-logit-a.pdf", width=10, height=6)
data(lymph)
col2 <- likertColor(2, colorFunctionOption = "default")[2:1]

useOuterStripsT2L1(
xyplot(age ~ acid.ph | grade * stage * X.ray, data=lymph,
       group=nodes, pch=levels(lymph$nodes), col=col2, cex=2.2,
       layout=c(4, 2),
       ## above is necessary, below makes it prettier
       main=list("age ~ acid.ph | grade * stage * X.ray, group=nodes", cex=1.6),
       aspect=1,
       between=list(x=c(.5, 1, .5), y=1),
       scales=list(cex=1, alternating=FALSE),
       xlab=list(cex=1.4), ylab=list(cex=1.4),
       par.strip.text=list(cex=1.6),
       key=list(space="right",
         text=list(levels(lymph$nodes), cex=1.5, adj=1, col=col2),
         columns=2,
         border=1,
         title="nodes", cex.title=1.25,
         cex=1))
)
## hhdev.off()


###################################################
### code chunk number 23: logi.tex:1234-1241
###################################################
old.width <- options(width=65)
## hhcapture("logit-j.Rout", '
lymph3.glm <- glm(nodes ~ X.ray + stage + grade + acid.ph,
                  data=lymph, family=binomial)
anova(lymph3.glm, test="Chisq")
## ')
options(old.width)


###################################################
### code chunk number 24: logi.tex:1267-1316
###################################################
## hhpdf("p8.pdf", width=8.5, height=6)
logit.p.hat <- predict.glm(lymph3.glm, type="link")
p.hat <- predict.glm(lymph3.glm, type="response")
odds.hat <- p.hat/(1-p.hat)
lhat <- cbind(lymph, p.hat=p.hat, odds.hat=odds.hat, logit.p.hat=logit.p.hat)
lhat.sort <- lhat[with(lhat, order(X.ray, stage, grade, acid.ph)),]
lhat.sort$Xsg <- with(lhat.sort, interaction(X.ray, stage, grade))
lhat.sort$Xsg <- factor(lhat.sort$Xsg, levels=unique(lhat.sort$Xsg))


p8d <-
xyplot(nodes.j ~ acid.ph | grade + stage + X.ray, data=lhat.sort,
       layout=c(4,2), between=list(x=c(.5, 1, .5), y=1.5),
       groups=lhat.sort$nodes,
       pch=levels(lhat.sort$nodes), col=col2, cex=2.2,
       main=list("jittered observed and predicted probability(nodes)", cex=1.6),
## above is necessary, below makes it prettier
       scales=list(cex=1, alternating=FALSE, y=list(at=c(0, .25, .5, .75, 1))),
       xlab=list(cex=1.4), ylab=list(cex=1.4),
       par.strip.text=list(cex=1.4),
       strip=strip.custom(strip.names=c(TRUE,TRUE)),
       key=list(space="right",
         text=list(levels(lymph$nodes), cex=1.5, adj=1, col=col2),
         columns=2,
         border=1,
         title="nodes", cex.title=1.25,
         cex=1))

## p8d

ul35 <- unique(lymph[, 3:5])
ul35 <- ul35[with(ul35, order(X.ray, stage, grade)),]
ul35$Xsg <- with(ul35, interaction(X.ray, stage, grade))
ul35$Xsg <- factor(ul35$Xsg, levels=unique(ul35$Xsg))
old.warn <- options(warn=-1)  ## "row names were found from a short variable and have been discarded"
tmp <- lapply(1:8, function(i) cbind(acid.ph=29:198, ul35[i,]))
options(old.warn)
tmp2 <- do.call("rbind", tmp)
tmp2$nodes.hat <- predict.glm(lymph3.glm, type="response", newdata=tmp2)

p8e <- xyplot(nodes.hat ~ acid.ph |  grade + stage + X.ray, data=tmp2,
              layout=c(4,2),
              type="l", col="black")
## p8e

Na <- layer(panel.abline(h=c(0,1), lty=2, col="gray30"))

useOuterStripsT2L1(p8d + p8e + Na)
## hhdev.off()


###################################################
### code chunk number 25: logi.tex:1367-1396
###################################################
## hhpdf("lymph-m.pdf", width=4, height=5.5)
col8 <- brewer.pal(8, 'Dark2')

## values plotting symbol based on stage and full fitted lines
DDg <-
xyplot(p.hat + odds.hat + logit.p.hat ~ acid.ph, data=lhat.sort,
       group=Xsg, type="p",
       layout=c(1, 3), scales=list(relation="free", rot=0),
       ## above is necessary, below makes it prettier
       pch=c("0","0","1","1","0","0","1","1"), cex=1.2, col=col8,
       key=list(space="right",
         border=1,
         text=list(levels(lhat.sort$Xsg)),
         points=list(pch=c("0","0","1","1","0","0","1","1"), col=col8, cex=1.2),
         title="X.ray.  \n  stage.\n     grade", cex.title=1, lines.title=1.2,
         cex=1,
         rev=TRUE))
DDgcl <-
combineLimits(as.matrix(update(
  DDg,
  ylab=list(c("p.hat", "odds.hat", "logit.p.hat"), rot=0),
  strip=FALSE, between=list(y=1))))
DDgcll <-
xyplot(nodes.hat + odds(nodes.hat) + logit(nodes.hat) ~ acid.ph, data=tmp2,
       group=Xsg, type="l", col=col8, lwd=2,
       layout=c(1, 3), scales=list(relation="free", rot=0))
DDGcla <- DDgcl + DDgcll
DDGcla
## hhdev.off()


###################################################
### code chunk number 26: logi.tex:1409-1459
###################################################
EEg0 <-
xyplot(p.hat + odds.hat + logit.p.hat ~ acid.ph, subset=(X.ray==0), data=lhat.sort,
       group=Xsg, type="p",
       layout=c(1, 3), scales=list(relation="free", rot=0),
       ## above is necessary, below makes it prettier
       pch=c("0","0","1","1","0","0","1","1"), cex=1.2, col=col8,
       key=list(space="right",
         border=1,
         text=list(levels(lhat.sort$Xsg)),
         points=list(pch=c("0","0","1","1","0","0","1","1"), col=col8, cex=1.2),
         title="X.ray.\n stage.\n grade", cex.title=1.2,
         cex=1,
         rev=TRUE))
EEg0l <-
xyplot(nodes.hat + odds(nodes.hat) + logit(nodes.hat) ~ acid.ph, subset=(X.ray==0), data=tmp2,
       group=Xsg, type="l", col=col8, lwd=2,
       layout=c(1, 3), scales=list(relation="free", rot=0))
EEg0a <- EEg0 + EEg0l
EEg0a

EEg1 <-
xyplot(p.hat + odds.hat + logit.p.hat ~ acid.ph, subset=(X.ray==1), data=lhat.sort,
       group=Xsg, type="p",
       layout=c(1, 3), scales=list(relation="free", rot=0),
       ## above is necessary, below makes it prettier
       pch=c("0","0","1","1","0","0","1","1"), cex=1.2, col=col8,
       key=list(space="right",
         border=1,
         text=list(levels(lhat.sort$Xsg)),
         points=list(pch=c("0","0","1","1","0","0","1","1"), col=col8, cex=1.2),
         title="X.ray.\n stage.\n grade", cex.title=1.2,
         cex=1,
         rev=TRUE))
EEg1l <-
xyplot(nodes.hat + odds(nodes.hat) + logit(nodes.hat) ~ acid.ph, subset=(X.ray==1), data=tmp2,
       group=Xsg, type="l", col=col8, lwd=2,
       layout=c(1, 3), scales=list(relation="free", rot=0))
EEg1a <- EEg1 + EEg1l
EEg1a

## hhpdf("lymph-n.pdf", width=6, height=5.5)
EEga <-
update(cbind(EEg0, EEg1),
       strip=FALSE, strip.left=FALSE,
       between=list(x=1, y=1),
       scales=list(x=list(at=pretty(lhat.sort$acid.ph), relation="same", alternating=1),
         y=list(relation="free")),
       xlab.top=c("X-ray: 0", "X-ray: 1"))
combineLimits(EEga) + cbind(EEg0l, EEg1l)
## hhdev.off()


###################################################
### code chunk number 27: logi.tex:1475-1489
###################################################
## hhpdf("lymph-mn.pdf", width=8, height=5.5)
DEEga <-
update(cbind(update(DDgcl, scales=list(at=NULL)), EEg0, EEg1),
       strip=FALSE, strip.left=FALSE,
       between=list(x=c(2, .5), y=1),
       scales=list(x=list(at=pretty(lhat.sort$acid.ph), relation="same", alternating=1),
         y=list(relation="free")),
       xlab.top=c("X-ray: both   ", "   X-ray: 0", "X-ray: 1"))
combineLimits(DEEga) + cbind(as.vector(DDgcll), EEg0l, EEg1l)
## hhdev.off()

## notice that latticeExtra:::`+.trellis` is is not commutative with HH:::`cbind.trellis`
## cbind(EEg0a, EEg1a)        ## doesn't work
## EEga + cbind(EEg0l, EEg1l) ## works


###################################################
### code chunk number 28: logi.tex:1525-1547
###################################################
## hhpdf("logi-c.pdf", width=6, height=5.5)
## ignore grade and stage, simplification for exposition
N0 <-
xyplot(nodes.j ~ acid.ph | X.ray, data=lymph,
       groups=nodes,
       pch=levels(lymph$nodes), col=col2, cex=2.2,
       ## above is necessary, below makes it prettier
       layout=c(1,2), between=list(y=1),
       scales=list(cex=1, alternating=FALSE, y=list(at=c(0, .25, .5, .75, 1))),
       xlab=list(cex=1.4), ylab=list(cex=1.4),
       par.strip.text=list(cex=1.4),
       strip=FALSE,
       strip.left=strip.custom(strip.names=c(TRUE,TRUE)),
       key=list(space="right",
         text=list(levels(lymph$nodes), adj=1, col=col2),
         columns=2,
         border=1,
         title="nodes", cex.title=1.25,
         cex=1))
N0.Na <- N0 + Na
print(N0.Na, panel.width=list(3.5,"in"))
## hhdev.off()


###################################################
### code chunk number 29: logi.tex:1568-1579
###################################################
## hhpdf("logit-cg.pdf", width=6, height=5.5)
apx12 <- seq(20,200,20)

nodes.freq <- table(cut(lymph$acid.ph, apx12), lymph$X.ray, lymph$nodes)
nodes.prop <- nodes.freq[,,"1"] / (nodes.freq[,,"0"]+nodes.freq[,,"1"])

Ns <- layer(panel.segments(apx12[-10], nodes.prop[,panel.number()],
                           apx12[-01], nodes.prop[,panel.number()]))
N0.Na.Ns <- N0 + Na + Ns
print(N0.Na.Ns, panel.width=list(3.5,"in"))
## hhdev.off()


###################################################
### code chunk number 30: logi.tex:1623-1628
###################################################
## hhcapture("logit-f.Rout", '
lymph1.glm <- glm(nodes ~ X.ray + acid.ph, data=lymph, family=binomial)
anova(lymph1.glm, test="Chisq")
summary(lymph1.glm)$coef
## ')


###################################################
### code chunk number 31: logi.tex:1643-1654
###################################################
col.xray <- brewer.pal(12, 'Paired')[c(10,12)] ## col8[6:7] ## "black" ## col8[1:2]
## hhpdf("logit-d.pdf", width=6, height=5.5)
tmpX <- lapply(c("0","1"), function(i) data.frame(acid.ph=29:198, X.ray=i))
tmpX2 <- do.call("rbind", tmpX)
tmpX2$nodes.hat <- predict.glm(lymph1.glm, type="response", newdata=tmpX2)
Np <- xyplot(nodes.hat ~ acid.ph | X.ray, data=tmpX2,
             groups=X.ray,
             type="l", layout=c(1,2), col=col.xray)
N0.Na.Np <- N0 + Na + Np
print(N0.Na.Np, panel.width=list(3.5,"in"))
## hhdev.off()


###################################################
### code chunk number 32: logi.tex:1668-1703
###################################################
## hhpdf("logit-e.pdf", width=6, height=7)
logit.p.hat <- predict.glm(lymph1.glm, type="link")
p.hat <- predict.glm(lymph1.glm, type="response")
odds.hat <- p.hat/(1-p.hat)

lhat2 <- cbind(lymph, p.hat=p.hat, odds.hat=odds.hat, logit.p.hat=logit.p.hat)
lhat2.srt <- lhat2[order(lhat2$X.ray, lhat2$acid.ph),]

DD2g <-
xyplot(p.hat + odds.hat + logit.p.hat ~ acid.ph, data=lhat2.srt,
       group=X.ray, type="p",
       layout=c(1, 3), scales=list(relation="free", rot=0, cex=1),
       ## above is necessary, below makes it prettier
       pch=c("0","1"), cex=2, col=col.xray,
       key=list(space="right",
         border=1, columns=2,
         points=list(pch=c("0","1"), col=col.xray, cex=1.5),
         title="X.ray", cex.title=1.2, lines.title=1.4,
         cex=1.4,
         rev=TRUE))
DD2gcl <-
combineLimits(as.matrix(update(
  DD2g, xlab=list(cex=1.2),
  ylab=list(c("p.hat", "odds.hat", "logit.p.hat"), rot=0, cex=1.2),
  strip=FALSE, between=list(y=1))))

DD2gcll <-
xyplot(nodes.hat + odds(nodes.hat) + logit(nodes.hat) ~ acid.ph, data=tmpX2,
       group=X.ray, type="l", col=col.xray,
       layout=c(1, 3), scales=list(relation="free", rot=0))
DD2Gcla <- DD2gcl + DD2gcll
## DD2Gcla
DD2Gcla.update <- update(DD2Gcla, ylim=list(c(0,1), c(0, 10), c(-2, 8)), layout=c(1,3))
print(DD2Gcla.update, panel.width=list(3.5,"in"))
## hhdev.off()


###################################################
### code chunk number 33: logi.tex:1732-1736
###################################################
## hhcapture("logit-k.Rout", '
lymph1Xa.glm <- glm(nodes ~ X.ray * acid.ph, data=lymph, family=binomial)
anova(lymph1Xa.glm, test="Chisq")
## ')


###################################################
### code chunk number 34: logi.tex:1750-1758
###################################################
## hhpdf("logit-f.pdf", width=6, height=5.5)
tmpX2$nodes.hatXa <- predict.glm(lymph1Xa.glm, type="response", newdata=tmpX2)
NpXa <- xyplot(nodes.hatXa ~ acid.ph | X.ray, data=tmpX2,
             groups=X.ray,
             type="l", layout=c(1,2), col=col.xray)
N0.Na.NpXa <- N0 + Na + NpXa
print(N0.Na.NpXa, panel.width=list(3.5,"in"))
## hhdev.off()


###################################################
### code chunk number 35: logi.tex:1786-1821
###################################################
## hhpdf("logit-g.pdf", width=6, height=7)
logit.p.hat <- predict.glm(lymph1Xa.glm, type="link")
p.hat <- predict.glm(lymph1Xa.glm, type="response")
odds.hat <- p.hat/(1-p.hat)

lhat2Xa <- cbind(lymph, p.hat=p.hat, odds.hat=odds.hat, logit.p.hat=logit.p.hat)
lhat2Xa.srt <- lhat2Xa[order(lhat2Xa$X.ray, lhat2Xa$acid.ph),]

DD2gXa <-
xyplot(p.hat + odds.hat + logit.p.hat ~ acid.ph, data=lhat2Xa.srt,
       group=X.ray, type="p",
       layout=c(1, 3), scales=list(relation="free", rot=0, cex=1),
       ## above is necessary, below makes it prettier
       pch=c("0","1"), cex=2, col=col.xray,
       key=list(space="right",
         border=1, columns=2,
         points=list(pch=c("0","1"), col=col.xray, cex=1.5),
         title="X.ray", cex.title=1, lines.title=1.4,
         cex=1.4,
         rev=TRUE))
DD2gXacl <-
combineLimits(as.matrix(update(
  DD2gXa,
  ylab=list(c("p.hat", "odds.hat", "logit.p.hat"), rot=0, cex=1.2),
  strip=FALSE, between=list(y=1))))

DD2gXacll <-
xyplot(nodes.hatXa + odds(nodes.hatXa) + logit(nodes.hatXa) ~ acid.ph, data=tmpX2,
       group=X.ray, type="l", col=col.xray,
       layout=c(1, 3), scales=list(relation="free", rot=0))
DD2GXacla <- DD2gXacl + DD2gXacll
## DD2GXacla
DD2GXacla.update <- update(DD2GXacla, ylim=list(c(0,1), c(0, 10), c(-2, 8)))
print(DD2GXacla.update, panel.width=list(3.5,"in"))
## hhdev.off()


###################################################
### code chunk number 36: logi.tex:2003-2025
###################################################
## hhpdf("logit-k.pdf", width=8, height=7)

tmp <-
update(DEEga[2:3, 3], layout=c(2,1),
       ylab=list(DEEga$ylab[[1]][3], rot=0),
       xlab.top=DEEga$xlab.top[2:3],
       scales=list(y=list(limits=c(-3, 4), at=-3:4, relation="same")),
       main="Common Scaling: Easy to Compare Panels",
       par.settings=list(layout.heights=NULL, layout.widths=NULL))

print(position=c(0,.49,1,1), more=TRUE,  ## top
tmp +
   cbind(as.vector(DDgcll), EEg0l, EEg1l)
)

print(position=c(0,0,1,.47), more=FALSE, ## bottom
update(tmp,
       scales=list(y=list(relation="free")),
       main="Separate Scaling: Difficult to Compare Panels") +
   cbind(as.vector(DDgcll), EEg0l, EEg1l)
)
## hhdev.off()



### R code from vignette source 'HE-examples.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
set.seed(1071)
old.opts <- options(width=85, digits=5, useFancyQuotes = FALSE, continue="  ")
library(heplots)
library(candisc)


###################################################
### code chunk number 2: plastic-mod
###################################################
plastic.mod <- lm(cbind(tear, gloss, opacity) ~ rate*additive, data=Plastic)
Anova(plastic.mod, test.statistic="Roy")


###################################################
### code chunk number 3: plastic-univar
###################################################
Anova(update(plastic.mod, tear ~ .))
Anova(update(plastic.mod, gloss ~ .))
Anova(update(plastic.mod, opacity ~ .))


###################################################
### code chunk number 4: plastic1a (eval = FALSE)
###################################################
## # Compare evidence and effect scaling 
## colors = c("red", "darkblue", "darkgreen", "brown")
## heplot(plastic.mod, size="evidence", col=colors, cex=1.25)
## heplot(plastic.mod, size="effect", add=TRUE, lwd=4, term.labels=FALSE, col=colors)


###################################################
### code chunk number 5: plastic1b (eval = FALSE)
###################################################
## ## add interaction means
## intMeans <- termMeans(plastic.mod, 'rate:additive', abbrev.levels=2)
## #rownames(intMeans) <- apply(expand.grid(c('Lo','Hi'), c('Lo', 'Hi')), 1, paste, collapse=':')
## points(intMeans[,1], intMeans[,2], pch=18, cex=1.2, col="brown")
## text(intMeans[,1], intMeans[,2], rownames(intMeans), adj=c(0.5,1), col="brown")
## lines(intMeans[c(1,3),1], intMeans[c(1,3),2], col="brown")
## lines(intMeans[c(2,4),1], intMeans[c(2,4),2], col="brown")


###################################################
### code chunk number 6: plastic1
###################################################
#      <<plastic1a>>
# Compare evidence and effect scaling 
colors = c("red", "darkblue", "darkgreen", "brown")
heplot(plastic.mod, size="evidence", col=colors, cex=1.25)
heplot(plastic.mod, size="effect", add=TRUE, lwd=4, term.labels=FALSE, col=colors)
#      <<plastic1b>>
## add interaction means
intMeans <- termMeans(plastic.mod, 'rate:additive', abbrev.levels=2)
#rownames(intMeans) <- apply(expand.grid(c('Lo','Hi'), c('Lo', 'Hi')), 1, paste, collapse=':')
points(intMeans[,1], intMeans[,2], pch=18, cex=1.2, col="brown")
text(intMeans[,1], intMeans[,2], rownames(intMeans), adj=c(0.5,1), col="brown")
lines(intMeans[c(1,3),1], intMeans[c(1,3),2], col="brown")
lines(intMeans[c(2,4),1], intMeans[c(2,4),2], col="brown")


###################################################
### code chunk number 7: plastic-mod
###################################################
plastic.mod


###################################################
### code chunk number 8: plastic-tests
###################################################
print(linearHypothesis(plastic.mod, c("rateHigh", "additiveHigh"), title="Main effects"), SSP=FALSE)

print(linearHypothesis(plastic.mod, c("rateHigh", "additiveHigh", "rateHigh:additiveHigh"), title="Groups"), SSP=FALSE)


###################################################
### code chunk number 9: plastic2
###################################################
heplot(plastic.mod, hypotheses=list("Group" = 
       c("rateHigh", "additiveHigh", "rateHigh:additiveHigh ")),
       col=c(colors, "purple"),
       lwd=c(2, 3, 3, 3, 2), cex=1.25)
heplot(plastic.mod, hypotheses=list("Main effects" = 
       c("rateHigh", "additiveHigh")), add=TRUE,
       col=c(colors, "darkgreen"), cex=1.25)


###################################################
### code chunk number 10: plastic1-HE3D (eval = FALSE)
###################################################
## colors = c("pink", "darkblue", "darkgreen", "brown")
## heplot3d(plastic.mod, col=colors)


###################################################
### code chunk number 11: MJdata
###################################################
str(MockJury)


###################################################
### code chunk number 12: MJdata1
###################################################
table(MockJury$Attr)
table(MockJury$Attr, MockJury$Crime)


###################################################
### code chunk number 13: jury.mod1
###################################################
(jury.mod1 <- lm( cbind(phyattr, happy, independent, sophisticated) ~ Attr, data=MockJury))
Anova(jury.mod1, test="Roy")


###################################################
### code chunk number 14: jury-mod1-HE
###################################################
heplot(jury.mod1, main="HE plot for manipulation check")


###################################################
### code chunk number 15: jury-mod1-pairs
###################################################
pairs(jury.mod1)


###################################################
### code chunk number 16: jury-can1a
###################################################
jury.can <- candisc(jury.mod1)
jury.can


###################################################
### code chunk number 17: jury-can1
###################################################
opar <- par(xpd=TRUE)
heplot(jury.can, prefix="Canonical dimension", main="Canonical HE plot")
par(opar)


###################################################
### code chunk number 18: jury-mod2
###################################################
# influence of Attr of photo and nature of crime on Serious and Years
jury.mod2 <- lm( cbind(Serious, Years) ~ Attr * Crime, data=MockJury)
Anova(jury.mod2, test="Roy")


###################################################
### code chunk number 19: jury-mod2-HE
###################################################
heplot(jury.mod2)


###################################################
### code chunk number 20: jury-mod3-HE
###################################################
# stepdown test (ANCOVA), controlling for Serious
jury.mod3 <- lm( Years ~ Serious + Attr * Crime, data=MockJury)
t(coef(jury.mod3))
Anova(jury.mod3)


###################################################
### code chunk number 21: jury-mod3-eff
###################################################
library(effects)
jury.eff <- allEffects(jury.mod3)
plot(jury.eff, ask=FALSE)


###################################################
### code chunk number 22: skulls1
###################################################
data(Skulls)
str(Skulls)
table(Skulls$epoch)


###################################################
### code chunk number 23: skulls2
###################################################
# make shorter labels for epochs
Skulls$epoch <- factor(Skulls$epoch, labels=sub("c","",levels(Skulls$epoch)))
# assign better variable labels
vlab <- c("maxBreadth", "basibHeight", "basialLength", "nasalHeight")


###################################################
### code chunk number 24: skulls3
###################################################
means <- aggregate(cbind(mb, bh, bl, nh) ~ epoch, data=Skulls, FUN=mean)[,-1]
rownames(means) <- levels(Skulls$epoch)
means


###################################################
### code chunk number 25: skulls4
###################################################
pairs(means, vlab,
      panel = function(x, y) {
          text(x, y, levels(Skulls$epoch))
          lines(x,y)
      })


###################################################
### code chunk number 26: skulls-bwplot
###################################################
library(lattice)
library(reshape2)
sklong <- melt(Skulls, id="epoch")

bwplot(value ~ epoch | variable, data=sklong, scales="free", 
	ylab="Variable value", xlab="Epoch",
	strip=strip.custom(factor.levels=paste(vlab, " (", levels(sklong$variable), ")", sep="")),
	panel = function(x,y, ...) {
		panel.bwplot(x, y, ...)
		panel.linejoin(x,y, col="red", ...)
	}) 


###################################################
### code chunk number 27: skulls5
###################################################
# fit manova model
sk.mod <- lm(cbind(mb, bh, bl, nh) ~ epoch, data=Skulls)
Manova(sk.mod)


###################################################
### code chunk number 28: skulls5a
###################################################
coef(sk.mod)


###################################################
### code chunk number 29: skulls6
###################################################
coef(sk.mod)["epoch.L",]
print(linearHypothesis(sk.mod, "epoch.L"), SSP=FALSE) # linear component


###################################################
### code chunk number 30: skulls6
###################################################
print(linearHypothesis(sk.mod, c("epoch.Q", "epoch.C", "epoch^4")), SSP=FALSE)


###################################################
### code chunk number 31: skulls-HE-pairs
###################################################
pairs(sk.mod, variables=c(1,4,2,3),
	hypotheses=list(Lin="epoch.L", NonLin=c("epoch.Q", "epoch.C", "epoch^4")), 
	var.labels=vlab[c(1,4,2,3)])


###################################################
### code chunk number 32: skulls-HE3D (eval = FALSE)
###################################################
## heplot3d(sk.mod, hypotheses=list(Lin="epoch.L", Quad="epoch.Q", 
##                             NonLin=c("epoch.Q", "epoch.C", "epoch^4")), 
## 	col=c("pink", "blue"))


###################################################
### code chunk number 33: skulls-can1
###################################################
library(candisc)
sk.can <- candisc(sk.mod)
sk.can


###################################################
### code chunk number 34: skulls-can2
###################################################
heplot(sk.can, prefix="Canonical dimension")


###################################################
### code chunk number 35: rohwer-some
###################################################
some(Rohwer,n=6)


###################################################
### code chunk number 36: rohwer-separate
###################################################
rohwer.ses1 <- lm(cbind(SAT, PPVT, Raven) ~ n + s + ns + na + ss, data=Rohwer, subset=SES=="Hi")
Anova(rohwer.ses1)

rohwer.ses2 <- lm(cbind(SAT, PPVT, Raven) ~ n + s + ns + na + ss, data=Rohwer, subset=SES=="Lo")
Anova(rohwer.ses2)


###################################################
### code chunk number 37: rohwer-coef
###################################################
coef(rohwer.ses1)
coef(rohwer.ses2)


###################################################
### code chunk number 38: rohwer-HE1
###################################################
heplot(rohwer.ses1, ylim=c(40,110),col=c("red", "black"), lwd=2, cex=1.2)
heplot(rohwer.ses2, add=TRUE, col=c("blue", "black"), grand.mean=TRUE, error.ellipse=TRUE, lwd=2, cex=1.2)
means <- aggregate(cbind(SAT,PPVT)~SES, data=Rohwer,  mean)
text(means[,2], means[,3], labels=means[,1], pos=3, cex=2, col=c("red", "blue"))


###################################################
### code chunk number 39: rohwer-mod
###################################################
# MANCOVA, assuming equal slopes
rohwer.mod <- lm(cbind(SAT, PPVT, Raven) ~ SES + n + s + ns + na + ss, 
                 data=Rohwer)
Anova(rohwer.mod)


###################################################
### code chunk number 40: rohwer-mod-test
###################################################
(covariates <- rownames(coef(rohwer.mod))[-(1:2)])
Regr<-linearHypothesis(rohwer.mod, covariates)
print(Regr, digits=5, SSP=FALSE)


###################################################
### code chunk number 41: rohwer-HE2
###################################################
colors <- c("red", "blue", rep("black",5), "#969696")
heplot(rohwer.mod, col=colors,
      hypotheses=list("Regr" = c("n", "s", "ns", "na", "ss")),
      cex=1.5, lwd=c(2, rep(3,5), 4),
      main="(SAT, PPVT, Raven) ~ SES + n + s + ns + na + ss")


###################################################
### code chunk number 42: rohwer-HE3
###################################################
heplot(rohwer.mod, col=colors,  variables=c(1,3),
      hypotheses=list("Regr" = c("n", "s", "ns", "na", "ss")),
      cex=1.5, lwd=c(2, rep(3,5), 4),
      main="(SAT, PPVT, Raven) ~ SES + n + s + ns + na + ss")


###################################################
### code chunk number 43: rohwer-HE3 (eval = FALSE)
###################################################
## pairs(rohwer.mod, col=colors,
##       hypotheses=list("Regr" = c("n", "s", "ns", "na", "ss")),
##       cex=1.3, lwd=c(2, rep(3,5), 4))


###################################################
### code chunk number 44: rohwer-HE3D (eval = FALSE)
###################################################
## colors <- c("pink", "blue", rep("black",5), "#969696")
## heplot3d(rohwer.mod, col=colors,
## 	hypotheses=list("Regr" = c("n", "s", "ns", "na", "ss")))


###################################################
### code chunk number 45: rohwer-mod2
###################################################
rohwer.mod2 <- lm(cbind(SAT, PPVT, Raven) ~ SES * (n + s + ns + na + ss),
                  data=Rohwer)
Anova(rohwer.mod2)


###################################################
### code chunk number 46: rohwer-mod2-test
###################################################
(coefs <- rownames(coef(rohwer.mod2)))
print(linearHypothesis(rohwer.mod2, coefs[grep(":", coefs)]), SSP=FALSE)


###################################################
### code chunk number 47: rohwer-HE4
###################################################
colors <- c("red", "blue", rep("black",5), "#969696")
heplot(rohwer.mod2, col=c(colors, "brown"), 
      terms=c("SES", "n", "s", "ns", "na", "ss"), 
      hypotheses=list("Regr" = c("n", "s", "ns", "na", "ss"),
                      "Slopes" = coefs[grep(":", coefs)]))


###################################################
### code chunk number 48: hern1
###################################################
Hern.mod <- lm(cbind(leave, nurse, los) ~ age + sex +  pstat +  build + cardiac + resp, 
               data=Hernior)
Anova(Hern.mod) 


###################################################
### code chunk number 49: hern.rsq
###################################################
Hern.summary <- summary(Hern.mod)
unlist(lapply(Hern.summary, function(x) x$r.squared))


###################################################
### code chunk number 50: hern2
###################################################
# test overall regression
Regr <- linearHypothesis(Hern.mod, c("age", "sexm", "pstat", "build", "cardiac", "resp"))
print(Regr, digits=5, SSP=FALSE)


###################################################
### code chunk number 51: hern-pairs
###################################################
clr <- c("red", "darkgray", "blue", "darkgreen", "magenta", "brown", "black")
vlab <- c("LeaveCondition\n(leave)", "NursingCare\n(nurse)", "LengthOfStay\n(los)")
hyp <- list("Regr" = c("age", "sexm", "pstat", "build", "cardiac", "resp"))
pairs(Hern.mod, hypotheses=hyp, col=clr, var.labels=vlab)


###################################################
### code chunk number 52: hern-can1
###################################################
Hern.canL <- candiscList(Hern.mod)


###################################################
### code chunk number 53: hern.can2 (eval = FALSE)
###################################################
## plot(Hern.canL)


###################################################
### code chunk number 54: hern-can-all (eval = FALSE)
###################################################
## plot(Hern.canL, term="pstat")
## plot(Hern.canL, term="build")
## plot(Hern.canL, term="age")
## plot(Hern.canL, term="cardiac")


###################################################
### code chunk number 55: hern-can-pstat
###################################################
plot(Hern.canL, term="pstat")


###################################################
### code chunk number 56: hern-can-build
###################################################
plot(Hern.canL, term="build")


###################################################
### code chunk number 57: hern-can-age
###################################################
plot(Hern.canL, term="age")


###################################################
### code chunk number 58: hern-can-cardiac
###################################################
plot(Hern.canL, term="cardiac")


###################################################
### code chunk number 59: grades1
###################################################
data(SocGrades)
grades.mod <- lm(cbind(midterm1, midterm2, final, eval) ~ 
	class + sex + gpa + boards + hssoc + pretest, data=SocGrades)
Anova(grades.mod, test="Roy")


###################################################
### code chunk number 60: grades1
###################################################
grades.mod2 <- update(grades.mod, . ~ .^2)
Anova(grades.mod2, test="Roy")


###################################################
### code chunk number 61: grades3
###################################################
grades.mod3 <- update(grades.mod, . ~ . + class:sex - hssoc - pretest)
Anova(grades.mod3, test="Roy")


###################################################
### code chunk number 62: grades-pairs
###################################################
pairs(grades.mod3)


###################################################
### code chunk number 63: grades-HE3D (eval = FALSE)
###################################################
## heplot3d(grades.mod3, wire=FALSE)


###################################################
### code chunk number 64: grades4
###################################################
# calculate canonical results for all terms
grades.can <- candiscList(grades.mod3)
# extract canonical R^2s
unlist(lapply(grades.can, function(x) x$canrsq))


###################################################
### code chunk number 65: grades-can-class
###################################################
# plot class effect in canonical space
 op <- par(xpd=TRUE)
 heplot(grades.can, term="class", scale=4, fill=TRUE, var.col="black", var.lwd=2)
 par(op)


###################################################
### code chunk number 66: grades-can-all (eval = FALSE)
###################################################
## plot(grades.can, term="sex")
## plot(grades.can, term="gpa")


###################################################
### code chunk number 67: grades-can-sex
###################################################
plot(grades.can, term="sex")


###################################################
### code chunk number 68: grades-can-gpa
###################################################
plot(grades.can, term="gpa")


###################################################
### code chunk number 69: tidyup
###################################################
options(old.opts )



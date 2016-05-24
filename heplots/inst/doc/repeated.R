### R code from vignette source 'repeated.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
set.seed(1071)
old.opts <- options(width=85, digits=6, useFancyQuotes = FALSE, continue="  ", prompt="R> ")
library(heplots)
library(lattice)
library(nlme)
library(car)
# avoid deprecated warning
if (car2 <- packageDescription("car")[["Version"]] >= 2) data.ellipse <- dataEllipse



###################################################
### code chunk number 2: voc0
###################################################
some(VocabGrowth,5)


###################################################
### code chunk number 3: voc1
###################################################
voc <- reshape(VocabGrowth, direction="long", varying=list(grade=1:4), timevar="Grade", v.names="Vocabulary")
boxplot(Vocabulary ~ Grade, data=voc, col="bisque",
	ylab="Vocabulary", main="Vocabulary Growth data")
abline(lm(Vocabulary ~ as.numeric(Grade), data=voc), col="red")
means <- tapply(voc$Vocabulary, voc$Grade, mean)
points(1:4, means, pch=7, col="blue")
lines(1:4, means, col="blue", lwd=2)


###################################################
### code chunk number 4: voc2
###################################################
(Vocab.mod <- lm(cbind(grade8,grade9,grade10,grade11) ~ 1, data=VocabGrowth))


###################################################
### code chunk number 5: voc3
###################################################
(Vocab.aov0 <- Anova(Vocab.mod, type="III"))


###################################################
### code chunk number 6: voc-Anova
###################################################
Vocab.aov0$SSP     # H matrix
Vocab.aov0$SSPE    # E matrix


###################################################
### code chunk number 7: markH0
###################################################
mark.H0 <- function(x=0, y=0, cex=2, pch=19, col="green3", lty=2, pos=2) {
  points(x,y, cex=cex, col=col, pch=pch)
  text(x,y, expression(H[0]), col=col, pos=pos)
  if (lty>0) abline(h=y, col=col, lty=lty)
  if (lty>0) abline(v=x, col=col, lty=lty)
}


###################################################
### code chunk number 8: voc4
###################################################
heplot(Vocab.mod, terms="(Intercept)", type="III")
mark.H0(0,0)
title(expression(paste("Multivariate test of ", H[0], " : ", bold(mu)==0)))


###################################################
### code chunk number 9: voc5
###################################################
idata <-data.frame(grade=ordered(8:11))
(Vocab.aov <- Anova(Vocab.mod, idata=idata, idesign=~grade, type="III"))


###################################################
### code chunk number 10: voc6
###################################################
trends <- as.matrix(VocabGrowth) %*% poly(8:11, degree=3)
colnames(trends)<- c("Linear", "Quad", "Cubic")

# test all trend means = 0 == Grade effect
within.mod <- lm(trends ~ 1)

Anova(within.mod, type="III")


###################################################
### code chunk number 11: voc8
###################################################
op <- par(mfrow=c(1,2))
data.ellipse(trends[,1:2], xlim=c(-4,8), ylim=c(-3,3), levels=0.68,
	main="(a) Data ellipse ")
mark.H0(0,0)

heplot(within.mod, terms="(Intercept)", col=c("red", "blue"), type="III",
	term.labels="Grade", , xlim=c(-4,8), ylim=c(-3,3),
	main="(b) HE plot for Grade effect")
mark.H0(0,0)
par(op)


###################################################
### code chunk number 12: obk1
###################################################
library("car")	# for OBrienKaiser data

# simplified analysis of OBrienKaiser data, collapsing over hour
OBK <- OBrienKaiser
OBK$pre <- rowMeans(OBK[,3:7])
OBK$post <- rowMeans(OBK[,8:12])
OBK$fup <- rowMeans(OBK[,13:17])
# remove separate hour scores
OBK <- OBK[,-(3:17)]


###################################################
### code chunk number 13: obk2
###################################################
table(OBK$gender, OBK$treatment)


###################################################
### code chunk number 14: obk3
###################################################
contrasts(OBK$gender)
contrasts(OBK$treatment)


###################################################
### code chunk number 15: obk-MANOVA
###################################################
# MANOVA model
mod.OBK <- lm(cbind(pre, post, fup) ~  treatment*gender,  data=OBK)
mod.OBK


###################################################
### code chunk number 16: obk4
###################################################
session <- ordered(c("pretest", "posttest", "followup"),
    levels=c("pretest", "posttest", "followup"))
contrasts(session) <- matrix(c(-1,  1, 0,
                                0, -1, 1), ncol=2)
session
idata <- data.frame(session)


###################################################
### code chunk number 17: obk4
###################################################
# Multivariate tests for repeated measures
aov.OBK <- Anova(mod.OBK, idata=idata, idesign=~session, test="Roy")
aov.OBK


###################################################
### code chunk number 18: obk5
###################################################
# All multivariate tests
summary(aov.OBK, univariate=FALSE)
# Univariate tests for repeated measures
summary(aov.OBK, multivariate=FALSE)


###################################################
### code chunk number 19: obk-HE1
###################################################
# HE plots for Between-S effects
heplot(mod.OBK, hypotheses=c("treatment1", "treatment2"),
  col=c("red", "black", "blue", "brown", "gray40", "gray40"),
  hyp.labels=c("(A,B)-Control", "A-B"),
  main="Between-S effects and contrasts"
	)
lines(c(3,7), c(3,7), col="green")


###################################################
### code chunk number 20: obk-3d (eval = FALSE)
###################################################
## heplot3d(mod.OBK, hypotheses=c("treatment1", "treatment2"),
## 	col=c("pink", "black", "blue", "brown", "gray40", "gray40"),
## 	hyp.labels=c("(A,B)-Control", "A-B")
## 	)
## # rotate around y axis
## play3d( rot8y <- spin3d(axis=c(0,1,0)),duration=12 )


###################################################
### code chunk number 21: obk-HE2
###################################################
pairs(mod.OBK, col=c("red", "black", "blue", "brown"))


###################################################
### code chunk number 22: obk6
###################################################
# Transform to profile contrasts for within-S effects
OBK$session.1 <- OBK$post - OBK$pre
OBK$session.2 <- OBK$fup - OBK$post
mod1.OBK <- lm(cbind(session.1, session.2) ~ treatment*gender,  data=OBK)
Anova(mod1.OBK, test="Roy", type="III")


###################################################
### code chunk number 23: obk-HE3
###################################################
# HE plots for Within-S effects
heplot(mod1.OBK,
  main="Within-S effects: Session * (Treat*Gender)",
  remove.intercept=FALSE, type="III",
  xlab="Post-Pre", ylab="Fup-Post",
  term.labels=c("session", "treatment:session", "gender:session",
                "treatment:gender:session"),
  col=c("red", "black", "blue", "brown"),
  xlim=c(-2,4), ylim=c(-2,3)
)
mark.H0(0,0)


###################################################
### code chunk number 24: obk2-1
###################################################
session <- factor(rep(c("pretest", "posttest", "followup"), c(5, 5, 5)),
    levels=c("pretest", "posttest", "followup"))
contrasts(session) <- matrix(c(-1,  1, 0,
		                            0, -1, 1), ncol=2)
hour <- ordered(rep(1:5, 3))
within <- data.frame(session, hour)


###################################################
### code chunk number 25: obk2-2
###################################################
str(within)
within


###################################################
### code chunk number 26: obk2-3
###################################################
mod.OBK2 <- lm(cbind(pre.1, pre.2, pre.3, pre.4, pre.5,
                     post.1, post.2, post.3, post.4, post.5,
                     fup.1, fup.2, fup.3, fup.4, fup.5) ~  treatment*gender,
                data=OBrienKaiser)
(aov.OBK2 <- Anova(mod.OBK2, idata=within, idesign=~session*hour, test="Roy"))


###################################################
### code chunk number 27: obk2-4
###################################################
M.session <- matrix(c(-1,  1, 0,
                       0, -1, 1), ncol=2)
rownames(M.session) <-c("pre", "post", "fup")		
colnames(M.session) <-paste("s", 1:2, sep="")		

M.hour <- matrix(c(-2, -1, 0, 1, 2,
                    2, -1, -2, -1, 1,
                    -1, 2, 0, -2, 1,
                    1, -4, 6, -4, 1), ncol=4)
rownames(M.hour) <- paste("hour", 1:5, sep="")	
colnames(M.hour) <- c("lin", "quad", "cubic", "4th")
M.hour

unit <- function(n, prefix="") {
	J <-matrix(rep(1, n), ncol=1)
	rownames(J) <- paste(prefix, 1:n, sep="")
	J
}

M.s <- kronecker( M.session, unit(5, "h"), make.dimnames=TRUE)

(M.h <- kronecker( unit(3, "s"), M.hour, make.dimnames=TRUE))

M.sh <- kronecker( M.session, M.hour, make.dimnames=TRUE)


###################################################
### code chunk number 28: obk2-5
###################################################
Y.hour <- as.matrix(OBrienKaiser[,3:17]) %*% M.h
mod.OBK2.hour <- lm(Y.hour ~  treatment*gender,  data=OBrienKaiser)


###################################################
### code chunk number 29: obk2-HE1
###################################################
labels <- c("hour", paste(c("treatment","gender","treatment:gender"),":hour", sep=""))
colors <- c("red", "black", "blue", "brown", "purple")
heplot(mod.OBK2.hour, type="III", remove.intercept=FALSE, term.labels=labels, col=colors)
mark.H0()


###################################################
### code chunk number 30: obk2-HE2
###################################################
pairs(mod.OBK2.hour, type="III", remove.intercept=FALSE, term.labels="hour", col=colors)


###################################################
### code chunk number 31: doubly
###################################################
M.measure <- diag(2)
rownames(M.measure)<- c("Y1", "Y2")
colnames(M.measure)<- c("Y1", "Y2")

kronecker(cbind(1, M.session), M.measure, make.dimnames=TRUE)


###################################################
### code chunk number 32: wl1
###################################################
table(WeightLoss$group)
some(WeightLoss)


###################################################
### code chunk number 33: wl-means
###################################################
op <- par(mfrow=c(1,2))
library("gplots")
WLlong <-reshape(WeightLoss, varying=list(WeightLoss=2:4, SelfEsteem=5:7), direction="long",
	timevar = "month",
	v.names=c("WeightLoss", "SelfEsteem"))
WLlong <- WLlong[ order(WLlong$id),]

plotmeans(WeightLoss ~ interaction(month, group), data=WLlong,
	ylab="Weight Loss", xlab="Month / Group",
	connect=list(1:3, 4:6, 7:9), pch=7, p=0.68,
	n.label=FALSE, main="Weight Loss: Group x Month",
	legends =  paste(rep(1:3, 3),c("", "Control", "", "", "Diet", "", "", "DietEx", ""),sep="\n")
	)
abline(v=c(3.5, 6.5), col="gray")

plotmeans(SelfEsteem ~ interaction(month, group), data=WLlong,
	ylab="Self Esteem", xlab="Month / Group",
	connect=list(1:3, 4:6, 7:9), pch=7, p=0.68,
	n.label=FALSE, main="Self Esteem: Group x Month",
	legends =  paste(rep(1:3, 3),c("", "Control", "", "", "Diet", "", "", "DietEx", ""),sep="\n")
	)
abline(v=c(3.5, 6.5), col="gray")
par(op)


###################################################
### code chunk number 34: wl2
###################################################
contrasts(WeightLoss$group) <- matrix(c(-2,1,1, 0, -1, 1),ncol=2)
(wl.mod<-lm(cbind(wl1,wl2,wl3,se1,se2,se3)~group, data=WeightLoss))


###################################################
### code chunk number 35: wl3
###################################################
Anova(wl.mod)


###################################################
### code chunk number 36: wl-HE1
###################################################
op <- par(mfrow=c(1,2))
heplot(wl.mod, hypotheses=c("group1", "group2"),
	xlab="Weight Loss, month 1", ylab="Weight Loss, month 2")
heplot(wl.mod, hypotheses=c("group1", "group2"), variables=4:5,
	xlab="Self Esteem, month 1", ylab="Self Esteem, month 2")
par(op)


###################################################
### code chunk number 37: wl4
###################################################
measure <- kronecker(diag(2), unit(3, 'M')/3, make.dimnames=TRUE)
colnames(measure)<- c('WL', 'SE')
measure

between <- as.matrix(WeightLoss[,-1]) %*% measure

between.mod <- lm(between ~ group, data=WeightLoss)
Anova(between.mod, test="Roy", type="III")


###################################################
### code chunk number 38: wl-HE2
###################################################
heplot(between.mod, hypotheses=c("group1", "group2"),
	xlab="Weight Loss", ylab="Self Esteem",
	col=c("red", "blue", "brown"),
	main="Weight Loss & Self Esteem: Group Effect")


###################################################
### code chunk number 39: wl5
###################################################
month <- kronecker(diag(2), poly(1:3, degree=2), make.dimnames=TRUE)
colnames(month)<- c('WL1', 'WL2', 'SE1', 'SE2')
round(month,digits=4)
trends <- as.matrix(WeightLoss[,-1]) %*% month
within.mod <- lm(trends ~ group, data=WeightLoss)
Anova(within.mod, test="Roy", type="III")


###################################################
### code chunk number 40: wl-HE3
###################################################
op <- par(mfrow=c(1,2))
heplot(within.mod, hypotheses=c("group1", "group2"), variables=c(1,3),
	xlab="Weight Loss - Linear", ylab="Self Esteem - Linear",
	type="III", remove.intercept=FALSE,
	term.labels=c("month", "group:month"),
	main="(a) Within-S Linear Effects")
mark.H0()

heplot(within.mod, hypotheses=c("group1", "group2"), variables=c(2,4),
	xlab="Weight Loss - Quadratic", ylab="Self Esteem - Quadratic",
	type="III", remove.intercept=FALSE,
	term.labels=c("month", "group:month"),
	main="(b) Within-S Quadratic Effects")
mark.H0()
par(op)


###################################################
### code chunk number 41: voc-new1
###################################################
(Vocab.mod <- lm(cbind(grade8,grade9,grade10,grade11) ~ 1, data=VocabGrowth))

idata <-data.frame(grade=ordered(8:11))

heplot(Vocab.mod, type="III", idata=idata, idesign=~grade, iterm="grade",
	main="HE plot for Grade effect")


###################################################
### code chunk number 42: obk-new0
###################################################
mod.OBK <- lm(cbind(pre, post, fup) ~  treatment*gender,  data=OBK)
session <- ordered(c("pretest", "posttest", "followup"),
    levels=c("pretest", "posttest", "followup"))
contrasts(session) <- matrix(c(-1,  1, 0,
                                0, -1, 1), ncol=2)


###################################################
### code chunk number 43: obk-new1
###################################################
idata <- data.frame(session)
heplot(mod.OBK, idata=idata, idesign=~session, iterm="session",
  col=c("red", "black", "blue", "brown"),
  main="Within-S effects: Session * (Treat*Gender)")


###################################################
### code chunk number 44: obk-new2 (eval = FALSE)
###################################################
## mod.OBK2 <- lm(cbind(pre.1, pre.2, pre.3, pre.4, pre.5,
##                      post.1, post.2, post.3, post.4, post.5,
##                      fup.1, fup.2, fup.3, fup.4, fup.5) ~  treatment*gender,
##                 data=OBrienKaiser)
## heplot(mod.OBK2, idata=within, idesign=~hour, iterm="hour")
## heplot(mod.OBK2, idata=within, idesign=~session*hour, iterm="session:hour")


###################################################
### code chunk number 45: ortho-xyplot1-code
###################################################
data("Orthodont", package="nlme")
library("lattice")
xyplot(distance ~ age|Sex, data=Orthodont, type='b', groups=Subject, pch=15:25,
    col=palette(), cex=1.3, main="Orthodont data")


###################################################
### code chunk number 46: ortho-xyplot1
###################################################
print(trellis.last.object())


###################################################
### code chunk number 47: ortho-xyplot2-code
###################################################
xyplot(distance ~ age | Sex, data = Orthodont, groups = Subject,
    main = "Pooled OLS and Individual linear regressions ~ age", type = c('g', 'r'),
    panel = function(x, y, ...) {
        panel.xyplot(x, y, ..., col = gray(0.5))
        panel.lmline(x, y, ..., lwd = 3, col = 'red')
                })


###################################################
### code chunk number 48: ortho-xyplot2
###################################################
plot(trellis.last.object())


###################################################
### code chunk number 49: Ortho.mix1
###################################################
Ortho <- Orthodont
Ortho$year <- Ortho$age - 8  # make intercept = initial status

Ortho.mix1 <- lme(distance ~ year * Sex, data=Ortho,
		random = ~ 1 + year | Subject, method="ML")
#Ortho.mix1
anova(Ortho.mix1)


###################################################
### code chunk number 50: Ortho.mix3
###################################################
Ortho.mix3 <- lme(distance ~ year*Sex + I(year^2) + I(year^3), data=Ortho,
		random = ~ 1 + year | Subject, method="ML")
anova(Ortho.mix3)


###################################################
### code chunk number 51: Ortho.anova
###################################################
anova(Ortho.mix1, Ortho.mix3)


###################################################
### code chunk number 52: Ortho-fm1-code
###################################################
grid <- expand.grid(year=0:6, Sex=c("Male", "Female"))
grid$age <- grid$year+8  # plot vs. age
fm.mix1 <-cbind(grid, distance = predict(Ortho.mix1, newdata = grid, level=0))

xyplot(distance ~ age, data=fm.mix1, groups=Sex, type="b",
    par.settings = list(superpose.symbol = list(cex = 1.2, pch=c(15,16))),
    auto.key=list(text=levels(fm.mix1$Sex), points = TRUE, x=0.05, y=0.9, corner=c(0,1)),
    main="Linear mixed model: predicted growth")


###################################################
### code chunk number 53: Ortho-fm1
###################################################
plot(trellis.last.object())


###################################################
### code chunk number 54: Ortho-fm3
###################################################
fm.mix3 <-cbind(grid, distance = predict(Ortho.mix3, newdata = grid, level=0))
print(xyplot(distance ~ age, data=fm.mix3, groups=Sex, type="b",
	par.settings = list(superpose.symbol = list(cex = 1.2, pch=c(15,16))),
	auto.key=list(text=levels(fm.mix3$Sex), points = TRUE, x=0.05, y=0.9, corner=c(0,1)),
	main="Cubic mixed model: predicted growth"))


###################################################
### code chunk number 55: ortho-prep
###################################################
library("nlme")
Orthowide <- reshape(Orthodont, v.names="distance", idvar=c("Subject", "Sex"),
	timevar="age", direction="wide")
some(Orthowide, 4)


###################################################
### code chunk number 56: ortho-mlm
###################################################
Ortho.mod <- lm(cbind(distance.8, distance.10, distance.12, distance.14) ~ Sex, data=Orthowide)

idata <- data.frame(age=ordered(seq(8,14,2)))
Ortho.aov <- Anova(Ortho.mod, idata=idata, idesign=~age)
Ortho.aov


###################################################
### code chunk number 57: ortho-HE
###################################################
op <- par(mfrow=c(1,2))
heplot(Ortho.mod, variables=c(1,4), asp=1, col=c("red", "blue"),
	xlim=c(18,30), ylim=c(18,30),
	main="Orthodont data: Sex effect")
abline(0,1, col="green")

heplot(Ortho.mod, idata=idata, idesign=~age, iterm="age", col=c("red", "blue", "brown"),
	main="Orthodont data: Within-S effects")
par(op)


###################################################
### code chunk number 58: ortho-pairs
###################################################
pairs(Ortho.mod, idata=idata, idesign=~age, iterm="age", col=c("red", "blue", "brown"))


###################################################
### code chunk number 59: ortho-nonlin-HE
###################################################
heplot(Ortho.mod, idata=idata, idesign=~age, iterm="age", col=c("red", "blue", "brown"),
	variables=c(2,3), main="Orthodont data: Nonlinear Within-S effects")


###################################################
### code chunk number 60: ortho-lh1
###################################################
linearHypothesis(Ortho.mod, hypothesis="SexFemale", idata=idata, idesign=~age, iterms="age", title="Sex:age effect")


###################################################
### code chunk number 61: ortho-lh2
###################################################
linear <- idata
contrasts(linear$age, 1) <- contrasts(linear$age)[,1]
print(linearHypothesis(Ortho.mod, hypothesis="(Intercept)",
	idata=linear, idesign=~age, iterms="age", title="Linear age"), SSP=FALSE)
print(linearHypothesis(Ortho.mod, hypothesis="SexFemale",
	idata=linear, idesign=~age, iterms="age", title="Linear Sex:age"), SSP=FALSE)


###################################################
### code chunk number 62: ortho-lh3
###################################################
nonlin <- idata
contrasts(nonlin$age, 2) <- contrasts(nonlin$age)[,2:3]
print(linearHypothesis(Ortho.mod, hypothesis="(Intercept)",
	idata=nonlin, idesign=~age, iterms="age", title="Nonlinear age"), SSP=FALSE)
print(linearHypothesis(Ortho.mod, hypothesis="SexFemale",
	idata=nonlin, idesign=~age, iterms="age", title="Nonlinear Sex:age"), SSP=FALSE)


###################################################
### code chunk number 63: tidyup
###################################################
options(old.opts )



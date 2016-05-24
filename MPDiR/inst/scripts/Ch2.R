### R code from vignette source 'LinMod/Models.Stex'

###################################################
### code chunk number 1: Models.Stex:8-27
###################################################
options(width=60, digits = 3, str = strOptions(strict.width = "cut"))
pkg <- search()[2] 
while (search()[2] !=  if(.Platform$GUI == "AQUA") "tools:RGUI" else "package:stats") {
#detach(pos = match(pkg, search()))
spkg <- strsplit( pkg, ":"  )[[1]][2] 
if (packageHasNamespace(spkg, .libPaths()[1]) )
unloadNamespace(spkg ) else detach(pos = match(pkg, search()))
 pkg <- search()[2]
 }
rm(list = ls())
gc()
library(MPDiR)
library(lattice) 
ltheme <- canonical.theme("pdf", color = FALSE) ## in-built B&W theme 
ltheme$strip.background$bg <- "grey85" ## change strip bg 
lattice.options(default.theme = ltheme) ## set as default 
formals(deparse)$width.cutoff <- 50L
assignInNamespace("deparse", deparse, "base")
rm(deparse)


###################################################
### code chunk number 2: Models.Stex:139-140 (eval = FALSE)
###################################################
## Y ~ X1 + X2 + 1


###################################################
### code chunk number 3: Models.Stex:148-149 (eval = FALSE)
###################################################
## Y ~ X1 + X2 


###################################################
### code chunk number 4: Models.Stex:156-157 (eval = FALSE)
###################################################
## Y ~ X1 + X2 - 1


###################################################
### code chunk number 5: Models.Stex:160-161 (eval = FALSE)
###################################################
## Y ~ X1 + X2 + 0


###################################################
### code chunk number 6: Models.Stex:190-191 (eval = FALSE)
###################################################
## Y ~ A + B + A:B


###################################################
### code chunk number 7: Models.Stex:198-199 (eval = FALSE)
###################################################
## Y ~ A * B


###################################################
### code chunk number 8: Models.Stex:262-263
###################################################
(A <- rep(factor(1:3), each = 2))


###################################################
### code chunk number 9: Models.Stex:271-272 (eval = FALSE)
###################################################
## Y ~ A


###################################################
### code chunk number 10: Models.Stex:277-278
###################################################
model.matrix( ~  A)


###################################################
### code chunk number 11: Models.Stex:300-301
###################################################
model.matrix( ~  A - 1)


###################################################
### code chunk number 12: Models.Stex:319-320
###################################################
model.matrix( ~  A, contrasts = list(A = "contr.sum"))


###################################################
### code chunk number 13: Models.Stex:335-336
###################################################
data(Motion)


###################################################
### code chunk number 14: Models.Stex:345-346
###################################################
str(Motion)


###################################################
### code chunk number 15: Models.Stex:359-369 (eval = FALSE)
###################################################
## library(lattice)
## xyplot(LnThresh ~ LnAge | Mtype + Sex, data = Motion,
##    xlab = "Log Age (months)",
##    ylab = "Log Threshold (contrast)",
##    panel = function(x, y, ...) {
##       panel.lmline(x, y, ...)
##       panel.loess(x, y, lty = 2, lwd = 2, col = "black",...)
##       panel.xyplot(x, y, col = "black", pch = 16, cex = 1.5, 
##       ...) },
## )


###################################################
### code chunk number 16: Models.Stex:406-419
###################################################
library(lattice)
print(
	xyplot(LnThresh ~ LnAge | Mtype + Sex, data = Motion,
	xlab = "Log Age (months)",
	ylab = "Log Threshold (contrast)", 
	panel = function(x, y, ...) {
		panel.lmline(x, y, ...)
		panel.loess(x, y, lty = 2, lwd = 3, col = "black",...)
		panel.xyplot(x, y, col = "black", pch = 16, ...)
		},
	strip = function(..., bg) strip.default(..., bg = "gray85"),
	)
)


###################################################
### code chunk number 17: Models.Stex:458-459
###################################################
Motion.lm0 <- lm(LnThresh ~ Mtype + Sex +  LnAge, data = Motion)


###################################################
### code chunk number 18: Models.Stex:487-488
###################################################
coef(summary(Motion.lm0))


###################################################
### code chunk number 19: Models.Stex:511-512
###################################################
Motion.lm1 <- update(Motion.lm0, . ~ .^2)


###################################################
### code chunk number 20: Models.Stex:524-525 (eval = FALSE)
###################################################
## anova(Motion.lm0, Motion.lm1)


###################################################
### code chunk number 21: Models.Stex:563-565
###################################################
Motion.lm2 <- update(Motion.lm0, . ~ . - Sex)
anova(Motion.lm2, Motion.lm0)


###################################################
### code chunk number 22: Models.Stex:578-581
###################################################
opar <- par(mfrow = c(2, 2))
plot(Motion.lm2)
par(opar)


###################################################
### code chunk number 23: Models.Stex:737-747
###################################################
 site <- file.path("http://vision.arc.nasa.gov/modelfest/", 
 	"data/modelfestbaselinedata.csv")
 ModelFest <- read.csv(site, header = FALSE)
 Obs <- ModelFest$V1
 ModelFest.df <- data.frame(LContSens = c(t(ModelFest[, -1])),
	Obs = rep(Obs, each = ncol(ModelFest) - 1),
	Stim = rep(paste("Stim", 1:43, sep = ""), each = 4)
	)
ModelFest.df$Stim <- with(ModelFest.df, factor(Stim, 
		levels = unique(Stim)))


###################################################
### code chunk number 24: Models.Stex:781-786
###################################################
mfGab <- subset(ModelFest.df, 
	Stim %in% paste("Stim", 1:10, sep = ""))
mfGab$Stim <- mfGab$Stim[, drop = TRUE]
SpatFreq <- c(1.12, 2^seq(1, 4.5, 0.5), 30)
mfGab$SpatFreq <- rep(SpatFreq, each = 4)


###################################################
### code chunk number 25: Models.Stex:826-829
###################################################
with(mfGab, interaction.plot(SpatFreq, Obs, LContSens, mean, 
	type = "b", xlab = "Spatial Frequency (c/deg)", 
	ylab = "Log Contrast Sensitivity", cex.lab = 1.5))


###################################################
### code chunk number 26: Models.Stex:852-854
###################################################
mfGab.lm <- lm(LContSens ~ Obs * Stim, mfGab)
anova(mfGab.lm)


###################################################
### code chunk number 27: Models.Stex:860-862
###################################################
sum(anova(mfGab.lm)[, 2]) 
sd(mfGab$LContSens)^2 * (nrow(mfGab) - 1)


###################################################
### code chunk number 28: Models.Stex:886-888
###################################################
mfGab.aov <- aov(LContSens ~ Stim + Error(Obs/Stim), mfGab)
summary(mfGab.aov)


###################################################
### code chunk number 29: Models.Stex:958-961
###################################################
library(lme4.0)
mfGab.lmer <- lmer(LContSens ~ Stim + (1 | Obs), mfGab)
print(mfGab.lmer, correlation = FALSE)


###################################################
### code chunk number 30: Models.Stex:984-986
###################################################
mfGab.lmer2 <-update(mfGab.lmer, . ~ . + (1 | Obs:Stim))
print(mfGab.lmer2, correlation = FALSE)


###################################################
### code chunk number 31: Models.Stex:994-995
###################################################
anova(mfGab.lmer2)


###################################################
### code chunk number 32: Models.Stex:1004-1005
###################################################
anova(mfGab.lmer, mfGab.lmer2)


###################################################
### code chunk number 33: Models.Stex:1068-1070
###################################################
 mfGab.mcmc <- mcmcsamp(mfGab.lmer2, n = 1000)
 HPDinterval(mfGab.mcmc)

###################################################
### mcmcsamp appears broken w/ respect to generating
### estimates of the variance components; 
###   conf. ints do not contain the estimate itself!
###   Here is an alternate method w/ lme from nlme
###################################################
unloadNamespace("lme4.0")
library(nlme)
mfGab.lme <- lme(LContSens ~ Stim, data = mfGab,
	random = ~ 1 | Obs/Stim)
intervals(mfGab.lme)

###################################################
### code chunk number 34: Models.Stex:1090-1092
###################################################
( mfGab.fe2 <- fixef(mfGab.lmer2) )
mfGab.fe2[-1] <- mfGab.fe2[1] + mfGab.fe2[-1]


###################################################
### code chunk number 35: Models.Stex:1105-1106
###################################################
str( mfGab.re2 <- ranef(mfGab.lmer2) )


###################################################
### code chunk number 36: Models.Stex:1112-1115
###################################################
mfGab.pred2 <- mfGab.fe2 + 
	rep(mfGab.re2[[2]][, 1], each = length(mfGab.fe2)) +
	ranef(mfGab.lmer2)[[1]][, 1]


###################################################
### code chunk number 37: Models.Stex:1119-1135
###################################################
Thr.mean <- t(with(mfGab, tapply(LContSens, list(Obs, Stim), mean)))
Thr.df <- stack(as.data.frame(Thr.mean))
names(Thr.df) <- c("Mean", "Obs")
Thr.df$SpatFreq <- SpatFreq
Thr.df$pred2 <- mfGab.pred2
print(
xyplot(Mean ~ SpatFreq | Obs, Thr.df, subscripts = TRUE, 
    id = Thr.df$pred2, scales = list(x = list(log = TRUE)),
    xlab = "Spatial Frequency (c/deg)",
    ylab = "Log Contrast Sensitivity",
    panel = function(x, y, subscripts, id,  ...){
	panel.xyplot(x, y)
	panel.lines(x, mfGab.fe2)
	panel.lines(x, id[subscripts], lty = 2, lwd = 2)
     } )
)


###################################################
### code chunk number 38: Models.Stex:1216-1217
###################################################
data(Chromatic)


###################################################
### code chunk number 39: Models.Stex:1233-1242 (eval = FALSE)
###################################################
## library(ggplot2)
## qplot(Age, Thresh, data = Chromatic, 
##     facets = .~ Axis, col = "black", log = "y", 
##     geom = c("point"),
##     xlab = "Age (years)", ylab = "Chromatic Threshold") + 
##     scale_x_log2(limits = 2^c(-3.5, 7.5),
##         breaks = 2^seq(-2, 6, 2), 
##         labels =   2^seq(-2, 6, 2) ) +
## 	geom_smooth(colour = "black")


###################################################
### code chunk number 40: Models.Stex:1251-1263
###################################################
library(ggplot2)
print(
qplot(Age, Thresh, data = Chromatic, 
	facets = .~ Axis,   log = "y", 
	geom =  c("point"),
	xlab = "Age (years)",
	ylab = "Chromatic Threshold") + 
	scale_x_log2(limits = 2^c(-3.5, 7.5),
		breaks = 2^seq(-2, 6, 2), 
		labels =   2^seq(-2, 6, 2) ) +
	geom_smooth(colour = "black", se = FALSE)
)


###################################################
### code chunk number 41: Models.Stex:1302-1308
###################################################
Chrom.nls <- nls(Thresh ~ 
		a[Axis] * Age^-alph + b[Axis] * Age^alph, 
		data = Chromatic, 
		start = list(a = c(8, 8, 11) * 10^-3,
			  b = c(3, 3, 8) * 10^-5, alph = 1)
			    )


###################################################
### code chunk number 42: Models.Stex:1332-1333
###################################################
summary(Chrom.nls)


###################################################
### code chunk number 43: Models.Stex:1344-1349 (eval = FALSE)
###################################################
## par(mfrow = c(1, 2))
## plot(fitted(Chrom.nls), resid(Chrom.nls))
## abline(0, 0, lty = 2)
## qqnorm(resid(Chrom.nls))
## qqline(resid(Chrom.nls))


###################################################
### code chunk number 44: Models.Stex:1361-1367
###################################################
opar <- par(mfrow = c(1, 2))
plot(fitted(Chrom.nls), resid(Chrom.nls))
abline(0, 0, lty = 2)
qqnorm(resid(Chrom.nls))
qqline(resid(Chrom.nls))
par(opar)


###################################################
### code chunk number 45: Models.Stex:1375-1376
###################################################
Chrom2.nls <- update(Chrom.nls, log10(.) ~ log10(.))


###################################################
### code chunk number 46: Models.Stex:1382-1388
###################################################
opar <- par(mfrow = c(1, 2))
plot(fitted(Chrom2.nls), resid(Chrom2.nls))
abline(0, 0, lty = 2)
qqnorm(resid(Chrom2.nls))
qqline(resid(Chrom2.nls))	
par(opar)


###################################################
### code chunk number 47: Models.Stex:1411-1415
###################################################
Chrom3.nls <- update(Chrom2.nls, . ~ 
  log10(a[Axis] * Age^-alph[Axis] + b[Axis] * Age^alph[Axis]),
  start = list(a = c(8, 8, 11) * 10^-3, 
  b = c(3, 3, 8) * 10^-5, alph = c(1, 1, 1)))


###################################################
### code chunk number 48: Models.Stex:1445-1446
###################################################
confint(Chrom2.nls)


###################################################
### code chunk number 49: Models.Stex:1462-1473
###################################################
Chromatic$fit <- 10^fitted(Chrom3.nls)
print(
qplot(Age, Thresh, data = Chromatic, 
	facets = .~ Axis, log = "y", geom =  c("point"),
	xlab = "Age (years)",
	ylab = "Chromatic Threshold") + 
	scale_x_log2(limits = 2^c(-3.5, 7.5),
		breaks = 2^seq(-2, 6, 2), 
		labels =   2^seq(-2, 6, 2) ) +
	geom_line(aes(x = Age, y = fit), size = 2, colour = rgb(0, 0, 0))
)


###################################################
### code chunk number 50: Models.Stex:1537-1541
###################################################
library(mgcv)
Chrom.gam <- gam(log10(Thresh) ~ s(Age, by = Axis) + Axis, 
	data = Chromatic)
summary(Chrom.gam)	


###################################################
### code chunk number 51: Models.Stex:1565-1567
###################################################
Chrom2.gam <- update(Chrom.gam, .  ~ s(Age) + Axis)
anova(Chrom2.gam, Chrom.gam, test = "F")	


###################################################
### code chunk number 52: Models.Stex:1572-1588
###################################################
opar <- par(mfrow = c(1, 3))
cc <- coef(Chrom2.gam)[1:3] + c(0, coef(Chrom2.gam)[1], 
	coef(Chrom2.gam)[1])
names(cc) <- c("Deutan", "Protan", "Tritan")
for (Ax in names(cc)) {
  plot(Chrom2.gam, rug = FALSE,  
	shift = cc[Ax], 
	ylim = c(-3.5, -1) - cc[Ax],
	main = Ax, log = "x",
	xlab = "Age (years)",
	ylab = "Log(Chromatic Threshold)")
  points(log10(Thresh) ~ Age, Chromatic,
	subset = Axis == Ax,
	col = rgb(0, 0, 0, 0.5))
}	
par(opar)


###################################################
### code chunk number 53: Models.Stex:1742-1744
###################################################
Chrom.glm <- glm(Thresh ~ Axis:(I(Age^-1) + Age), 
	family = Gamma(link = "identity"), data = Chromatic)


###################################################
### code chunk number 54: Models.Stex:1763-1765
###################################################
Chrom2.glm <- update(Chrom.glm, . ~ Axis + (I(Age^-1) + Age))
anova(Chrom2.glm, Chrom.glm, test = "F")


###################################################
### code chunk number 55: Models.Stex:1777-1796
###################################################
opar <- par(mfrow = c(1, 3))
xx <- 10^seq(-0.6, 2, len = 200)
for (Ax in levels(Chromatic$Axis)){
  plot(log10(Thresh) ~ Age, Chromatic, 
    subset = Axis == Ax,
    xlab = "Age (years)", 
    ylab = expression(paste(plain(Log)[10], " (Chromatic Threshold)")),
    main = Ax, log = "x", 
    ylim = c(-3.5, -1))
  Chrom.pred <- predict(Chrom.glm,
    newdata = data.frame(Age = xx, Axis = Ax),
    se.fit = TRUE)
  polygon(c(xx, rev(xx)), with(Chrom.pred, log10(c(fit + 1.96 * se.fit, 
    rev(fit - 1.96 * se.fit)))), col = "grey")
  lines(xx, log10(Chrom.pred$fit))
  points(log10(Thresh) ~ Age, Chromatic, subset = Axis == Ax,
    pch = 21, bg = "white")
}
par(opar)


###################################################
### code chunk number 56: Models.Stex:1870-1880
###################################################
NumSuccess <- rbinom(5, 30, 0.8)
pSuccess <- NumSuccess/30
resp.mat <- cbind(NumSuccess, 30 - NumSuccess)
pSuccess.glm <- glm(resp.mat ~ 1, binomial)
# glm estimate
plogis(coef( pSuccess.glm ))
# mean proportion Success
mean(pSuccess) 
# mean of logits on response scale
plogis(mean(qlogis( pSuccess )))  



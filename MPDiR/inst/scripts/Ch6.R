### R code from vignette source 'ClassImage/ClassImage.Stex'

###################################################
### code chunk number 1: ClassImage.Stex:9-30
###################################################
options(width=65, digits = 3,  str = strOptions(strict.width = "cut"))
#pkg <- search()[2] 
#while (search()[2] !=  if(.Platform$GUI == "AQUA") "tools:RGUI" else "package:stats") {
##detach(pos = match(pkg, search()))
#spkg <- strsplit( pkg, ":"  )[[1]][2] 
#if (packageHasNamespace(spkg, .libPaths()[1]) )
#unloadNamespace(spkg ) else detach(pos = match(pkg, search()))
# pkg <- search()[2]
 #}
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
library(cacheSweave)
setCacheDir("ClassImage/")


###################################################
### code chunk number 2: ClassImage.Stex:54-89
###################################################
opar <- par(mfrow = c(2, 2))
N <- 32
Stimulus <- matrix(0.5, N, N)
N2 <- N/2
Align <- as.matrix(rbind(expand.grid(N2 - 1:5, (N2-1):(N2-1)), 
				expand.grid(N2 + 1:5, (N2-1):(N2-1))))
Stimulus[Align] <- 0
image(seq(N), seq(N), Stimulus, axes = FALSE, 
	col = grey(c(0.35, 0.5)),
	xlab = "", ylab = "", asp = 1,
	main = "Offset Absent", cex.main = 1.75)

Stimulus <- matrix(0.5, N, N)
NonAlign <-as.matrix(rbind(expand.grid(N2 - 1:5, (N2-1):(N2-1)), 
				expand.grid(N2 + 1:5, N2:N2)))
Stimulus[NonAlign] <- 0
image(seq(N), seq(N), Stimulus, axes = FALSE, 
	col = grey(c(0.35, 0.5)),
	xlab = "", ylab = "", asp = 1,
	main = "Offset Present",  cex.main = 1.75)
	
Stimulus <- matrix(runif(N^2, 0.4, 0.6), N, N)
Stimulus[Align] <- Stimulus[Align] - 0.0975
image(seq(N), seq(N), Stimulus, axes = FALSE,
	xlab = "", ylab = "", asp = 1, 
	col = grey(seq(0.3, 0.7, len = 256)),
	main = "Offset Absent\nin Uniform Noise",  cex.main = 1.75)

Stimulus <- matrix(runif(N^2, 0.4, 0.6), N, N)
Stimulus[NonAlign] <- Stimulus[Align] - 0.0975
image(seq(N), seq(N), Stimulus, axes = FALSE,
	xlab = "", ylab = "", asp = 1, 
	col = grey(seq(0.3, 0.7, len = 256)),
	main = "Offset Present\nin Uniform Noise",  cex.main = 1.75)
par(opar)


###################################################
### code chunk number 3: ClassImage.Stex:183-193
###################################################
N <- 32
w <- 0.65
Contrast <- 0.7
x <- seq(-3, 3, len = N)
Trials <- 10000
Signal <- dnorm(x, sd = w) 
Template <- Signal * cos(2 * pi * x * w)
Stimulus <- matrix(Noise <- rnorm(N * Trials),  N, Trials) + 
	cbind(matrix(rep(Contrast * Signal, Trials/2), 32),
	matrix(0, 32, Trials/2))


###################################################
### code chunk number 4: ClassImage.Stex:213-223
###################################################
opar <- par(mfrow = c(1, 1))
plot(x, rowMeans(Stimulus[, 1:(Trials/2)]), 
	type = "p", ylim = c(-0.25, 0.45), cex = 1.5,
	xlab = "Physical Dimension", ylab = "Average Noise Modulation")
abline(0, 0, col = "grey", lwd = 2)
points(x, rowMeans(Stimulus[, -(1:(Trials/2))]), 
	cex = 1.5, pch = 2)
lines(x, Contrast * Signal, lwd = 2)
lines(x, Contrast * Template, lty = 2, lwd = 2)
par(opar)


###################################################
### code chunk number 5: ClassImage.Stex:248-257 (eval = FALSE)
###################################################
## plot(x, rowMeans(Stimulus[, 1:(Trials/2)]), 
## 	type = "p", ylim = c(-0.25, 0.45), cex = 1.5,
## 	xlab = "Physical Dimension", 
## 	ylab = "Average Stimulus Modulation")
## abline(0, 0, col = "grey", lwd = 2)
## points(x, rowMeans(Stimulus[, -(1:(Trials/2))]), 
## 	cex = 1.5, pch = 2)
## lines(x, Contrast * Signal, lwd = 2)
## lines(x, Contrast * Template, lty = 2, lwd = 2)


###################################################
### code chunk number 6: ClassImage.Stex:262-266
###################################################
DecisionVariable <- drop(Template %*% Stimulus > 0)
Classification <- factor(4 - (rep(c(1, 0), each = Trials/2) + 
	2 * DecisionVariable), levels = 1:4, 
	labels = c("Hit", "FA", "Miss", "CR"))


###################################################
### code chunk number 7: ClassImage.Stex:278-287 (eval = FALSE)
###################################################
## ExpNoise <- matrix(Noise, N, Trials)
## CompImages <- aggregate(t(ExpNoise), 
## 	list(Classification), mean)[, -1]
## matplot(as.matrix(t(CompImages)), type = "l", 
## 	lwd = 2, col = "black", ylim = c(-0.5, 0.5),
## 	xlab = "Physical Stimulus", 
## 	ylab = "Average Noise Modulation")
## legend(2.5, 0.5, legend = c("H", "FA", "M", "CR"), 
## 	lty = 1:4, cex = 1.2, bty = "n")


###################################################
### code chunk number 8: ClassImage.Stex:307-309 (eval = FALSE)
###################################################
## CIWTS <- c(1, 1, -1, -1)
## SimClassIm <-crossprod(CIWTS, as.matrix(CompImages))


###################################################
### code chunk number 9: ClassImage.Stex:318-335
###################################################
opar <- par(mfrow= c(1, 2))
ExpNoise <- matrix(Noise, N, Trials)
CompImages <- aggregate(t(ExpNoise), list(Classification), mean)[, -1]
matplot(as.matrix(t(CompImages)), type = "l", 
	lwd = 2, col = "black", ylim = c(-0.5, 0.5),
	xlab = "Physical Stimulus", ylab = "Average Noise Modulation")
legend(1, 0.5, legend = c("H", "FA", "M", "CR"), lty = 1:4, cex = 1.2,
	bty = "n")
mtext("a", 3, at = 0, adj = 0, cex = 2, line = 1)
SimClassIm <-crossprod(as.matrix(CompImages), c(1, 1, -1, -1))
plot(x, SimClassIm,  cex = 1.75, 
	xlab = "Physical Stimulus",
	ylab = "")
lines(x,   4 * Contrast * Signal, lwd = 2)
lines(x, 4 * Contrast * Template, lty = 2, lwd = 2)
mtext("b", 3, at = -3, adj = 0, cex = 2, line = 1)
par(opar)


###################################################
### code chunk number 10: ClassImage.Stex:390-396
###################################################
idObs <- matrix(0.5, 24, 14)
uppos <- as.matrix(expand.grid(14:19, 8:8))
downpos <- as.matrix(expand.grid(14:19, 7:7))
idObs[uppos] <- 0
idObs[downpos] <- 1
image(idObs, col = grey(c(0, 0.7, 1)), axes = FALSE)


###################################################
### code chunk number 11: ClassImage.Stex:428-430
###################################################
f  <- c(0.7, 1, 0.7)
f %o% f


###################################################
### code chunk number 12: ClassImage.Stex:538-557
###################################################
opar <- par(mfrow = c(1, 2))
x <- seq(1:32) * 20/1000
f <- 0.5 * dnorm(x, 0.32, 0.16) * sin(2 * pi *2 * x/0.64)
plot(c(0, x, 0.66) ,c(0,  f, 0), type = "h", 
		xlab = "Time (s)", 
		ylab = "Relative Luminance Modulation",
		cex.lab = 1.35)
points(x, rep(1.1, length(x)), pch = 16, cex = 1.4,
	col = grey((f - min(f))/diff(range(f))))
mtext("a", 3, 1, at = 0, cex = 2)
f <-  runif(32, -0.3, 0.3) + 0.1 * f
plot(c(0, x, 0.66) ,c(0,  f, 0), type = "h", 
		xlab = "Time (s)", 
		ylab = "",#Relative Luminance Modulation",
		cex.lab = 1.35)
points(x, rep(max(f), length(x)), pch = 16, cex = 1.4,
	col = grey((f - min(f))/diff(range(f))))
mtext("b", 3, 1, at = 0, cex = 2)
par(opar)


###################################################
### code chunk number 13: ClassImage.Stex:590-593
###################################################
data(Gabor)
GaborClsIm.clim <- drop(with(Gabor, tapply(N, list(time, resp), 
		mean)) %*% c(1, 1, -1, -1))


###################################################
### code chunk number 14: ClassImage.Stex:606-607
###################################################
Gabor.lm1 <- lm(N ~ time/resp - 1, Gabor)


###################################################
### code chunk number 15: ClassImage.Stex:622-624
###################################################
GaborClsIm.lm1 <- crossprod(matrix(coef(Gabor.lm1), nrow = 4, byrow = TRUE),
	c(0, 1, -1, -1))


###################################################
### code chunk number 16: ClassImage.Stex:643-654
###################################################
TrueSignal <- with(Gabor, dnorm(unique(Time), mean = 0.32, 
					sd = 0.16) * sin(2 * pi * unique(Time)/0.32))
TrueSignal <- TrueSignal * coef(lm(GaborClsIm.clim ~ TrueSignal - 1))
plot(unique(Gabor$Time), GaborClsIm.clim,
	ylim = c(-0.1, 0.05), xlab = "TIme (s)",
	ylab = "Relative Luminance Modulation",
	cex = 1.5
	)
lines(unique(Gabor$Time), GaborClsIm.lm1, lwd = 2)
lines(unique(Gabor$Time), TrueSignal, col = rgb(0.7, 0.7, 0.7, 0.5), lwd = 5)
abline(0, 0, lty = 2)


###################################################
### code chunk number 17: ClassImage.Stex:668-669 (eval = FALSE)
###################################################
## lm(N ~ time:resp - 1, Gabor)


###################################################
### code chunk number 18: ClassImage.Stex:679-682 (eval = FALSE)
###################################################
## library(gmodels)
## lm(N ~ time:resp - 1, Gabor,
##     contrasts = list(resp = make.contrasts(c(1, 1, -1, -1))))


###################################################
### code chunk number 19: ClassImage.Stex:701-703 (eval = FALSE)
###################################################
## Gabor.lm0 <- update(Gabor.lm1, . ~ . - time:resp)
## anova(Gabor.lm0, Gabor.lm1)


###################################################
### code chunk number 20: ClassImage.Stex:705-707
###################################################
Gabor.lm0 <- update(Gabor.lm1, . ~ . - time:resp)
print(anova(Gabor.lm0, Gabor.lm1), signif.stars = FALSE)


###################################################
### code chunk number 21: ClassImage.Stex:715-716 (eval = FALSE)
###################################################
## drop1(Gabor.lm1, test = "F")


###################################################
### code chunk number 22: ClassImage.Stex:718-719
###################################################
print(drop1(Gabor.lm1, test = "F"), signif.stars = FALSE)


###################################################
### code chunk number 23: ClassImage.Stex:780-785
###################################################
Resp <- factor(unclass(Gabor$resp) < 3, labels = 0:1)[seq(1, nrow(Gabor), 32)]
Stim <- factor(unclass(Gabor$resp)  %%  2, labels = 0:1)[seq(1, nrow(Gabor), 32)]
Gabor.wide <- data.frame(matrix(Gabor$N, ncol = 32, byrow = TRUE))
names(Gabor.wide) <- paste("t", 1:32, sep = "")
Gabor.wide  <- cbind(Resp, Stim,  Gabor.wide)


###################################################
### code chunk number 24: ClassImage.Stex:788-791
###################################################
GaborClsIm.glm1 <- glm(Resp ~ Stim/. - 1,
		family = binomial(link = "probit"),
		data = Gabor.wide)


###################################################
### code chunk number 25: ClassImage.Stex:807-818 (eval = FALSE)
###################################################
## ClsIm <- coef(GaborClsIm.glm1)[seq(4, 66, 2)]
## plot(unique(Gabor$Time), ClsIm, 
## 	xlab = "TIme (s)", pch = 16, col = "black",
## 	ylab = "Classification Image Weights",
## 	cex = 1.5, ylim = c(-2, 1.5)
## 	)
## se <- summary(GaborClsIm.glm1)$coefficients[seq(4, 66, 2), 2]
## polygon(c(unique(Gabor$Time), rev(unique(Gabor$Time))),
## 	c(ClsIm + se, rev(ClsIm - se)), 
## 	col = rgb(0.25, 0.25, 0.25, 0.25))
## abline(0, 0, lty = 3)


###################################################
### code chunk number 26: ClassImage.Stex:831-840
###################################################
TmpStim <- with(Gabor, matrix(N, nrow = 32))	
FTmpStim <- apply(TmpStim, 2, fft)
FTmpStim <- Mod(FTmpStim)
f.mat <- data.frame(t(FTmpStim))
names(f.mat) <- paste("f", 1:32, sep = "")
f.mat <- cbind(Resp, f.mat)
GaborFreqIm.glm <- glm(Resp ~ Stim/. - 1,
		family = binomial(link = "probit"),
		data = f.mat)


###################################################
### code chunk number 27: ClassImage.Stex:844-900
###################################################
opar <- par(mfrow = c(2, 2))
AClsIm <- coef(GaborClsIm.glm1)[seq(3, 66, 2)]
plot(unique(Gabor$Time),AClsIm, 
	xlab = "TIme (s)", pch = 16, col = "black",
	ylab = "Classification Image Weights",
	cex = 1.5, ylim = c(-2, 1.5),
	main = "Signal Absent", cex.lab = 1.5, cex.main = 2
	)
Ase <- summary(GaborClsIm.glm1)$coefficients[seq(3, 66, 2), 2]
polygon(c(unique(Gabor$Time), rev(unique(Gabor$Time))),
	c(AClsIm + Ase, rev(AClsIm - Ase)), border = NA,
	col = rgb(0.25, 0.25, 0.25, 0.25))
abline(0, 0, lty = 3, lwd = 2)
mtext("a", 3, at = 0, adj = 0, cex = 2, line = 1)
PClsIm <- coef(GaborClsIm.glm1)[seq(4, 66, 2)]
plot(unique(Gabor$Time),PClsIm, 
	xlab = "TIme (s)", pch = 16, col = "black",
	ylab = "",
	cex = 1.5, ylim = c(-2, 1.5),
	main = "Signal Present", cex.lab = 1.5, cex.main = 2
	)
Pse <- summary(GaborClsIm.glm1)$coefficients[seq(4, 66, 2), 2]
polygon(c(unique(Gabor$Time), rev(unique(Gabor$Time))),
	c(PClsIm + Pse, rev(PClsIm - Pse)),  border = NA,
	col = rgb(0.25, 0.25, 0.25, 0.25))
abline(0, 0, lty = 3, lwd = 2)
mtext("b", 3, at = 0, adj = 0, cex = 2, line = 1)
AbsFreqIm <- coef(GaborFreqIm.glm)[seq(3, 36, 2)]
TF <- (0:16)/(32 * 0.02)
plot(TF, AbsFreqIm, 
	pch = 16, col = "black",
	xlab = "Temporal Frequency (Hz)",
	ylab = "Classification Image Weights",
	cex = 1.5, ylim = c(-0.2, 1), cex.lab = 1.5
	)
AFse <- summary(GaborFreqIm.glm)$coefficients[seq(3, 36, 2), 2]
polygon(c(TF, rev(unique(TF))),
	c(AbsFreqIm + AFse, rev(AbsFreqIm - AFse)),  border = NA,
	col = rgb(0.25, 0.25, 0.25, 0.25))
abline(0, 0, lty = 3, lwd = 2)
mtext("c", 3, at = 0, adj = 0, cex = 2, line = 1)
PresFreqIm <- coef(GaborFreqIm.glm)[seq(4, 36, 2)]
TF <- (0:16)/(32 * 0.02)
plot(TF, PresFreqIm, 
	pch = 16, col = "black",
	xlab = "Temporal Frequency (Hz)",
	ylab = "",
	cex = 1.5, ylim = c(-0.2, 1), cex.lab = 1.5
	)
PFse <- summary(GaborFreqIm.glm)$coefficients[seq(4, 36, 2), 2]
polygon(c(TF, rev(unique(TF))),
	c(PresFreqIm + PFse, rev(PresFreqIm - PFse)), border = NA,
	col = rgb(0.25, 0.25, 0.25, 0.25))
abline(0, 0, lty = 3, lwd = 2)
mtext("d", 3, at = 0, adj = 0, cex = 2, line = 1)
par(opar)


###################################################
### code chunk number 28: ClassImage.Stex:955-964 (eval = FALSE)
###################################################
## TmpStim <- with(Gabor, matrix(N, nrow = 32))	
## FTmpStim <- apply(TmpStim, 2, fft)
## FTmpStim <- Mod(FTmpStim)
## f.mat <- data.frame(t(FTmpStim))
## names(f.mat) <- paste("f", 1:32, sep = "")
## f.mat <- cbind(Resp, f.mat)
## GaborFreqIm.glm <- glm(Resp ~ Stim/. - 1,
## 		family = binomial(link = "probit"),
## 		data = f.mat)


###################################################
### code chunk number 29: ClassImage.Stex:1061-1068
###################################################
library(mgcv)
Gabor$Stim <- factor(with(Gabor, unclass(resp) %% 2))
Gabor$Resp <- factor(with(Gabor, unclass(resp) < 3), 
	labels = c("0", "1"))
Imat <- model.matrix(~Stim:N - 1, data = Gabor)
Gabor$by1 <- Imat[, 1]
Gabor$by2 <- Imat[, 2]


###################################################
### code chunk number 30: ClassImage.Stex:1080-1083
###################################################
Gab.gam <- gam(Resp ~ s(Time, bs = "ts", k = 25, by = by1) +
	s(Time, bs = "ts", k = 25, by = by2), 
	family = binomial("probit"), data = Gabor)


###################################################
### code chunk number 31: ClassImage.Stex:1126-1127
###################################################
summary(Gab.gam)


###################################################
### code chunk number 32: ClassImage.Stex:1163-1183
###################################################
nd <- data.frame(Time = unique(Gabor$Time),
				N = rep(1, 32),
				by1 = rep(0, 32),
				by2 = rep(1, 32)
				)
opar <- par(mfrow = c(1, 2))
plot(Gab.gam, rug = FALSE, select = 1, shade = TRUE, 
	n = 50, xlab = "Time (s)" , main = "Signal Absent", 
	cex.lab = 1.5, cex.main = 2)
abline(0, 0, lty = 2)
mtext("a", 3, at = 0, adj = 0, cex = 2, line = 1)
plot(Gab.gam, rug = FALSE, select = 2, shade = TRUE, 
	n = 50, xlab = "Time (s)" , main = "Signal Present", 
	cex.lab = 1.5, cex.main = 2)
lines(unique(Gabor$Time), 
	coef(lm(predict(Gab.gam, newdata = nd) ~ TrueSignal - 1)) * TrueSignal, 
	lwd = 3, col = "black", lty = 2)
abline(0, 0, lty = 2)
mtext("b", 3, at = 0, adj = 0, cex = 2, line = 1)
par(opar)


###################################################
### code chunk number 33: ClassImage.Stex:1222-1227
###################################################
Gab.gam2 <- gam(Resp ~ s(Time, bs = "ts", by = by1, k = 25) +
	s(Time, bs = "ts", by = by2, k = 25) +
	s(Time, bs = "ts", by = by1^2, k = 25) +
	s(Time, bs = "ts", by = by2^2, k = 25), 
	family = binomial("probit"), data = Gabor)


###################################################
### code chunk number 34: ClassImage.Stex:1234-1235
###################################################
anova(Gab.gam, Gab.gam2, test = "Chisq")


###################################################
### code chunk number 35: ClassImage.Stex:1254-1269
###################################################
PA <- c("Absent, ", "Present, ")
Ord <- c("1", "2")
Suf <- c("st", "nd")
sh <- c(0, 0, 7.25, -10.5)
opar <- par(mfrow = c(2, 2))
for (ix in 1:4){
        pix <- PA[2 - ix %% 2]
        oix <- Ord[2 - (ix < 3)]
        six <- Suf[2 - (ix < 3)]
        plot(Gab.gam2, rug = FALSE, shade = TRUE, 
                main = bquote(paste("Signal ", .(pix), .(oix)^{.(six)}, " order smooth")),
                select = ix, shift = sh[ix], ylim = c(-1.75, 1.75) - sh[ix], cex.main = 1)
        abline(0, 0, lty = 2)
}
par(opar)


###################################################
### code chunk number 36: ClassImage.Stex:1391-1417 (eval = FALSE)
###################################################
## Trials <- c(100, 500, 1000, 5000, 10000)
## NumSims <- 100
## glc <- glm.control(maxit = 1000)	
## 
## err <- lapply(Trials, function(NumTrials){
##    Noise <- array(rnorm(N * NumTrials * NumSims), 
## 		c(N, NumTrials, NumSims))
##    Stimulus <- Noise + rep(cbind(matrix(rep(Contrast * Signal, 
## 	NumTrials/2), nrow = N),
## 	matrix(0, N, NumTrials/2)), NumSims)
##    resp <- apply(Stimulus, 3, function(x, t) crossprod(x, t), 
## 	t = Template) > 0
##    Stim <- factor(rep(c(1, 0), each = NumTrials/2))
##    sapply(seq(dim(resp)[2]), function(x, y, z) {
## 	Respdf <- data.frame(resp = factor(y[, x], 
## 	   labels = c(0, 1)), 
## 	   N = t(z[, , x]))
## 	 CI.glm <- glm(resp ~ Stim/. - 1, binomial("probit"), 
## 	    Respdf,  control = glc)
## 	 CI <- coef(CI.glm)[seq(4, 2 + 2 * N, 2)] + 
## 	   coef(CI.glm)[seq(3, 2 + 2 * N, 2)]
##    res <- var(Template - coef(lm(Template ~ CI - 1)) * CI)
##    res
## 	   }, y = resp, z = Noise)
## 	})
## names(err) <- Trials


###################################################
### code chunk number 37: ClassImage.Stex:1422-1447
###################################################
load("~/projets/svn/mpdir/Book/examples/ClassIm/GLMvsLM.Rdata")
err.q <- lapply(err, function(x){
				apply(x, 1, quantile, c(0.025, 0.5, 0.975))
	}
)
Trials <- c(100, 500, 1000, 5000, 10000)
load("~/Documents/articles/KMClglm/GAMres.Rdata")
plot(Trials, sapply(err.GAM, median), log = "xy",
	type = "b", ylim = c(5e-7, 5e-2), xlim = c(75, 15000),
	ylab = "Residual Variance", cex = 1.75, lwd = 2,
	cex.lab = 1.375, pch = 24, col = "black", bg = "white")
abline(h=10^seq(-6, -2), col = "grey")
qq <- sapply(err.GAM, quantile, c(0.025, 0.975))
#polygon(c(Trials, rev(Trials)), c(qq[1, ], rev(qq[2, ])), 
#	border = "white", col = rgb(0.5, 0.5, 0.5, 0.25) )
lines(Trials, sapply(err.q, function(x) x[2, 1]), lwd = 2)
lines(Trials, sapply(err.q, function(x) x[2, 2]), lwd = 2)
points(Trials, sapply(err.q, function(x) x[2, 1]), cex = 1.5,
	pch = 21, col = "black", bg = "white")
points(Trials, sapply(err.q, function(x) x[2, 2]), cex = 1.5, 
	pch = 16, col = "black")
points(Trials, sapply(err.GAM, median), col = "black", bg = "white", pch = 24, cex = 1.5)
rect(100, 1e-6, 500, 1e-4, col = "white", border = "white")
legend(100, 1.5e-4, c("LM", "GLM", "GAM"), pch = c(21, 21, 24), 
	col = rep("black", 3), pt.bg = c("black", "white", "white"),  bty = "n", cex = 1.75)


###################################################
### code chunk number 38: ClassImage.Stex:1471-1525
###################################################
opar <- par(mfrow = c(1, 2))
load("~/Documents/articles/KMClglm/KMClglm/GAM500.Rdata")
load("~/projets/svn/mpdir/Book/examples/ClassIm/IntNoise500.Rdata")
errIN.q <- lapply(errIN, function(x){
				apply(x, 1, quantile, c(0.025, 0.5, 0.975))
	}
)
IntNoise <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2)
plot(IntNoise, sapply(err.GAM500, median), log = "xy",
	type = "b", ylim = c(1e-5, 0.1), xlim = c(0.0005, 5),
	ylab = "Residual Variance", cex = 1.75, lwd = 2,
	xlab = expression(paste("Internal Noise (", sigma, ")")),
	cex.lab = 1.375, main = "500 Trials", cex.main = 2,
	pch = 24, col = "black", bg = "white")
abline(h=10^seq(-5, -1), col = "grey")
qq <- sapply(err.GAM500, quantile, c(0.025, 0.975))
#polygon(c(IntNoise, rev(IntNoise)),  c(qq[1, ], rev(qq[2, ])),
#	border = "white", col = rgb(0.5, 0.5, 0.5, 0.25) )
lines(IntNoise, sapply(errIN.q, function(x) x[2, 1]), lwd = 2)
lines(IntNoise, sapply(errIN.q, function(x) x[2, 2]), lwd = 2)
points(IntNoise, sapply(errIN.q, function(x) x[2, 2]), pch = 16, 
	col = "black", cex = 1.75)
points(IntNoise, sapply(errIN.q, function(x) x[2, 1]), pch = 21, 
	col = "black", cex = 1.75, bg = "white")
points(IntNoise, sapply(err.GAM500, median), cex = 1.75,
	pch = 24, col = "black", bg = "white")
mtext("a", 3, 1, at = 0.001, cex = 2)

load("~/Documents/articles/KMClglm/err5000.Rdata")
load("~/projets/svn/mpdir/Book/examples/ClassIm/IntNoise.Rdata")
errIN.q <- lapply(errIN, function(x){
				apply(x, 1, quantile, c(0.025, 0.5, 0.975))
	}
)
plot(IntNoise, sapply(err.GAM5000, median), , log = "xy",
	type = "b", ylim = c(1e-5, 0.1), xlim = c(0.0005, 5),
	ylab = "Residual Variance", cex = 1.75, lwd = 2,
	xlab = expression(paste("Internal Noise (", sigma, ")")),
	cex.lab = 1.375, main = "5000 Trials", cex.main = 2,
	pch = 24, col = "black", bg = "white")
abline(h=10^seq(-5, -1), col = "grey")
qq <- sapply(err.GAM5000, quantile, c(0.025, 0.975))
#polygon(c(IntNoise, rev(IntNoise)),  c(qq[1, ], rev(qq[2, ])),
#	border = "white", col = rgb(0.5, 0.5, 0.5, 0.25) )
lines(IntNoise, sapply(errIN.q, function(x) x[2, 1]), lwd = 2)
lines(IntNoise, sapply(errIN.q, function(x) x[2, 2]), lwd = 2)
points(IntNoise, sapply(errIN.q, function(x) x[2, 2]), pch = 16, 
	col = "black", cex = 1.75)	
points(IntNoise, sapply(errIN.q, function(x) x[2, 1]), pch = 21, 
	col = "black", cex = 1.75, bg = "white")
points(IntNoise, sapply(err.GAM5000, median), cex = 1.75,
	pch = 24, col = "black", bg = "white")
mtext("b", 3, 1, at = 0.001, cex = 2)
par(opar)



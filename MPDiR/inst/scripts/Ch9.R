### R code from vignette source 'MEMods/MEmods.Stex'
### Encoding: UTF-8

###################################################
### code chunk number 0: Code for polmer function  
###################################################
polmer <- function(formula, data, lnk = "logit", which.lme4 = "lme4.0",  ...){
	if(require(which.lme4, character.only = TRUE)) invisible() else stop("need lme4 package installed!\n")
# Get fixed and random terms
   fxForm <- lme4.0:::nobars(formula)
   fxForm[[2]] <- NULL
   rTerms <- lme4.0:::findbars(formula)
   ranTerms <- sapply(rTerms, deparse)
   fixTerms <- labels(terms(fxForm))
	ranNames <- as.vector(sapply(ranTerms, function(x)
		strsplit(x, "\\| ")[[1]][2]
		))
	LeftRanNames <- as.vector(sapply(ranTerms, function(x)
		strsplit(x, " | ", fixed = TRUE)[[1]][1]
		))
	NoRanEff <- !(length(ranTerms) > 0)

		
# Make fixed-effect model matrix
	Resp <- ordered(data[[as.character(formula[[2]])]])
	Cuts <- as.numeric(levels(Resp)[-length(levels(Resp))])
	cumRat <- as.vector(sapply(Cuts, 
		function(x) Resp <= x))
	fRat <- gl(length(Cuts), nrow(data), 
		nrow(data) * length(Cuts))
	X <- model.matrix(~ fRat - 1)
	labs <- sapply(Cuts, function(x)
		paste(x, "|", x+1, sep = "")
		)
	colnames(X) <- labs
	
	fX <- -model.matrix(fxForm, data)[, -1]
	fX.names <- if (inherits(fX, "matrix")) colnames(fX) else paste(fixTerms, levels(data[[fixTerms]])[2], sep = "")
	fX <- kronecker(matrix(rep(1, length(Cuts)), ncol = 1),
		 fX)
	X <- cbind(X, fX)
	colnames(X)[-seq(length(Cuts))] <- fX.names
	p.df <- data.frame(cumRat = cumRat, X = X)
	names(p.df) <- c("cumRat", colnames(X))

# Make random-effect model vectors
	frm <- if(!NoRanEff){
		tmp <- sapply(seq_len(length(ranNames)), 
			function(x)
				rep(data[[ranNames[x]]], length(Cuts)))

		for(ix in seq_len(ncol(tmp))) 
			assign(ranNames[ix], tmp[, ix])
		rxForm <- paste(paste("(", ranTerms, ")", 
			sep = "", collapse = " + "), " - 1")
	as.formula(paste("cumRat ~ .  + ", rxForm))
	} else as.formula("cumRat ~ . - 1")
	for(ix in LeftRanNames)
	  if(ix != "1") 
	  	assign(ix, rep(data[[ix]], each = length(Cuts)))
	  	
	CLL <- if (NoRanEff)
	  substitute(
	  glm(FRM, data = p.df, 
	  	family =  binomial(LNK), ...), 
	  		list(FRM = frm, LNK = lnk))  else
	  substitute(
	  glmer(FRM, data = p.df, 
	  	family =  binomial(LNK), ...),
	  		list(FRM = frm, LNK = lnk))	
	res <- eval(CLL)
#	if (inherits(res, "mer")) print(res, cor = cor) else
#		print(res)
	res
}


###################################################
### code chunk number 1: MEmods.Stex:8-27
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
### code chunk number 2: MEmods.Stex:190-191
###################################################
str(Faces2)


###################################################
### code chunk number 3: MEmods.Stex:211-212
###################################################
coef(summary(glm(Resp ~ Stim, binomial(probit), Faces2)))


###################################################
### code chunk number 4: MEmods.Stex:242-246
###################################################
library(lme4.0)
F2.glmm1 <- glmer(Resp ~ Stim + (1 | Obs) + (1 | Image), 
	data =  Faces2, family = binomial(probit))
F2.glmm1


###################################################
### code chunk number 5: MEmods.Stex:312-314
###################################################
( m <- matrix(c(-0.5, -0.5,
                -1, 1), 2, 2, byrow = TRUE) )


###################################################
### code chunk number 6: MEmods.Stex:318-319
###################################################
solve(m)


###################################################
### code chunk number 7: MEmods.Stex:324-327
###################################################
Stm <- ifelse(Faces2$Stim == "P", 0.5, -0.5)
Stm <- cbind(-1, Stm)
colnames(Stm) <- c("_Criterion", "_d'")


###################################################
### code chunk number 8: MEmods.Stex:333-335
###################################################
 glmer(Resp ~ Stm + (1 | Obs) + (1 | Image) - 1, 
	Faces2,binomial(probit))


###################################################
### code chunk number 9: MEmods.Stex:373-374 (eval = FALSE)
###################################################
##  (Stim | Obs)


###################################################
### code chunk number 10: MEmods.Stex:383-389
###################################################
F2.glmm2 <- glmer(Resp ~ Stim + (Stim | Obs) + (1 | Image), 
	Faces2, binomial(probit))
F2.glmm3 <- glmer(Resp ~ Stim + (1 | Obs) + (Stim | Image), 
	Faces2, binomial(probit))
F2.glmm4 <- glmer(Resp ~ Stim + (Stim | Obs) + (Stim | Image), 
	Faces2, binomial(probit))


###################################################
### code chunk number 11: MEmods.Stex:393-395
###################################################
anova(F2.glmm1, F2.glmm2, F2.glmm4)
anova(F2.glmm1, F2.glmm3, F2.glmm4)


###################################################
### code chunk number 12: MEmods.Stex:440-443
###################################################
library(ordinal)
clmm(factor(SimRating) ~ sibs + (1 | Obs), data = Faces, 
	link = "probit", Hess = TRUE)


###################################################
### code chunk number 13: MEmods.Stex:494-511
###################################################
SimRating <- ordered(Faces$SimRating)
Cuts <- seq(0, length(levels(SimRating)) - 2)
cumRat <- as.vector(sapply(Cuts, 
	function(x) Faces$SimRating <= x))
fRat <- gl(10, nrow(Faces), nrow(Faces) * length(Cuts))
X <- model.matrix(~ fRat - 1)
X <- cbind(X, -(Faces$sibs == "1"))
ColNames <- c("Resp", sapply(Cuts, function(x) 
 	paste(x, '|', x + 1, sep = "")),  "sibs")
Rat.df <- data.frame(Resp = cumRat, X = X)
names(Rat.df) <- ColNames

Obs <- rep(Faces$Obs, length(Cuts))
Image <- rep(Faces$Image, length(Cuts)) 

g1 <- glmer(Resp ~ . - 1 + (sibs - 1 | Obs), 
	Rat.df, binomial("probit"))


###################################################
### code chunk number 14: MEmods.Stex:532-535
###################################################
Faces.glmm <- polmer(SimRating ~ sibs + (sibs1 - 1 | Obs) +
	(sibs1 - 1| Image), data = Faces, lnk = "probit")
print(Faces.glmm, cor = FALSE)


###################################################
### code chunk number 15: MEmods.Stex:619-625
###################################################
NLevs <- 7
Intensity <- 10^seq(-2, -0.5, len = NLevs)
Ntrials <- 50
NRun <- 10
mu <- 0.1
sigma <- 0.05


###################################################
### code chunk number 16: MEmods.Stex:631-637
###################################################
beta0 <- -mu/sigma
beta1 <- 1/sigma
LinPred <- function(Intensity, mu, sigma)
	-mu/sigma + (1/sigma) * Intensity
xx <- 10^seq(-2, -0.5, len = 200)
TrueFunc <- pnorm(LinPred(xx, mu, sigma))


###################################################
### code chunk number 17: MEmods.Stex:641-642
###################################################
pnorm(-beta0, lower.tail = FALSE)


###################################################
### code chunk number 18: MEmods.Stex:649-652
###################################################
DV0 <- rep(beta0, NLevs * NRun * Ntrials) +  
  rep(beta1, NLevs * NRun * Ntrials) *
  rep(Intensity, NRun * Ntrials)


###################################################
### code chunk number 19: MEmods.Stex:678-689
###################################################
SimResp0 <- rbinom(NLevs * NRun * Ntrials, 1, 
	pnorm(DV0))
Sim0.df <- data.frame(Intensity = Intensity, 
	Resp = SimResp0, 
	Run = rep(paste("R", seq(NRun), sep = ""), 
		each = NLevs * Ntrials)) 
Prop0.mat <- with(Sim0.df, 
	table(Resp, round(Intensity, 3), Run))[2, , ]/Ntrials
Prop0.df <- stack(as.data.frame(Prop0.mat))
Prop0.df$Intensity <- Intensity
names(Prop0.df) <- c("PCorrect", "Run", "Intensity")


###################################################
### code chunk number 20: MEmods.Stex:703-704
###################################################
 set.seed(19051960)


###################################################
### code chunk number 21: MEmods.Stex:706-712
###################################################
SD1 <- 0.5
b0 <- rnorm(NRun, sd = SD1)
DV1 <- rep(beta0, NLevs * NRun * Ntrials) + 
  rep(b0, each = NLevs * Ntrials) + 
  rep(beta1, NLevs * NRun * Ntrials) *
  rep(Intensity, NRun * Ntrials)


###################################################
### code chunk number 22: MEmods.Stex:723-724
###################################################
pnorm(-beta0 + c(-SD1, SD1), lower.tail = FALSE)


###################################################
### code chunk number 23: MEmods.Stex:736-741
###################################################
SimResp1 <- rbinom(NLevs * NRun * Ntrials, 1, pnorm(DV1))
Sim1.df <- data.frame(Intensity = Intensity, 
	Resp = SimResp1, 
	Run = rep(paste("R", seq(NRun), sep = ""), 
		each = NLevs * Ntrials)) 


###################################################
### code chunk number 24: MEmods.Stex:748-753
###################################################
Prop1.mat <- with(Sim1.df, 
	table(Resp, round(Intensity, 3), Run))[2,, ]/Ntrials
Prop1.df <- stack(as.data.frame(Prop1.mat))
Prop1.df$Intensity <- Intensity
names(Prop1.df) <- c("PCorrect", "Run", "Intensity")


###################################################
### code chunk number 25: MEmods.Stex:774-780
###################################################
SD2 <- 10
b1 <- rnorm(NRun, sd = SD2)
DV2 <- rep(beta0, NLevs * NRun * Ntrials) + 
  ( rep(beta1, NLevs * NRun * Ntrials) +
    rep(b1, each = NLevs * Ntrials) ) *
  rep(Intensity, NRun * Ntrials)


###################################################
### code chunk number 26: MEmods.Stex:786-796
###################################################
SimResp2 <- rbinom(NLevs * NRun * Ntrials, 1, pnorm(DV2))
Sim2.df <- data.frame(Intensity = Intensity, 
	Resp = SimResp2, 
	Run = rep(paste("R", seq(NRun), sep = ""), 
		each = NLevs * Ntrials)) 
Prop2.mat <- with(Sim2.df, 
	table(Resp, round(Intensity, 3), Run))[2,, ]/Ntrials
Prop2.df <- stack(as.data.frame(Prop2.mat))
Prop2.df$Intensity <- Intensity
names(Prop2.df) <- c("PCorrect", "Run", "Intensity")


###################################################
### code chunk number 27: MEmods.Stex:809-840
###################################################
opar <- par(mfrow = c(1, 3))	
plot(xx, TrueFunc, type = "l", lwd = 3, log = "x",
	xlab = "Intensity", ylab = "Pr(Correct)")
for (ix in levels(Prop0.df$Run)){
  lines(PCorrect ~ Intensity, Prop0.df, subset = Run == ix, col = "grey")
  }
points(Intensity, t(with(Sim0.df, 
	table(Resp, round(Intensity, 3))))[, 2] / 
		(Ntrials * NRun), pch = 21, bg = "white", 
		cex = 1.2)	
plot(xx, TrueFunc, type = "l", lwd = 3, log = "x",
	xlab = "Intensity", ylab = "Pr(Correct)",
	main = "Random Intercept")
for (ix in levels(Prop1.df$Run)){
  lines(PCorrect ~ Intensity, Prop1.df, subset = Run == ix, col = "grey")
  }
points(Intensity, t(with(Sim1.df, 
	table(Resp, round(Intensity, 3))))[, 2] / 
		(Ntrials * NRun), pch = 21, bg = "white", 
		cex = 1.2)
plot(xx, TrueFunc, type = "l", lwd = 3, log = "x",
	xlab = "Intensity", ylab = "Pr(Correct)",
	main = "Random Slope")
for (ix in levels(Prop2.df$Run)){
  lines(PCorrect ~ Intensity, Prop2.df, subset = Run == ix, col = "grey")
  }
points(Intensity, t(with(Sim2.df, 
	table(Resp, round(Intensity, 3))))[, 2] / 
		(Ntrials * NRun), pch = 21, bg = "white", 
		cex = 1.2)
par(opar)


###################################################
### code chunk number 28: MEmods.Stex:860-866
###################################################
Sim0.lst <- lmList(Resp ~ Intensity | Run, 
	Sim0.df, binomial(probit))
Sim1.lst <- lmList(Resp ~ Intensity | Run, 
	Sim1.df, binomial(probit))
Sim2.lst <- lmList(Resp ~ Intensity | Run, 
	Sim2.df, binomial(probit))


###################################################
### code chunk number 29: MEmods.Stex:885-897
###################################################
print(  plot(confint(Sim0.lst)), 
	more = TRUE,
	split = c(1, 1, 3, 1)
	)
print(  plot(confint(Sim1.lst), main = "Random Intercept"), 
	more = TRUE,
	split = c(2, 1, 3, 1)
	)
print(  plot(confint(Sim2.lst), main = "Random Slope"), 
	more = FALSE,
	split = c(3, 1, 3, 1)
	)


###################################################
### code chunk number 30: MEmods.Stex:938-947
###################################################
Sim0.glmm <- glmer(Resp ~ Intensity + (1 | Run) +
	(Intensity - 1| Run), 
	data = Sim0.df,  family = binomial(probit))
Sim1.glmm <- glmer(Resp ~ Intensity + (1 | Run) +
	(Intensity - 1| Run), 
	data = Sim1.df,  family = binomial(probit))
Sim2.glmm <- glmer(Resp ~ Intensity + (1 | Run) +
	(Intensity - 1| Run), 
	data = Sim2.df,  family = binomial(probit))


###################################################
### code chunk number 31: MEmods.Stex:954-955 (eval = FALSE)
###################################################
## (Intensity | Run)


###################################################
### code chunk number 32: MEmods.Stex:963-964
###################################################
summary(Sim0.glmm)


###################################################
### code chunk number 33: MEmods.Stex:989-996
###################################################
sapply(Sim.lst <- list("Sim0" = Sim0.glmm, 
	"Sim1" = Sim1.glmm, "Sim2" = Sim2.glmm), 
	function(x){
		VC <- VarCorr(x)
		sapply(VC, attr, which = "stddev")
	}
 )


###################################################
### code chunk number 34: MEmods.Stex:1017-1018
###################################################
sapply(Sim.lst, fixef)


###################################################
### code chunk number 35: MEmods.Stex:1028-1033
###################################################
lapply(list(Sim0 = Sim0.df, 
   Sim1 = Sim1.df, Sim2 = Sim2.df), 
   function(x) summary( 
      glm(Resp ~ Intensity, binomial(probit), x) )$coef
	)


###################################################
### code chunk number 36: MEmods.Stex:1062-1063
###################################################
str(Context)


###################################################
### code chunk number 37: MEmods.Stex:1073-1111
###################################################
opar <- par(mfrow = c(1, 2))
plot(c(0.1, 0.9), c(0.1, 0.9), type = "n", axes = FALSE,
	xlab = "", ylab = "")
rect(0.19, 0.5 - 1/20, 0.21, 0.5 + 1/20, col = "white",
	lty = 2, lwd = 2)
for (ix in c(1, 3, 5, 9, 11, 13) + 0.5)
   rect(0.39, (ix - 0.75)/15, 0.41, (ix + 0.75)/15,
   		col = "lightgrey", border = NA)
rect(0.39, 0.5 - 1/20, 0.41, 0.5 + 1/20, col = "white",
	lty = 2, lwd = 2)
for (ix in c(1, 3, 5, 9, 11, 13) + 0.5)
   rect(0.59, (ix - 0.75)/15, 0.61, (ix + 0.75)/15,
   		col = "darkgrey", border = NA)
rect(0.59, 0.5 - 1/20, 0.61, 0.5 + 1/20, col = "white",
	lty = 2, lwd = 2)
for (ix in c(1, 3, 5, 9, 11, 13) + 0.5)
   rect(0.79, (ix - 0.75)/15, 0.81, (ix + 0.75)/15,
   		col = "black", border = NA)
rect(0.79, 0.5 - 1/20, 0.81, 0.5 + 1/20, col = "white",
	lty = 2, lwd = 2)
Context <- within(Context, Pc <- NumYes/(NumYes + NumNo))
Pc.sd <- with(Context, tapply(Pc, list(TargCntr, ContCntr), sd))
Pc.sem <- Pc.sd/sqrt(6)
Pc.mean <- with(Context, tapply(Pc, list(TargCntr, ContCntr), mean))
with(Context, 
  interaction.plot(TargCntr, ContCntr, 
	Pc, mean, pch = 1:4, bg = "white", cex = 1,
	type = "b", ylim = c(0, 1), 
	fixed = TRUE, 
	trace.label = "\n   Flanker \n   Contrast",
	xlab = expression(paste("Target Contrast ", C[t])),
	ylab = expression(paste("Yes rate P(yes | ", C[t], ")")), 
	)
)
for(ix in 1:4)
	segments(1:5, (Pc.mean + Pc.sem)[, ix], 1:5,
		(Pc.mean - Pc.sem[, ix])[, ix])
par(opar)


###################################################
### code chunk number 38: MEmods.Stex:1133-1150 (eval = FALSE)
###################################################
## Context <- within(Context, Pc <- NumYes/(NumYes + NumNo))
## Pc.sd <- with(Context, tapply(Pc, list(TargCntr, ContCntr), sd))
## Pc.sem <- Pc.sd/sqrt(6)
## Pc.mean <- with(Context, tapply(Pc, list(TargCntr, ContCntr), mean))
## with(Context, 
##   interaction.plot(TargCntr, ContCntr, 
## 	Pc, mean, pch = 1:4, bg = "white", cex = 1,
## 	type = "b", ylim = c(0, 1), 
## 	fixed = TRUE, 
## 	trace.label = "\n   Flanker \n   Contrast",
## 	xlab = expression(paste("Target Contrast ", C[t])),
## 	ylab = expression(paste("Yes rate P(yes | ", C[t], ")")), 
## 	)
## )
## for(ix in 1:4)
## 	segments(1:5, (Pc.mean + Pc.sem)[, ix], 1:5,
## 		(Pc.mean - Pc.sem[, ix])[, ix])


###################################################
### code chunk number 39: MEmods.Stex:1170-1181
###################################################
print(
xyplot(Pc ~ TargCntr | Obs, data = Context, 
	groups = ContCntr,	type = c("l", "p"),
	xlab = expression(paste("Target Contrast ", C[t])),
	ylab = expression(paste("Yes rate P(yes | ", C[t], ")")),
	auto.key =	list(space = "right", 
		title = "Flanker \nContrast ", cex = 0.75), 
	par.settings = 
		list(strip.background = list(col = "lightgrey"), 
		superpose.symbol = list(col = "black", pch = 1:4),
		superpose.line = list(col = "black", lty = 4:1)))  )		


###################################################
### code chunk number 40: MEmods.Stex:1223-1233
###################################################
Context <- within(Context,  {
	TrgCnt1000 <- 1000 * TargCntr
	FCtxCnt <- factor(ContCntr)
	resp <- cbind(NumYes, NumNo)
	} )
Context.lst <- lmList(resp ~ TrgCnt1000 * FCtxCnt | Obs, 
	Context, binomial(probit))
ord <- order(coef(Context.lst)[[1]])
print(
  plot(confint(Context.lst), order = ord) )


###################################################
### code chunk number 41: MEmods.Stex:1295-1298
###################################################
Context.glmm1 <- glmer(cbind(NumYes, NumNo) ~ 
	FCtxCnt * TrgCnt1000 + (TrgCnt1000 | Obs/FCtxCnt),
	data = Context, family = binomial("probit"))


###################################################
### code chunk number 42: MEmods.Stex:1320-1337
###################################################
getREmat <- function(obj, ...){
	CV <- VarCorr(obj)
	d <- do.call(rbind, 
		lapply(names(CV), function(nm, m){
			V <- diag(m[[nm]])
			SD <- attr(m[[nm]], "stddev")
			Cor <- attr(m[[nm]], "correlation")[2]
			Corc <- as.character(round(Cor, 3))
			if (Cor < 0) Corc <- paste(" ", Corc, sep = "")
			d <- data.frame(Groups = c(nm, ""), Name = names(V), Variance = V,
				Std.Dev = SD, Corr = c(NA, Corc))
		}, m = CV) )
	print(d, na.print = "", row.names = FALSE, ...)
	d$Corr <- as.numeric(d$Corr)
	invisible(d)
}	
getREmat(Context.glmm1, digits = 3)


###################################################
### code chunk number 43: MEmods.Stex:1350-1360
###################################################
Context.glmm2a <- glmer(cbind(NumYes, NumNo) ~ 
	FCtxCnt * TrgCnt1000 + (TrgCnt1000 | Obs) + 
	(1 | Obs:FCtxCnt) + (TrgCnt1000 + 0 | Obs:FCtxCnt),
	data = Context, family = binomial("probit"))
Context.glmm2b <- glmer(cbind(NumYes, NumNo) ~
	FCtxCnt * TrgCnt1000 + (1 | Obs) + 
	(TrgCnt1000 + 0 | Obs) + 
	(TrgCnt1000 | Obs:FCtxCnt),
	data = Context, family = binomial("probit"))
anova(Context.glmm1, Context.glmm2a)


###################################################
### code chunk number 44: MEmods.Stex:1362-1363
###################################################
anova(Context.glmm1, Context.glmm2b)


###################################################
### code chunk number 45: MEmods.Stex:1414-1424
###################################################
Context.glmm3a <- glmer(cbind(NumYes, NumNo) ~ 
	FCtxCnt * TrgCnt1000 + 
	(1 | Obs) + (TrgCnt1000 | Obs:FCtxCnt),
	data = Context, family = binomial("probit"))
Context.glmm3b <- glmer(cbind(NumYes, NumNo) ~ 
	FCtxCnt * TrgCnt1000 + 
	(TrgCnt1000 + 0 | Obs) + (TrgCnt1000 | Obs:FCtxCnt),
	data = Context, family = binomial("probit"))
anova(Context.glmm2b, Context.glmm3a)
anova(Context.glmm2b, Context.glmm3b)


###################################################
### code chunk number 46: MEmods.Stex:1469-1474
###################################################
Context.glmm4 <-  glmer(cbind(NumYes, NumNo) ~
	 FCtxCnt + TrgCnt1000 + 
	(TrgCnt1000 + 0 | Obs) + (TrgCnt1000 | Obs:FCtxCnt),
	data = Context, family = binomial("probit"))
anova(Context.glmm3b, Context.glmm4)


###################################################
### code chunk number 47: MEmods.Stex:1544-1561
###################################################
xx <- seq(0, 0.01, len = 200)
gm4.fe <- fixef(Context.glmm4)
intc <- c(gm4.fe[1], gm4.fe[1] + gm4.fe[-c(1, 5)])
plot(c(0, 0.011), c(0, 1), type = "n",
	xlab = expression(paste("Target Contrast ", C[t])),
	ylab = expression(paste("Yes rate P(yes | ", C[t], ")"))
	)
for (ix in 1:4) {
  points(unique(Context$TargCntr), Pc.mean[, ix], 
  	ylim = c(0, 0.01), cex = 1.4,
  	pch = ix)
  lines(xx, pnorm(intc[ix] + gm4.fe[5] * xx * 1000), 
  	lty = 5 - ix, lwd = 3)
 }
legend(0.007, 0.45, c(0, 0.01, 0.05, 0.4), pch = 21:24, 
	lty = 1:4, bty = "n",
	 title = "Flanker of \n Contrast ")	


###################################################
### code chunk number 48: MEmods.Stex:1586-1604
###################################################
gm4.re <- ranef(Context.glmm4)
opar <- par(mfrow = c(2, 3))
for (obs in levels(Context$Obs)[c(4:6, 1:3)]){
	NObs <- which(LETTERS %in% obs) - 1 
	plot(Pc ~ TargCntr, Context, 
		pch = rep(1:4, each = 5),
		xlim = c(0, 0.01), ylim = c(0, 1),
		subset = Obs == obs, cex = 1.4,
		xlab = "Target Contrast",
		ylab = "Yes rate",
		main = paste("Obs: ",obs))
	for (Cond in 1:4) {
	  nn <- 4 * NObs + Cond
	  lines(xx, pnorm(intc[Cond] + gm4.re[[1]][nn, 1] +
	  	 (gm4.fe[5] +  gm4.re[[1]][nn, 2] +
	  	 	gm4.re[[2]][obs, 1]) * xx * 1000), 
		lty = 5 - Cond, lwd = 3) } }
par(opar)


###################################################
### code chunk number 49: MEmods.Stex:1688-1707
###################################################
library(MLDS)
kk.lst <- list(R1 = kk1, R2 = kk2, R3 = kk3)
kk.mlds.lst <- lapply(kk.lst, mlds)
opar <- par(mfrow = c(1, 2))
plot(0:1, c(0, 10), type = "n", 
	xlab = "r", 
	ylab = "Difference Scale")
for(ix in 1:3) lines(kk.mlds.lst[[ix]], type = "b")
lines(kk.mlds.lst[[1]]$stimulus, 
	rowMeans(sapply(kk.mlds.lst, "[[", "pscale")),
	lwd = 3)
lines(mlds(do.call(rbind, kk.lst)), lty = 2, lwd = 2)
plot(0:1, c(0, 1), type = "n", 
	xlab = "r", 
	ylab = "Difference Scale")
for(ix in 1:3) lines(kk.mlds.lst[[ix]], standard.scale = TRUE)
x <- seq(0, 0.98, len = 100)
lines(x, (x/0.98)^2, lwd = 3)
par(opar)


###################################################
### code chunk number 50: MEmods.Stex:1787-1795
###################################################
Stim <- attr(kk1, "stimulus")
kk <- do.call(rbind, kk.lst)
kk$Run <- factor(rep(names(kk.lst), each = nrow(kk1)))
kk$dv <- matrix(Stim[unlist(kk[, -c(1, 6)])]^2, 
	ncol = 4) %*% c(1, -1, -1, 1)
kk.glmm <- glmer(resp ~ dv + (dv + 0 | Run) - 1, kk,
	binomial("probit"))
kk.glmm


###################################################
### code chunk number 51: MEmods.Stex:1856-1858
###################################################
data(Transparency)
str(Transparency)


###################################################
### code chunk number 52: MEmods.Stex:1866-1877
###################################################
Stim <- attr(Transparency, "stimulus")
Tr.mlds.lst <- lapply(levels(Transparency$Obs), function(obs) {
	tmp.df <- subset(Transparency, Obs == obs, select = 1:5)
	mlds(as.mlds.df(tmp.df, st = Stim))
	})
Tr.mlds.cb <- mlds(as.mlds.df(Transparency[, -6], st = Stim))
plot(Stim, rowMeans(sapply(Tr.mlds.lst, "[[", "pscale")), 
	type = "l", lwd = 3,  ylim = c(0, 20),
	xlab = "Index of Refraction", ylab = "Difference Scale")
for(md in Tr.mlds.lst) lines(md, type = "b")
lines(Tr.mlds.cb, lty = 2, lwd = 3)


###################################################
### code chunk number 53: MEmods.Stex:1941-1955
###################################################
stanScale <- sapply(Tr.mlds.lst, "[[", "pscale")
stanScale <- apply(stanScale, 2, function(x) x/x[length(x)])
colnames(stanScale) <- levels(Transparency$Obs)
Tr.DV.lst <- lapply(levels(Transparency$Obs), function(obs) {
     tmp.df <- subset(Transparency, Obs == obs, select = 1:5)
     RX <- make.ix.mat(tmp.df)
     MX <- RX[, -1]
     MX <- cbind(stim.1 = ifelse(rowSums(MX), 1, 0), MX) 
     as.matrix(MX) %*% stanScale[, obs]
     })	
Transparency$DV <- do.call(rbind, Tr.DV.lst)
Tr.glmm <- glmer(resp ~ DV + (DV + 0  | Obs) - 1, 
     Transparency, binomial("probit"))	
Tr.glmm


###################################################
### code chunk number 54: MEmods.Stex:1973-1982
###################################################
re <- ranef(Tr.glmm)$Obs$DV
fe <- fixef(Tr.glmm)
pref.df <- data.frame(Stimulus = rep(Stim, 6),
	FEpsc = fe * c(stanScale),
	REpsc = (rep(re, each = 10) + fe) * c(stanScale),
	pscale = as.vector(sapply(Tr.mlds.lst, function(x) 
		x$pscale)),
	Obs = rep(paste("O", 1:6, sep = ""), each = 10)
	)


###################################################
### code chunk number 55: MEmods.Stex:1992-2004
###################################################
print(
xyplot(pscale ~ Stimulus | Obs, pref.df,
	xlab = "Index of Refraction", ylab = "Difference Scale",
	subscripts = TRUE, ID = pref.df$Obs,
	panel = function(x, y, subscripts, groups, ...){
		panel.xyplot(x, y, col = "black")
		panel.lines(x, pref.df$REpsc[subscripts], 
			col = "black", lwd = 2)
		panel.lines(x, pref.df$FEpsc[subscripts], 
			lty = 2, lwd = 2, col = "black")
		} )
)


###################################################
### code chunk number 56: MEmods.Stex:2055-2061
###################################################
cStim <- mean(Stim[-1])
Tr.df <- data.frame(
    Stim = Stim[-1] - cStim,
    Obs = rep(factor(paste("O", 1:6, sep = "")), each = 9),
    logDS = c(log10(sapply(Tr.mlds.lst, coef)))
    )


###################################################
### code chunk number 57: MEmods.Stex:2084-2090
###################################################
plot(c(0, 20), c(0, 2.5), type = "n",
	xlab = "Difference Scale",
	ylab = "Standard Error")
for(ix in 1:6)
	points(coef(summary(Tr.mlds.lst[[ix]]$obj))[, 1:2],
		pch =as.character(ix))


###################################################
### code chunk number 58: MEmods.Stex:2118-2128
###################################################
P3 <- lmer(logDS ~ poly(Stim, degree = 6) + 
	(Stim + I(Stim^2) + I(Stim^3) | Obs), Tr.df, 
	REML = FALSE)
P2 <- lmer(logDS ~ poly(Stim, degree = 6) + 
	(Stim + I(Stim^2) | Obs), Tr.df, REML = FALSE)
P1 <- lmer(logDS ~ poly(Stim, degree = 6) + 
	(Stim | Obs), Tr.df, REML = FALSE)
P0 <- lmer(logDS ~ poly(Stim, degree = 6) + 
	(1 | Obs), Tr.df, REML = FALSE)
anova(P0, P1, P2, P3)


###################################################
### code chunk number 59: MEmods.Stex:2144-2152
###################################################
P2.a <- lmer(logDS ~ poly(Stim, degree = 6) + 
	(Stim | Obs) + (I(Stim^2) + 0 | Obs), Tr.df, 
	REML = FALSE)
P2.b <-  lmer(logDS ~ poly(Stim, degree = 6) + 
	(1 | Obs) + 
	(Stim + 0 | Obs) + (I(Stim^2) + 0 | Obs), Tr.df, 
	REML = FALSE)
anova(P2, P2.a, P2.b)


###################################################
### code chunk number 60: MEmods.Stex:2205-2214
###################################################
P2.2 <- lmer(logDS ~ poly(Stim, degree = 2) + 
	(Stim + I(Stim^2) | Obs), Tr.df)
P2.3 <- lmer(logDS ~ poly(Stim, degree = 3) + 
	(Stim + I(Stim^2) | Obs), Tr.df)
P2.4 <- lmer(logDS ~ poly(Stim, degree = 4) + 
	(Stim + I(Stim^2) | Obs), Tr.df)
P2.5 <- lmer(logDS ~ poly(Stim, degree = 5) + 
	(Stim + I(Stim^2) | Obs), Tr.df)
anova(P2, P2.5, P2.4, P2.3, P2.2)


###################################################
### code chunk number 61: MEmods.Stex:2230-2233
###################################################
feP23 <- fixef(P2.3)
fePred23 <- cbind(1, poly(Tr.df$Stim, degree = 3)) %*% feP23
Tr.df$fit23 <- fitted(P2.3)	


###################################################
### code chunk number 62: MEmods.Stex:2245-2263
###################################################
print(
xyplot(10^logDS ~ I(Stim + cStim) | Obs, Tr.df, 
  ylim = c(0, 20), xlab = "Index of Refraction",
  ylab = "Difference Scale", subscripts = TRUE, ID = Tr.df$Obs,
  panel = function(x, y, subscripts, groups, ...){
   panel.xyplot(x, y, col = "black")
   panel.lines(x, 10^(Tr.df$fit23[subscripts]), col = "black", lwd = 2)
   panel.lines(x, 10^fePred23[subscripts], lty = 2, lwd = 2, col = "black")
   }), more = TRUE, split = c(1, 1, 2, 1) )
print(
xyplot(logDS ~ I(Stim + cStim) | Obs, Tr.df,
  xlab = "Index of Refraction", ylab = "Log Difference Scale",
  subscripts = TRUE, ID = Tr.df$Obs,
  panel = function(x, y, subscripts, groups, ...){
   panel.xyplot(x, y, col = "black")
   panel.lines(x, Tr.df$fit23[subscripts], col = "black", lwd = 2)
   panel.lines(x, fePred23[subscripts], lty = 2, lwd = 2, col = "black")
   }), more = FALSE, split = c(2, 1, 2, 1))



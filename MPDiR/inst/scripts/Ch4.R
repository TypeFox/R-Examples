### R code from vignette source 'PsychFun/PsychFun.Stex'
### Encoding: UTF-8

###################################################
### code chunk number 1: PsychFun.Stex:9-30
###################################################
options(width=65, digits = 3, str = strOptions(strict.width = "cut"))
pkg <- search()[2] 
while (search()[2] !=  "package:cacheSweave") {
#if(.Platform$GUI == "AQUA") "tools:RGUI" else "package:stats") {
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
#library(cacheSweave)
ltheme <- canonical.theme("pdf", color = FALSE) ## in-built B&W theme 
ltheme$strip.background$bg <- "grey85" ## change strip bg 
lattice.options(default.theme = ltheme) ## set as default 
formals(deparse)$width.cutoff <- 50L
assignInNamespace("deparse", deparse, "base")
rm(deparse)


###################################################
### code chunk number 2: PsychFun.Stex:118-128
###################################################
n <- 200
M <- 7
Intensity <- seq(-3, 3, len = M)
p <- pnorm(Intensity)
r <- rbinom(M, n, p)
pObs <- r/n
plot(Intensity, pObs, cex = 1.75, 
	cex.lab = 1.37, cex.axis = 1.2,
	xlab = "Stimulus Level",
	ylab = "P[Yes]")


###################################################
### code chunk number 3: PsychFun.Stex:138-147 (eval = FALSE)
###################################################
## n <- 200
## M <- 7
## Intensity <- seq(-3, 3, len = M) # make a  vector of levels
## p <- pnorm(lntensity)  
## nYes <- rbinom(M, n, p)
## pObs <- nYes/n
## plot(Intensity, pObs, 
## 	xlab = "Stimulus Level",
## 	ylab = "P[Yes]")


###################################################
### code chunk number 4: PsychFun.Stex:268-285
###################################################
opar <- par(mfrow = c(1, 2))  # set up a 1x2 array of plots

# plot three location-scale psychometric function families
Intensity <- seq(-3, 3, len = 100)
plot(Intensity, pnorm(Intensity), type = "l", lwd = 2,
	ylab = "P[Intensity]") 
lines(Intensity,plogis(Intensity, scale = sqrt(3)/pi), lty = 2, lwd = 3)
lines(Intensity, pcauchy(Intensity, scale = 5/8), lty = 3, lwd = 3)
mtext("a", adj = 0, cex = 2, line = 1)

# plot two other psychometric function families (not location-scale families)
Intensity2 <- seq(0.1, 10, len = 200)
plot(Intensity2, plnorm(Intensity2),  type = "l", lwd = 2, log = "x",
	xlab = "Intensity", ylab = "P[Intensity]")
lines(Intensity2, pweibull(Intensity2, shape = 1, scale = 1.5), lty = 2, lwd = 3)
mtext("b", adj = 0, cex = 2, line = 1)
par(opar)


###################################################
### code chunk number 5: PsychFun.Stex:447-454
###################################################
library(MPDiR)
data(HSP)
HSP$nyes <- round(with(HSP, N * p/100))
HSP$nno <- with(HSP, N - nyes)
HSP <- transform(HSP, id = Obs:Run)
HSP$id <- factor(HSP$id[, drop = TRUE], levels = unique(HSP$id))
HSP$Obs <- factor(HSP$Obs, levels = c("SH", "SS", "MHP"))


###################################################
### code chunk number 6: PsychFun.Stex:492-499
###################################################
lnorm <- function(p, d) {
   mu <- p[1]
   sigma <- p[2]
   I <- log(d$Q);
   pr <- pnorm(I, mean = mu, sd = sigma)  # Gaussian cdf
   -sum( d$nyes * log(pr) +  d$nno * log( 1-pr ))
   }


###################################################
### code chunk number 7: PsychFun.Stex:510-520
###################################################
lnorm <- function(p, d) {
   mu <- p[1]
   sigma <- p[2]
   I <- log(d$Q);
   logpyes <- pnorm(I, mean = mu, sd = sigma, 
      lower.tail = TRUE, log.p = TRUE)
   logpno  <- pnorm(I, mean = mu, sd = sigma, 
      lower.tail = FALSE, log.p = TRUE)
   -sum( d$nyes * logpyes +  d$nno * logpno)
	}


###################################################
### code chunk number 8: PsychFun.Stex:548-551
###################################################
HSP$nyes <- round(with(HSP, N * p/100))
HSP$nno <- with(HSP, N - nyes)
SHR2 <- subset(HSP, Obs == "SH" & Run == "R2")


###################################################
### code chunk number 9: PsychFun.Stex:556-557
###################################################
SHR2.norm <- optim(par = c(5, 0.35), lnorm, d = SHR2)


###################################################
### code chunk number 10: PsychFun.Stex:597-601 (eval = FALSE)
###################################################
## qq <- seq(3, 6, len = 200)	
## pred.nrm <- pnorm(qq, mean = SHR2.norm$par[1], 
##     sd = SHR2.norm$par[2])
## lines(qq, pred.nrm)


###################################################
### code chunk number 11: PsychFun.Stex:611-646
###################################################
opar <- par(mfrow = c(1, 2))
with(SHR2, plot(p/100 ~ log(Q), xlab = "Log(Quanta/flash)",
	ylab = "Proportion seen", type = "n")
	)
qq <- seq(3, 6, len = 200)	
pred.nrm <- pnorm((qq - SHR2.norm$par[1])/SHR2.norm$par[2])
lines(qq, pred.nrm, lwd = 2)
lpois <- function(p, d) {
   -sum(d$nyes * ppois(p[2], d$Q/p[1], lower.tail = FALSE, 
       log.p = TRUE) + d$nno * ppois(p[2], d$Q/p[1], 
       lower.tail = TRUE, log.p = TRUE))
	}
SHR2.pois <- optim(c(18, 5), lpois, d = SHR2)
pred.poi <- ppois(SHR2.pois$par[2], exp(qq - log(SHR2.pois$par[1])), 
	lower.tail = FALSE)
lines(qq, pred.poi, lty = 2, lwd = 3)
with(SHR2, points(p/100 ~ log(Q), pch = 21, bg = "white"))
title(main =  "Obs:  SH;  Run:  2")
mtext("a", cex = 2, adj = 0, line = 0.5)
###
nn <- seq(1, 20)
lpois1 <- function(q, p, d) lpois(c(q, p), d)
SHR2.prof <- sapply(nn, function(iy) {
	unlist(
		optimize(lpois1, c(1, 1000), p = iy,
		d = SHR2) )#$objective
})
plot(nn, SHR2.prof[2, ], xlim = c(0, 21),
	xlab = expression(paste("h", nu)),
	ylab = "-Log Likelihood")
points(which.min(SHR2.prof[2, ]), 
	SHR2.prof[2, which.min(SHR2.prof[2, ])],  
	pch = 16, col = "black")
mtext("b", cex = 2, adj = 0, line = 0.5)
par(opar)


###################################################
### code chunk number 12: PsychFun.Stex:684-706
###################################################
opar <- par(mar = c(5, 6, 4, 2) + 0.1)
x <- exp(seq(-1.3, 3.3, len = 100))
plot(x, ppois(9, x, lower.tail = FALSE),
	 log = "x", type = "n",
	 xlab = expression(paste(
	 		italic("Average Number h"), nu, 
	 		italic(" per flash"))),
	 ylab = expression(
	 		paste(italic("Probability of n or more h"), nu, 
	 		italic(" per flash"))),
	 cex.lab = 1.8, cex.axis = 1.35)
abline(h = seq(0, 1, len = 6), v = c(0.5, 1, 2, 5, 10, 20), 
	col = "grey")
polygon(c(x -0.25, rev(x)), c(ppois(1, x, lower.tail = FALSE), 
		rev(ppois(10, x, lower.tail = FALSE))), col = "white", 		border = "white")
for (y in 1:9) 
	lines(x, ppois(y, x, lower.tail = FALSE))	
text(c(1.6, 2.45, 3.2, 3.8, 4.55, 5.3, 6.15, 6.8, 7.5),
	 c(0.55, 0.5, 0.45, 0.4, 0.35, 0.325, 0.3, 0.28, 0.25), 
	 c("n = 1", "2", "3", "4", "5", "6", "7", "8", "9"),
	 cex = 1.1, adj = 1)
par(opar)


###################################################
### code chunk number 13: PsychFun.Stex:738-743
###################################################
lpois <- function(p, d) { -sum(d$nyes * ppois(p[2], d$Q/p[1], 
   lower.tail = FALSE,  log.p = TRUE) + 
   d$nno * ppois(p[2], d$Q/p[1], 
      lower.tail = TRUE, log.p = TRUE))
   }


###################################################
### code chunk number 14: PsychFun.Stex:759-760
###################################################
SHR2.pois <- optim(c(18, 5), lpois, d = SHR2)


###################################################
### code chunk number 15: PsychFun.Stex:771-774 (eval = FALSE)
###################################################
## pred.poi <- ppois(SHR2.pois$par[2], 
##   exp(qq - log(SHR2.pois$par[1])), lower.tail = FALSE)
## lines(qq, pred.poi, lty = 2, lwd = 2)


###################################################
### code chunk number 16: PsychFun.Stex:778-779
###################################################
SHR2.pois$par


###################################################
### code chunk number 17: PsychFun.Stex:791-804 (eval = FALSE)
###################################################
## nn <- seq(1, 20)
## lpois1 <- function(q, p, d) lpois(c(q, p), d)
## SHR2.prof <- sapply(nn, function(iy) {
##   unlist(
##      optimize(lpois1, interval = c(1, 1000), p = iy,
##         d = SHR2) )
## })
## plot(nn, SHR2.prof[2, ], xlim = c(0, 21),
## 	xlab = expression(paste("h", nu)),
## 	ylab = "-Log Likelihood")
## points(which.min(SHR2.prof[2, ]), 
## 	SHR2.prof[2, which.min(SHR2.prof[2, ])],  
## 	pch = 16, col = "black")


###################################################
### code chunk number 18: PsychFun.Stex:832-844
###################################################
with(SHR2, plot(p/100 ~ log(Q), xlab = "Log(Quanta/flash)",
	ylab = "Proportion seen", cex = 1.5, type = "n")
	)
qq <- seq(3, 6, len = 200)	
for (ix in 5:7) {
pred.poi <- ppois(ix, exp(qq - log(SHR2.prof[1, ix])), 
	lower.tail = FALSE)
lines(qq, pred.poi, lty = ix - 4, lwd = 3)
}
with(SHR2, points(p/100 ~ log(Q), pch = 21, 
	bg = "white", cex = 1.5))
legend(3.65, 0.9, 5:7, lty = 1:3, lwd = 3, bty = "n")


###################################################
### code chunk number 19: PsychFun.Stex:1044-1046
###################################################
data(Vernier)
names(Vernier)


###################################################
### code chunk number 20: PsychFun.Stex:1091-1100
###################################################
panel.psyfun <- function(x, y, n, lnk = "logit", ...) {
	xy.glm <- glm(cbind(n * y, n * (1 - y)) ~ x, 
			binomial(lnk))
	rr <- current.panel.limits()$xlim
	xx <- seq(rr[1], rr[2], len = 100)
	yy <- predict(xy.glm, data.frame(x = xx),  
			type = "response")
	panel.lines(xx, yy,  ...)
	}


###################################################
### code chunk number 21: PsychFun.Stex:1106-1115
###################################################
print(
	xyplot(Pc ~ Phaseshift | Direction * WaveForm * TempFreq, 
	data = Vernier, layout = c(4, 2), aspect = "xy",
	xlab = "Phase Shift (deg)", ylab = "Proportion Upward",
	panel = function(x, y,  ...){
		panel.xyplot(x, y, ...)
		panel.psyfun(x, y, 20, lnk = "probit", lwd = 2, ...)
		}
	))


###################################################
### code chunk number 22: PsychFun.Stex:1146-1155
###################################################
v <- list()
v[["logit"]] <- glm(cbind(NumUpward, NumDownward) ~ 
	Direction * WaveForm * TempFreq * Phaseshift, 
	binomial, Vernier)
v[["probit"]] <- update(v[[1]], family = binomial(probit))
v[["cauchit"]] <- update(v[[1]], family = binomial(cauchit))
v[["weibull"]] <- update(v[[1]], family = binomial(cloglog),
  . ~ Direction * WaveForm * TempFreq * 
       log(Phaseshift + 50.1))


###################################################
### code chunk number 23: PsychFun.Stex:1173-1174
###################################################
sapply(v, AIC)


###################################################
### code chunk number 24: PsychFun.Stex:1210-1214
###################################################
library(MASS)
dropterm(v[["probit"]], test = "Chisq")
v.prob2 <- update(v[["probit"]], . ~ . - 
   Direction:WaveForm:TempFreq:Phaseshift)


###################################################
### code chunk number 25: PsychFun.Stex:1225-1229
###################################################
v.prob2.anova <- anova(v.prob2, test = "Chisq")
data.frame(terms =  rownames(v.prob2.anova), 
	P = v.prob2.anova[, 5] )
dropterm(v.prob2, test = "Chisq")


###################################################
### code chunk number 26: PsychFun.Stex:1274-1276
###################################################
v.prob3 <- update(v.prob2, . ~ . - 
  Direction:WaveForm:(TempFreq + Phaseshift))


###################################################
### code chunk number 27: PsychFun.Stex:1283-1284
###################################################
dropterm(v.prob3, test = "Chisq")


###################################################
### code chunk number 28: PsychFun.Stex:1306-1307
###################################################
v.prob4 <- update(v.prob3, . ~ . - Direction:WaveForm)


###################################################
### code chunk number 29: PsychFun.Stex:1342-1346
###################################################
v.prob5 <- update(v.prob4, . ~ 
	Phaseshift:((TempFreq + Direction + WaveForm) +  
	 + TempFreq:(Direction + WaveForm)) - 1)
anova(v.prob5, v.prob4, test = "Chisq")


###################################################
### code chunk number 30: PsychFun.Stex:1370-1371 (eval = FALSE)
###################################################
## Pc ~ (WaveForm + Direction):TempFreq:PhaseShift


###################################################
### code chunk number 31: PsychFun.Stex:1384-1385
###################################################
class(v.prob5)


###################################################
### code chunk number 32: PsychFun.Stex:1403-1412
###################################################
opar <- par(mfrow = c(1, 2))
levs <- psyF <- vector("list", 2)
names(psyF) <- names(levs) <- c("5", "20")
for (nlevs in names(levs)) {
	levs[[nlevs]] <- round(10^seq(-2, 0, len = as.integer(nlevs)), 3)
	psyF[[nlevs]] <- pnorm(levs[[nlevs]], mean = 0.2, sd = 0.2)
	plot(levs[[nlevs]], psyF[[nlevs]], log = "x")
}
par(opar)


###################################################
### code chunk number 33: PsychFun.Stex:1426-1434
###################################################
Ntrials <- 100
indiv.lst <- indiv.glm <-  vector("list", 2)
names(indiv.lst) <- names(indiv.glm) <- names(levs)
for (nlevs in names(levs)){
	Resp <- rbinom(Ntrials * length(psyF[[nlevs]]), 1, psyF[[nlevs]])
	indiv.lst[[nlevs]] <- data.frame(resp = Resp, levs = levs[[nlevs]])
	indiv.glm[[nlevs]] <- glm(resp ~ levs, binomial(probit), indiv.lst[[nlevs]])
}


###################################################
### code chunk number 34: PsychFun.Stex:1443-1447
###################################################
opar <- par(mfrow = c(2, 2))
for(nlevs in names(levs)) plot(indiv.glm[[nlevs]], which = 1:2, 
	main = paste("Number of levels = ", nlevs))
par(opar)


###################################################
### code chunk number 35: PsychFun.Stex:1477-1484
###################################################
GrpResp <- Grp.glm <- vector("list", 2)
names(GrpResp) <- names(Grp.glm) <- names(levs)
for (nlevs in names(levs)) {
	GrpResp[[nlevs]] <- t(with(indiv.lst[[nlevs]], 
		table(resp, levs)))
	Grp.glm[[nlevs]] <- glm(GrpResp[[nlevs]] ~ 
		levs[[nlevs]], binomial(probit))}


###################################################
### code chunk number 36: PsychFun.Stex:1496-1500
###################################################
opar <- par(mfrow = c(2, 2))
for(nlevs in names(levs)) plot(Grp.glm[[nlevs]], which = 1:2, 
	main = paste("Number of levels = ", nlevs))
par(opar)


###################################################
### code chunk number 37: PsychFun.Stex:1531-1534
###################################################
indiv.diags <- lapply(names(levs), 
	function(x) binom.diagnostics(indiv.glm[[x]], 
		nsim = 10000))


###################################################
### code chunk number 38: PsychFun.Stex:1536-1537
###################################################
names(indiv.diags) <- c("5", "20")


###################################################
### code chunk number 39: PsychFun.Stex:1542-1543
###################################################
str(indiv.diags[[1]])


###################################################
### code chunk number 40: PsychFun.Stex:1557-1558
###################################################
          plot(indiv.diags[[1]], cex = 0.5)


###################################################
### code chunk number 41: PsychFun.Stex:1568-1569
###################################################
           plot(indiv.diags[[2]], cex = 0.5)


###################################################
### code chunk number 42: PsychFun.Stex:1595-1596
###################################################
sapply(indiv.diags, "[[", "p")


###################################################
### code chunk number 43: PsychFun.Stex:1631-1638
###################################################
x <- 10^seq(-2.5, -0.5, len = 6)
pr <- pweibull(x, 3, 0.075)
Trials <- 30
set.seed(16121952)
ny <- rbinom(length(pr), 30, pr)
nn <- 30 - ny
res <- glm(cbind(ny, nn) ~ log10(x), binomial(cloglog))


###################################################
### code chunk number 44: PsychFun.Stex:1663-1667
###################################################
plot(x, pr, log = "x")
xx <- 10^seq(-2.5, -0.5, len = 100)
pred <- predict(res, newdata = data.frame(x = xx), type = "response")
lines(xx, pred)


###################################################
### code chunk number 45: PsychFun.Stex:1713-1721
###################################################
nd <- expand.grid(Phaseshift = seq(-50, 50, len = 200),
	WaveForm = c("Sine", "Square"),
	TempFreq = factor(c(2, 8)), 
	Direction = c("Downward", "Upward"))
v.pred <- predict(v.prob5, newdata = nd, type = "response")
nd$pred <- v.pred
nd$ID <- with(nd, WaveForm:TempFreq:Direction)
Vernier$ID <- with(Vernier, WaveForm:TempFreq:Direction)


###################################################
### code chunk number 46: PsychFun.Stex:1736-1748
###################################################
print(
xyplot(Pc ~ Phaseshift | Direction + WaveForm + TempFreq, 
	subscripts = TRUE, ID = Vernier$ID, aspect = "xy",
	data = Vernier, layout = c(4, 2), 
	panel = function(x, y, subscripts, ID,  ...){
		panel.xyplot(x, y, ...)
		panel.abline(v = 0, h = 0.5, col = "grey", ...)
		which <- unique(ID[subscripts])
		llines(nd$Phaseshift[nd$ID == which], 
			nd$pred[nd$ID == which], lwd = 2, ...)		
		}
	) )



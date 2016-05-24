### R code from vignette source 'SDT/SDT.Stex'

###################################################
### code chunk number 1: SDT.Stex:7-26
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
### code chunk number 2: SDT.Stex:29-30
###################################################
library(MPDiR)


###################################################
### code chunk number 3: SDT.Stex:115-125
###################################################
N <- 150
pH <- 0.9
pFA <- 0.1
H <- rbinom(1, N, pH)
M <- N - H
FA <- rbinom(1, N, pFA)
CR <- N - FA
(SDT_Resp.tab <-  as.table(matrix(c(H, M, FA, CR), nc = 2, 
    dimnames = list(Resp = c("Yes", "No"), 
    Stim = c("Signal", "Not Signal")))) )


###################################################
### code chunk number 4: SDT.Stex:206-219
###################################################
opar <- par(mfrow = c(1, 2), lwd = 3)
x <- seq(-3, 5, len = 300)
plot(x, dnorm(x), type = "l", ylim = c(0, 0.45), cex.axis = 1.5, cex.lab = 1.3)
lines(x, dnorm(x, mean = 2))
arrows(0, dnorm(0) + 0.01, 2, dnorm(0) + 0.01, angle = 90, code = 3, length = 0.05)
text(1, dnorm(0) + 0.0375, "d'", cex = 1.8)
segments(0, 0, 0, dnorm(0), lty = 3)
segments(2, 0, 2, dnorm(0), lty = 3)
mtext("a", side = 3, cex = 2, adj = 0, line = 0.5)
plot(x, dnorm(x), type = "l", ylim = c(0, 0.45), cex.axis = 1.5, cex.lab = 1.3)
lines(x, dnorm(x, mean = 2, sd = 1.4) * 1.4)  # LTM: it was hard to se the difference with 1.2 so I changed to 1.4
mtext("b", side = 3, cex = 2, adj = 0, line = 0.5)
par(opar)


###################################################
### code chunk number 5: SDT.Stex:318-326
###################################################
dp <- 1 ; crit <- 0.1
nS <- nN <- 1000
 pFA  <- 1 - pnorm(crit)
 pH <- 1 - pnorm(crit - dp)
 FA <- rbinom(1,nN, pFA)
 H <- rbinom(1, nS, pH)
 CR <- nN - FA
 M <- nS - H


###################################################
### code chunk number 6: SDT.Stex:340-352
###################################################
simulate.SDT <- function(dp, crit, nS, nN = NULL){
   if(missing(nN)) nN <- nS
   pFA <- 1 - pnorm(crit)
   pH <- 1 - pnorm(crit - dp)
   FA <- rbinom(1,nN, pFA)
   H <- rbinom(1, nS, pH)
   CR <- nN - FA
   M <- nS - H
   res <- matrix(c(H, M, FA, CR), nc = 2,
      dimnames = list(c("Yes", "No"), c("S", "NS")))
   as.table(res)
}


###################################################
### code chunk number 7: SDT.Stex:366-368
###################################################
( SDTsim.lst <- lapply(rep(1, 3), simulate.SDT,
	 crit = 0.5, nS = 1000) )


###################################################
### code chunk number 8: SDT.Stex:389-390
###################################################
 ( PrSDTsim.lst <-  lapply(SDTsim.lst, "/", 1000) )


###################################################
### code chunk number 9: SDT.Stex:399-402
###################################################
c(c = qnorm(1 - PrSDTsim.lst[[1]][1, 2]), 
  dp = qnorm(PrSDTsim.lst[[1]][1, 1]) - 
       qnorm( PrSDTsim.lst[[1]][1, 2]) )


###################################################
### code chunk number 10: SDT.Stex:453-455
###################################################
N <- 1000
Stim <- factor(rep(c("A", "B"), each = N/2))


###################################################
### code chunk number 11: SDT.Stex:466-467
###################################################
dc <- c(1.3, 0.5)


###################################################
### code chunk number 12: SDT.Stex:483-486
###################################################
z2dc <- matrix(c(1, 0, -1, -1), 2, 2)
dc2z <- solve(z2dc) 
p <- pnorm(dc2z %*% dc)


###################################################
### code chunk number 13: SDT.Stex:493-497
###################################################
Resp <- factor( ifelse(Stim == "B", rbinom(N, 1, p[1]), 
	rbinom(N, 1, p[2])), labels = LETTERS[1:2])
Resp.tab <- table(Resp, Stim)[2:1, 2:1]
z2dc %*% qnorm(Resp.tab[1, ]/(N/2))


###################################################
### code chunk number 14: SDT.Stex:517-524
###################################################
dpsim <- function(dc, N, n = 10000) {
  z2dc <- matrix(c(1, -0.5, -1, -0.5), 2, 2)
  dc2z <- solve(z2dc)
  p <- pnorm(dc2z %*% dc)
  PH <- rbinom(n, N/2, p[2])/N
  PF <- rbinom(n, N/2, p[1])/N
  apply(qnorm(cbind(PH, PF)), 1, diff)}


###################################################
### code chunk number 15: SDT.Stex:536-542
###################################################
dp.qu <- sapply(c(100, 1000, 10000), function(x) { 
   dpsim(c(0, 0), x, 10000)
   })
colnames(dp.qu) <- c(100, 1000, 10000)
dhists <- data.frame(dp = as.vector(dp.qu), 
   Trials = rep(c(100, 1000, 10000), each = 10000))


###################################################
### code chunk number 16: SDT.Stex:555-556
###################################################
apply(dp.qu, 2, quantile, c(0.025, 0.975))


###################################################
### code chunk number 17: SDT.Stex:567-573
###################################################
print(
    histogram(~ dp | factor(Trials), data = dhists,
    	layout = c(3, 1), breaks = 25,
	strip = strip.custom(style = 5),
	xlab = expression("d'"), type = "count")
)


###################################################
### code chunk number 18: SDT.Stex:619-620 (eval = FALSE)
###################################################
## Resp ~ Stim


###################################################
### code chunk number 19: SDT.Stex:629-632
###################################################
tea.glm <- glm(factor(Resp) ~ factor(Stim), 
  family = binomial("probit"))
coef(summary(tea.glm))


###################################################
### code chunk number 20: SDT.Stex:677-684
###################################################
X <- model.matrix(~ Stim)
ll <- function(b, X) {
	p <- pnorm(X %*% b)
	-sum(ifelse(Resp == "B", log(p), log(1 - p)))
	}
tea.opt <- optim(c(0, 1), ll, X = X)
tea.opt$par


###################################################
### code chunk number 21: SDT.Stex:691-704
###################################################
Recognition <- read.table(zz <- textConnection(
"Stim	Number	Resp	Cond
Old		69		Yes		Normal
New		31		Yes		Normal
Old		31		No		Normal
New		69		No		Normal
Old		89		Yes		Hypnotized
New		59		Yes		Hypnotized
Old		11		No		Hypnotized
New		41		No		Hypnotized"), TRUE)
Recognition$Cond <- relevel(Recognition$Cond, "Normal")
close(zz)
rm(zz)


###################################################
### code chunk number 22: SDT.Stex:711-712
###################################################
Recognition


###################################################
### code chunk number 23: SDT.Stex:742-745
###################################################
Recog.glm1 <- glm(Resp ~ Stim, binomial(probit),
	Recognition, weights = Number)
coef(summary(Recog.glm1))


###################################################
### code chunk number 24: SDT.Stex:765-767
###################################################
Recog.glm2 <- update(Recog.glm1, . ~ . + Cond)
coef(summary(Recog.glm2))


###################################################
### code chunk number 25: SDT.Stex:786-796
###################################################
res <- c(0.497, -0.497, 1.224, 0.229)
Stim <- rep(c("Old", "New"), 2)
Cond <- relevel(factor(rep(c("Norm", "Hypno"), 
	each = 2)), "Norm")
interaction.plot(Stim, Cond, res, type = "b", ylab = expression(italic(z)), 
	cex.lab = 1.5, pch = 16:17, cex = 1.3)
arrows(2, res[2], 2, res[1], code = 3, angle = 90)
text(2.1, 0, -diff(res[3:4]))
arrows(1, res[2], 1, res[4], code = 3, angle = 90)
text(1.1, sum(res[c(2, 4)])/2, diff(res[c(2, 4)]))


###################################################
### code chunk number 26: SDT.Stex:809-811
###################################################
Recog.glm3 <- update(Recog.glm2, . ~ . + Cond:Stim)
coef(summary(Recog.glm3))


###################################################
### code chunk number 27: SDT.Stex:831-832
###################################################
anova(Recog.glm1, Recog.glm2, Recog.glm3, test = "Chisq")


###################################################
### code chunk number 28: SDT.Stex:836-837
###################################################
 AIC(Recog.glm1, Recog.glm2, Recog.glm3)


###################################################
### code chunk number 29: SDT.Stex:858-860
###################################################
CmpltSep.df <- data.frame(Present = c(10000, 100),
	Absent = c(0, 10000), Stim = factor(c(1, 0)))


###################################################
### code chunk number 30: SDT.Stex:869-871
###################################################
cmpsep.1 <- glm(cbind(Present, Absent) ~ Stim, binomial, 
	CmpltSep.df)


###################################################
### code chunk number 31: SDT.Stex:878-879
###################################################
coef(summary(cmpsep.1))


###################################################
### code chunk number 32: SDT.Stex:883-886
###################################################
cmpsep.2 <- glm(cbind(Present, Absent) ~ Stim, 
	binomial(probit), CmpltSep.df)
coef(summary(cmpsep.2))


###################################################
### code chunk number 33: SDT.Stex:899-902
###################################################
cmpsep.3 <- glm(cbind(Present, Absent + 1) ~ Stim, 
	binomial(probit), CmpltSep.df) 
coef(summary(cmpsep.3))


###################################################
### code chunk number 34: SDT.Stex:929-943
###################################################
x <- seq(0, 1, len = 500)
xy <- expand.grid(pFA = x, pH = x)
xy$dp <- with(xy, qnorm(pH) - qnorm(pFA))
print(
contourplot(dp ~ pFA + pH, data = xy, 
	subscripts = TRUE,
	panel = function(x, y, z, subscripts, ...){
		panel.levelplot(x, y, z, subscripts, at = 0:3,
			labels = paste("d' =", 0:3),
			contour = TRUE, region = FALSE)
		panel.points(FA/N, H/N,  
			pch = 16, col = "black", cex = 1.25)
		panel.abline(v = 0.22, h = 0.9, lty = 2)
	}) ) 


###################################################
### code chunk number 35: SDT.Stex:1029-1046
###################################################
dp <- 1 ; sig <- 2;  logbeta <- 2
X <- seq(-5, 5, len = 200)
lamX <-  function(X, dp, sig)
   -((1 - sig^2) * X^2 - 2 * dp * X +
        dp^2 + 2 * sig^2 * log(sig))/(2 * sig^2)
plot(X, lamX(X, dp, sig), type = "n", ylab = "-log Likelihood")
rts <- as.real( polyroot(
  -c( dp * dp + 2 * sig^2 * log(sig) + logbeta * 2 * sig^2, 
  -2 * dp, 1 - sig^2) / (2 * sig^2) ))
rect(-6, 2, rts[2], 11, col = "lightgrey", border = NA)
rect(rts[1], 2, 6, 11, col = "lightgrey", border = NA)
lines(X, lamX(X, dp, sig))
abline(2, 0, lty = 2)
abline(v =c (rts[2], rts[1]), lty = 2)
text(3.5, 8, expression(paste(plain(log), beta > 2)))
text(-4.5, 8, expression(paste(plain(log), beta > 2)))
text(4, 2.3, expression(paste(plain(log), beta == 2)))


###################################################
### code chunk number 36: SDT.Stex:1098-1110
###################################################
 UVGSDTcrit <- function(dp, sig, logbeta){
    TwoSigSq <- 2 * sig^2
    minLam <- optimize(function(X, dp, sig)
      -((1 - sig^2) * X^2 - 2 * dp * X +
        dp^2 + TwoSigSq * log(sig))/TwoSigSq, c(-10, 10),
        dp = dp, sig = sig)$objective
   if (logbeta < minLam) warning("complex roots") 	
   cf <- -c(dp^2 + TwoSigSq * log(sig) + logbeta * TwoSigSq,
      -2 * dp, 1 - sig^2)/TwoSigSq
   proot <- polyroot(cf)
   sort(Re(proot))
 }


###################################################
### code chunk number 37: SDT.Stex:1112-1125
###################################################
rUVGSDT <- function(dp, sig, logbeta, Nn, Ns = NULL, 
		aggregate = TRUE){
	if (is.null(Ns)) Ns <- Nn
	crit <- UVGSDTcrit(dp, sig, logbeta)
	Rn <- rnorm(Nn) 
	Rs <- rnorm(Ns, dp, sig)
	RespN <- (Rn < crit[1]) | (Rn > crit[2])
	RespS <- (Rs < crit[1]) | (Rs > crit[2])
	if (aggregate)
	    c(pFA = sum(RespN)/Nn, pH = sum(RespS)/Ns) else
	    data.frame(Resp = c(RespN, RespS),
	    	Stim = factor(rep(0:1, c(Nn, Ns))))
}


###################################################
### code chunk number 38: SDT.Stex:1139-1172
###################################################
dp <- 1; sig <- 1.5; logbeta <- 1
xx <- seq(-5, 5, len = 200)
opar <- par(mfrow = c(1, 3), pty = "s" )
plot(xx, dnorm(xx), type = "l", xlab = "X")
lines(xx, dnorm(xx, mean = dp, sd = sig))
c12 <-UVGSDTcrit(dp, sig, logbeta)
abline(v = c12, lty = 2)
mtext(c(expression(c[1]), expression(c[2])), 3, at = c12, line = 0.3)
mtext("a", side = 3, line  = 1, adj = 0)
minLam <- optimize(lamX, c(-10, 10), 
	dp = dp, sig = sig)$objective + 1e-6
p <- t(sapply(seq(minLam, 10, len = 100), function(LBeta){
        c12 <- UVGSDTcrit(dp, sig, LBeta)
        1 - c(diff(pnorm(c12)), diff(pnorm(c12,dp,sig)))   }
) )
SimPt1 <- rUVGSDT(dp, sig, 0, 1000)
SimPt2 <- rUVGSDT(dp, sig, -0.5, 1000)
plot(p, type = "l",  xlim = c(0, 1), ylim = c(0, 1),
	xlab = expression(p[FA]), ylab = expression(p[H]))
points(c(SimPt1[1], SimPt2[1]), c(SimPt1[2], SimPt2[2]),
	pch = c(16, 1), col = "black")
abline(0, 1, lty = 2)
abline(v = 0.5, h = 0.5, col = "grey")
mtext("b", side = 3, line  = 1, adj = 0)
plot(qnorm(p), type = "l", xlim = c(-1, 1), ylim = c(-1, 1),
	xlab = expression(z[FA]), ylab = expression(z[H]))  
points(qnorm(c(SimPt1[1], SimPt2[1])), 
	qnorm(c(SimPt1[2], SimPt2[2])),
	pch = c(16, 1), col = "black")
abline(0, 1, lty = 2)
abline(v = 0, h = 0, col = "grey") 
mtext("c", side = 3, line  = 1, adj = 0)
par(opar)


###################################################
### code chunk number 39: SDT.Stex:1195-1208 (eval = FALSE)
###################################################
## rUVGSDT <- function(dp, sig, logbeta, Nn, Ns = NULL, 
## 		aggregate = TRUE){
## 	if (is.null(Ns)) Ns <- Nn
## 	crit <- UVGSDTcrit(dp, sig, logbeta)	
## 	Rn <- rnorm(Nn) 
## 	Rs <- rnorm(Ns, dp, sig)
## 	RespN <- (Rn < crit[1]) | (Rn > crit[2])
## 	RespS <- (Rs < crit[1]) | (Rs > crit[2])
## 	if (aggregate)
## 	    c(pFA = sum(RespN)/Nn, pH = sum(RespS)/Ns) else
## 	    data.frame(Resp = c(RespN, RespS),
## 	    	Stim = factor(rep(0:1, c(Nn, Ns))))
## }


###################################################
### code chunk number 40: SDT.Stex:1264-1266
###################################################
dp <- 1; sig <- 1.5; N <- 5000
lb <- c(1.5, 1)


###################################################
### code chunk number 41: SDT.Stex:1274-1275
###################################################
set.seed(16121952)


###################################################
### code chunk number 42: SDT.Stex:1277-1281
###################################################
UVGsim.df <- do.call(rbind, lapply(lb, function(LB) 
	rUVGSDT(dp, sig, LB, N, agg = FALSE))) 
UVGsim.df$Cond <- factor(rep(paste("C", seq(length(lb)), 
	sep = ""), each = 2 * N))


###################################################
### code chunk number 43: SDT.Stex:1297-1319
###################################################
lUVG <- function(p, d){ #p <-c(dp, sig, logbeta_1,_2,...,_n)
	sig <- p[2] 
	if (length(p[-(1:2)]) != length(levels(d[[3]])))
		stop("Number of initial estimates of 
		log beta not equal to number of levels of 
		conditions in the data!")
	Pr <- sapply(p[-(1:2)], function(lb, pr) {
		crit <- UVGSDTcrit(pr[1], sig, lb)
		1 - c(diff(pnorm(crit)), 
			diff(pnorm(crit, pr[1], pr[2])))
		}, pr = c(p[1], sig))
	ll <- sapply(seq_along(levels(d[[3]])), function(Cd){
		with(subset(d, Cond == levels(d$Cond)[Cd]),  
			ifelse(Stim == "1", 
			      Resp * log(Pr[2, Cd]) + 
			      (1 - Resp) * log(1 - Pr[2, Cd]),
			      Resp * log(Pr[1, Cd]) + 
			      (1 - Resp) * log(1 - Pr[1, Cd])
			))
		})
	-sum(unlist(ll))
}


###################################################
### code chunk number 44: SDT.Stex:1339-1351
###################################################
minLam <- optimize(function(X, dp, sig)
      -((1 - sig^2) * X^2 - 2 * dp * X +
        dp^2 + 2 * sig^2 * log(sig))/(2 * sig^2), 
        c(-20, 20), dp = dp, sig = sig)$objective
UVG.opt <- optim(c(1.1, 1.6, 1.4, 1.1), lUVG, d = UVGsim.df,
	method = "L-BFGS-B", hessian = TRUE,
	lower = c(0, 1, rep(minLam, length(lb))))
est <- UVG.opt$par
names(est) <- c("dprime", "sigma", 
	paste("log beta", seq(length(lb)), sep = ""))
est.se <- sqrt(diag(solve(UVG.opt$hessian)))
cbind(Estimate = est, SE = est.se, z = est/est.se)


###################################################
### code chunk number 45: SDT.Stex:1393-1397
###################################################
xx = seq(0, 5, len = 100)
plot(xx, dexp(xx), type = "l",
	xlab = "T", lwd = 2, ylim = c(0, 2))
lines(xx, dexp(xx, rate = 2), lty = 2, lwd = 2)


###################################################
### code chunk number 46: SDT.Stex:1469-1471
###################################################
N <- 1000
Stim <- factor(rep(Stim, 1000), labels = 1:2)


###################################################
### code chunk number 47: SDT.Stex:1475-1480
###################################################
dc <- c(1.4, 0.7)
p <- pnorm(dc2z %*% dc)
Resp <- factor( ifelse(Stim == "1", rbinom(N, 1, p[1]), 
	rbinom(N, 1, p[2])))
Resp.tab <- table(Resp, Stim)


###################################################
### code chunk number 48: SDT.Stex:1483-1486
###################################################
Pc <- sum(Resp.tab[2:3])/sum(Resp.tab)
library(psyphy)
dprime.mAFC(Pc, m = 2)


###################################################
### code chunk number 49: SDT.Stex:1490-1491
###################################################
(z2dc %*% qnorm(Resp.tab[c(2, 4)]/(2 * N)))[1]/sqrt(2) 


###################################################
### code chunk number 50: SDT.Stex:1510-1513
###################################################
 data(ecc2)
ID12 <- subset(ecc2, task == "ID" & Size == 12.4)
Pc <- with(ID12, Correct/(Correct + Incorrect))


###################################################
### code chunk number 51: SDT.Stex:1517-1518
###################################################
sapply(Pc, dprime.mAFC, m = 4)


###################################################
### code chunk number 52: SDT.Stex:1528-1530
###################################################
dpvec.mAFC <- Vectorize(dprime.mAFC, "Pc")
dpvec.mAFC(Pc, 4)


###################################################
### code chunk number 53: SDT.Stex:1629-1637
###################################################
N <- 1000 
dp <- 1 
Noise <- rnorm(N)  
IntResp <- rep(c(dp, 0), each = N/2) + Noise
Stim <- factor(rep(1:0, each = N/2))
IntCrit <- c(-Inf, -2:2, Inf)
Ratings <- cut(IntResp, IntCrit)
levels(Ratings) <- 1:6 


###################################################
### code chunk number 54: SDT.Stex:1648-1653
###################################################
X <- seq(-3, 4, len = 1000)
plot(X, dnorm(X), type = "l")
lines(X, dnorm(X, mean = 1))
abline(v = -2:2, lty = 2)
text(seq(-2.5, 2.5), 0.385, 1:6, cex = 1.3)


###################################################
### code chunk number 55: SDT.Stex:1672-1675
###################################################
( Rat.tab <- table(Stim, Ratings) )
( Rat.proptab <- apply(Rat.tab[, 6:1], 1, cumsum)[6:1, ]  / 
	(N/2) )


###################################################
### code chunk number 56: SDT.Stex:1686-1690
###################################################
Rat.lm <- lm(qnorm(Rat.proptab[-1, 2]) ~ 
	offset(qnorm(Rat.proptab[-1, 1])))
zfa <- seq(-4, 4, len = 100)
RatROC.fit <- pnorm(cbind(1, zfa) %*% c(coef(Rat.lm), 1) )


###################################################
### code chunk number 57: SDT.Stex:1693-1704
###################################################
ordRat <- ordered(Ratings)
cumRat <- as.vector(sapply(1:5, function(x) ordRat <= x))
X <- matrix(
c(rep(c(1, 0, 0, 0, 0), each = N),
rep(c(0, 1, 0, 0, 0), each = N),
rep(c(0, 0, 1, 0, 0), each = N),
rep(c(0, 0, 0, 1, 0), each = N),
rep(c(0, 0, 0, 0, 1), each = N)), ncol = 5)
X <- cbind(X, -rep(Stim == "1", 5))
Rat.df <- data.frame(Resp = cumRat, X = X)
Rat.glm <- glm(Resp ~ X - 1, binomial(probit), Rat.df)


###################################################
### code chunk number 58: SDT.Stex:1710-1723
###################################################
opar <- par(mfrow = c(1, 2))
plot(Rat.proptab, xlim = c(0, 1), ylim = c(0, 1),
	xlab = expression(P[FA]),
	ylab = expression(P[H]))
lines(pnorm(zfa), RatROC.fit)
lines(pnorm(zfa), pnorm(zfa + coef(Rat.glm)[6]), lty = 2, lwd = 2)
plot(qnorm(Rat.proptab), xlim = c(-3, 3), ylim = c(-3, 3),
	xlab = expression(z[FA]),
	ylab = expression(z[H]))
abline(coef(Rat.lm), 1)
abline(v = 0, h = 0, col = "grey")
abline(coef(Rat.glm)[6], 1, lty = 2, lwd = 2)
par(opar)	


###################################################
### code chunk number 59: SDT.Stex:1780-1790 (eval = FALSE)
###################################################
## ordRat <- ordered(Ratings)
## cumRat <- as.vector(sapply(1:5, function(x) ordRat <= x))
## X <- matrix(
## c(rep(c(1, 0, 0, 0, 0), each = N),
## rep(c(0, 1, 0, 0, 0), each = N),
## rep(c(0, 0, 1, 0, 0), each = N),
## rep(c(0, 0, 0, 1, 0), each = N),
## rep(c(0, 0, 0, 0, 1), each = N)), ncol = 5)
## X <- cbind(X, -rep(Stim == "1", 5))
## Rat.df <- data.frame(Resp = cumRat, X = X)


###################################################
### code chunk number 60: SDT.Stex:1812-1813
###################################################
( Rat.glm <- glm(Resp ~ . - 1, binomial(probit), Rat.df) )


###################################################
### code chunk number 61: SDT.Stex:1841-1843
###################################################
library(MASS)
polr(ordRat ~ Stim, method = "probit")


###################################################
### code chunk number 62: SDT.Stex:1862-1866
###################################################
library(ordinal)
Rat.clm0 <- clm(ordRat ~ Stim, link = "probit")
Rat.clm1 <- update(Rat.clm0, threshold = "symmetric")
Rat.clm2 <- update(Rat.clm1, threshold = "equidistant")


###################################################
### code chunk number 63: SDT.Stex:1869-1870
###################################################
sapply(list(Rat.clm0, Rat.clm1, Rat.clm2), "[[", "edf")


###################################################
### code chunk number 64: SDT.Stex:1875-1876
###################################################
anova(Rat.clm2, Rat.clm1, Rat.clm0)


###################################################
### code chunk number 65: SDT.Stex:1911-1913
###################################################
Rat.clm.UV <- update(Rat.clm0, scale = ~ Stim)
summary(Rat.clm.UV)


###################################################
### code chunk number 66: SDT.Stex:1919-1920
###################################################
exp(Rat.clm.UV$zeta)


###################################################
### code chunk number 67: SDT.Stex:1927-1928
###################################################
anova(Rat.clm.UV, Rat.clm0)


###################################################
### code chunk number 68: SDT.Stex:1956-1960
###################################################
data(Faces)
rat.tab <- with(Faces, table(sibs, SimRating))
cum.tab <- apply(rat.tab, 1, cumsum)
( cum.prop <- 1 - cum.tab/rowSums(rat.tab) )


###################################################
### code chunk number 69: SDT.Stex:1976-1991
###################################################
par(mfrow = c(1, 2))
plot(cum.prop, xlim = c(0, 1), ylim = c(0, 1),
  xlab = expression(P[FA]), ylab = expression(P[H]))
abline(0, 1, col = "grey")
dp <- mean(apply(qnorm(cum.prop[-11, ]), 1, diff))
xx <- seq(-3, 3, len = 100)
lines(pnorm(xx, lower.tail = FALSE), 
  pnorm(xx, mean = dp, lower.tail = FALSE))
mtext("a", line = 0.5, cex = 2, adj = 0)
plot(qnorm(cum.prop[-11, ]), xlim = c(-2.5, 2.5), 
  ylim = c(-2.5, 2.5), xlab = expression(z[FA]), 
  ylab = expression(z[H]))
abline(v = 0, h = 0, col = "grey")
lines(xx, xx + dp)
mtext("b", line = 0.5, cex = 2, adj = 0)


###################################################
### code chunk number 70: SDT.Stex:2011-2012
###################################################
-qnorm(cum.prop[-11, 1])


###################################################
### code chunk number 71: SDT.Stex:2024-2027
###################################################
Faces.clm1 <- clm(ordered(SimRating) ~ sibs,
	data = Faces, link = "probit")
Faces.clm2 <- update(Faces.clm1, scale = ~ sibs)


###################################################
### code chunk number 72: SDT.Stex:2030-2031
###################################################
anova(Faces.clm1, Faces.clm2)


###################################################
### code chunk number 73: SDT.Stex:2035-2036
###################################################
summary(Faces.clm1)


###################################################
### code chunk number 74: SDT.Stex:2061-2063
###################################################
addterm(Faces.clm1, scope = ~ sibs + gendiff + agediff, 
  test = "Chisq")


###################################################
### code chunk number 75: SDT.Stex:2068-2069
###################################################
summary(update(Faces.clm1, location = . ~ . + agediff))



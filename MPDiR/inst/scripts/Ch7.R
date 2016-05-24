### R code from vignette source 'MLDS/MLDS.Stex'

###################################################
### code chunk number 1: MLDS.Stex:9-28
###################################################
options(width=60, digits = 3, str = strOptions(strict.width = "cut"))
#pkg <- search()[2] 
#while (search()[2] !=  if(.Platform$GUI == "AQUA") "tools:RGUI" else "package:stats") {
##detach(pos = match(pkg, search()))
#spkg <- strsplit( pkg, ":"  )[[1]][2] 
#if (packageHasNamespace(spkg, .libPaths()[1]) ) unloadNamespace(spkg ) else #detach(pos = match(pkg, search()))
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


###################################################
### code chunk number 2: MLDS.Stex:33-34
###################################################
library(MLDS)


###################################################
### code chunk number 3: MLDS.Stex:107-128
###################################################
rc <- DefineMyScale()
N <- 100
ptSize <- 1
opar <- par(mfrow = c(3, 4))
for (ix in 1:11) {
	covm <- matrix(c(1, rep(rc[ix], 2), 1), 2, 2)
		xy <- MASS:::mvrnorm(n = N, rep(0, 2), covm, 
			empirical = TRUE)
		plot(xy, axes = FALSE, xlab = "", ylab = "", pty = "s",
			cex = ptSize, pch = 16, col = "black",
			xlim = c(-4, 4), ylim = c(-4, 4), 
			main = substitute(paste(r, " = ", rc), 
				list(rc = round(rc[ix], 2))), cex.main = 3  ) 
			}
data(kk3)
par(mar = c(7, 5, 4, 2) + 0.1)
plot(mlds(SwapOrder(kk3), method = "optim", opt.init = c(seq(0, 1, len = 11), 0.2)), type = "b",
		xlab = "r", ylab = "Difference Scale Value",
		cex.axis = 1.1, cex = 1.25, cex.lab = 2, 
		las = 1)	
par(opar)


###################################################
### code chunk number 4: MLDS.Stex:245-248
###################################################
opar <- par(mfrow = c(2, 2))
DisplayOneQuad(c(0.6, 0.9, 0.1, 0.3))
par(opar)


###################################################
### code chunk number 5: MLDS.Stex:332-333
###################################################
head(kk3)


###################################################
### code chunk number 6: MLDS.Stex:487-498
###################################################
ll.mlds <- function(p, d){
	psi <- c(0, plogis(p[-length(p)]), 1)
	sig <- exp(p[length(p)])
	resp <- d[, 1]
	stim <- as.matrix(d[, -1])
	del <- matrix(psi[stim], ncol = 4) %*%
		c(1, -1, -1, 1)/sig
	-sum(pnorm(del[resp == 1], log.p = TRUE)) -
	   sum(pnorm(del[resp == 0], lower.tail = FALSE, 
	   		log.p = TRUE))
}


###################################################
### code chunk number 7: MLDS.Stex:519-521 (eval = FALSE)
###################################################
## p <- c(qlogis(seq(0.1, 0.9, 0.1)), log(0.2)) # transformed init. values
## res <- optim(p, ll.mlds, d = kk3, method = "BFGS")


###################################################
### code chunk number 8: MLDS.Stex:615-616 (eval = FALSE)
###################################################
## glm(Resp ~ . - 1, family = binomial(probit), data = d)


###################################################
### code chunk number 9: MLDS.Stex:685-688 (eval = FALSE)
###################################################
## library(MLDS)
## Obs.out <- runQuadExperiment(DisplayTrial = "DisplayOneQuad", 
##    DefineStimuli = "DefineMyScale")


###################################################
### code chunk number 10: MLDS.Stex:740-742
###################################################
data(kk1, kk2, kk3)
kk <- SwapOrder(rbind(kk1, kk2, kk3))


###################################################
### code chunk number 11: MLDS.Stex:756-761
###################################################
kk.mlds <- mlds(kk)
summary(kk.mlds)
kkopt.mlds <- mlds(kk, method = "optim", 
	opt.init = c(seq(0, 1, len = 11), 0.2))
summary(kkopt.mlds)


###################################################
### code chunk number 12: MLDS.Stex:800-811
###################################################
opar <- par(mfrow = c(1, 2), pty = "s")
plot(kk.mlds, standard.scale = TRUE, 
	xlab = expression(r), ylab = "Difference Scale Value",
	cex.lab = 1.5, cex.axis = 1.5, las = 1)
lines(kkopt.mlds$stimulus, kkopt.mlds$pscale)
plot(kk.mlds$stimulus^2, kk.mlds$pscale/max(kk.mlds$pscale),
	xlab = expression(r^2), ylab = "Difference Scale Value", 
	cex.lab = 1.5, cex.axis = 1.5, las = 1, xlim = c(0, 1))
abline(0, 1, lty = 2, lwd = 2)
lines(kkopt.mlds$stimulus^2, kkopt.mlds$pscale)
par(opar)


###################################################
### code chunk number 13: MLDS.Stex:847-848
###################################################
kk.frm <- mlds(~ (sx/0.98)^p[1], p = c(2, 0.2), data = kk)


###################################################
### code chunk number 14: MLDS.Stex:876-878
###################################################
c(p = kk.frm$par, sigma = kk.frm$sigma)
sqrt(diag(solve(kk.frm$hess)))


###################################################
### code chunk number 15: MLDS.Stex:892-896
###################################################
ddf <- diff(c(attr(logLik(kk.frm), "df"), 
	attr(logLik(kk.mlds), "df")))
pchisq(-2 * c(logLik(kk.frm) - logLik(kk.mlds)), 
	ddf, lower.tail = FALSE)


###################################################
### code chunk number 16: MLDS.Stex:913-918
###################################################
plot(kk.mlds, standard.scale = TRUE,
	xlab = "r",
	ylab = "Difference Scale Value")
xx <- seq(0, 0.98, len = 100)
lines(xx, kk.frm$func(kk.frm$par, xx))


###################################################
### code chunk number 17: MLDS.Stex:959-969 (eval = FALSE)
###################################################
## RespFun <- function(x) x^2
## x <- seq(0, 1, len = 200)
## xy <- expand.grid(x = x, y = x)
## xy$Diff <- with(xy, ifelse(x - y >=  -0.01, 
## 	RespFun(x) - RespFun(y), NA))
## contourplot(Diff ~ x * y, xy, cut = 5,
## 	aspect = "iso",
## 	xlab = "Correlation Stimulus 1",
## 	ylab = "Correlation Stimulus 2",
## 	main = "Equi-Response Differences")


###################################################
### code chunk number 18: MLDS.Stex:990-1049
###################################################
RespFun <- function(x) x^2
x <- seq(0, 1, len = 200)
xy <- expand.grid(x = x, y = x)
xy$Diff <- with(xy, ifelse(x - y >= -0.01, 
	RespFun(x) - RespFun(y), NA))
print(
contourplot(Diff ~ x * y, xy, cut = 6,
	aspect = "iso",
	xlab = "Correlation Stimulus 1",
	ylab = "Correlation Stimulus 2",
	main = "Equi-Response Differences"), 
	more = TRUE, split = c(1, 1, 2, 2)
	)	

xy$Pdiff <-  with(xy, ifelse(x - y >= -0.01, 
pnorm(RespFun(x) - RespFun(y) - 0.4, sd = 0.1), NA))
print(
levelplot(Pdiff ~ x * y, xy, region = TRUE,
 aspect = "iso",
 xlab = "Correlation Stimulus 1",
 ylab = "Correlation Stimulus 2",
 panel = panel.levelplot.raster,
 interpolate = TRUE,
 main = expression(sigma == 0.1),
 col.regions = rev(grey(seq(0.2, 1, len = 100)))),
 more = TRUE, split = c(2, 1, 2, 2)
 )

xy$Pdiff <-  with(xy, ifelse(x - y >= -0.01, 
	pnorm(RespFun(x) - RespFun(y) - 0.4, sd = 0.2), NA))
print(
levelplot(Pdiff ~ x * y, xy, region = TRUE,
 aspect = "iso",
 xlab = "Correlation Stimulus 1",
 ylab = "Correlation Stimulus 2",
 panel = panel.levelplot.raster,
 interpolate = TRUE,
 xlim = c(0, 1), ylim = c(0, 1),
 main = expression(sigma == 0.2),
 col.regions = rev(grey(seq(0.2, 1, len = 100)))),
 more = TRUE, split = c(1, 2, 2, 2))
 
xy$Pdiff <-  with(xy, ifelse(x - y >= -0.01, 
	pnorm(RespFun(x) - RespFun(y) - 0.4, sd = 0.5), NA))
print(
levelplot(Pdiff ~ x * y, xy, region = TRUE,
 aspect = "iso",
 xlab = "Correlation Stimulus 1",
 ylab = "Correlation Stimulus 2",
 panel = panel.levelplot.raster,
 interpolate = TRUE,
 xlim = c(0, 1), ylim = c(0, 1),
 main = expression(sigma == 0.5),
 col.regions = rev(grey(seq(0.2, 1, len = 100))),
  colorkey = list(col = rev(grey(seq(0.2, 1, len = 17))),
 	 at = seq(-0.0667, 1.067, len = 18), 
	 labels = list(labels = seq(0, 1, 0.2),
 	 	at = seq(0, 1, 0.2)))),
 more = FALSE, split = c(2, 2, 2, 2))


###################################################
### code chunk number 19: MLDS.Stex:1101-1109 (eval = FALSE)
###################################################
## xy$Pdiff <-  with(xy, ifelse(x - y >= -0.01, 
## 	pnorm(RespFun(x) - RespFun(y) - 0.4, sd = 0.1), NA))
## levelplot(Pdiff ~ x * y, data = xy, 
## 	aspect = "iso",
## 	xlab = "Correlation Stimulus 1",
## 	ylab = "Correlation Stimulus 2",
## 	main = expression(sigma == 0.1),
## 	col.regions = rev(grey(seq(0.2, 1, len = 100))))


###################################################
### code chunk number 20: MLDS.Stex:1220-1233
###################################################
x <- seq(0,1, len = 11)   # 11 equally-spaced steps
psi <- x^2                       # true difference scale
sigma <- 0.2
n <- length(psi)
trials.quads <- t(combn(n, 4))    # 330 trials
N <- 1000
fit.quads <- lapply(SimMLDS(trials.quads, psi, sigma, N), 
	mlds)
cc.quads <- sapply(fit.quads, coef)
psi.quads <- apply(cc.quads, 2, 
	function(x) c(0, x/x[length(x)]))
mean.quads <- rowMeans(psi.quads)
sd.quads <- sd(t(psi.quads))


###################################################
### code chunk number 21: MLDS.Stex:1270-1285
###################################################
x <- seq(0,1, len =11)   # 11 equally-spaced steps
psi <- x^2   # true difference scale
sigma <- 0.2
n <- length(psi)
trials.triads <- t(combn(n, 3))    # 165 trials
trials.triads <- kronecker(matrix(c(1, 1), ncol = 1), 
	trials.triads)   # double number of trials
N <- 1000
fit.triads <- lapply(SimMLDS(trials.triads, psi, sigma, N), 
	mlds)
cc.triads <- sapply(fit.triads, coef)
psi.triads <- apply(cc.triads, 2, 
	function(x) c(0, x/x[length(x)]))
mean.triads <- rowMeans(psi.triads)
sd.triads <- sd(t(psi.triads))


###################################################
### code chunk number 22: MLDS.Stex:1290-1302
###################################################
opar <- par(mfrow=c(1,2))
plot(x, mean.quads, main = "Quads",
	xlab = expression(r^2),
	ylab = "Difference Scale Value")
segments(x, mean.quads + 1.96 * sd.quads, 
	x, mean.quads -  1.96 * sd.quads)
plot(x, mean.triads, main = "Triads",
	xlab = expression(r^2),
	ylab = "Difference Scale Value")
segments(x, mean.triads +  1.96 * sd.triads, 
	x, mean.triads -  1.96 * sd.triads)
par(opar)


###################################################
### code chunk number 23: MLDS.Stex:1340-1341
###################################################
kk.bt <- boot.mlds(kk.mlds, nsim = 10000)


###################################################
### code chunk number 24: MLDS.Stex:1343-1344
###################################################
str(kk.bt)


###################################################
### code chunk number 25: MLDS.Stex:1371-1379
###################################################
plot(kk.mlds, standard.scale=TRUE,
	xlab = expression(r^2),
	ylab = "Difference Scale Value")
kk.bt.res <- summary(kk.bt)
kk.mns <- kk.bt.res[, 1] 
kk.95ci <- qnorm(0.975) * kk.bt.res[, 2]
segments(kk.mlds$stimulus, kk.mns + kk.95ci, 
	kk.mlds$stimulus,  kk.mns - kk.95ci)


###################################################
### code chunk number 26: MLDS.Stex:1440-1441
###################################################
AIC(kk.mlds)


###################################################
### code chunk number 27: MLDS.Stex:1447-1449
###################################################
with(kk.mlds, (obj$null.deviance - deviance(obj)) / 
	obj$null.deviance)


###################################################
### code chunk number 28: MLDS.Stex:1479-1480
###################################################
kk.diag.prob <- binom.diagnostics(kk.mlds, 10000)


###################################################
### code chunk number 29: MLDS.Stex:1484-1485
###################################################
str(kk.diag.prob)


###################################################
### code chunk number 30: MLDS.Stex:1496-1510
###################################################
plot.mlds.diag <- function (x, alpha = 0.025, breaks = "Sturges", ...) 
{
    nsim <- dim(x$resid)[1]
    n <- dim(x$resid)[2]
    opar <- par(mfrow = c(1, 2))
    plot(sort(x$Obs.resid), (1:n - 0.5)/n, ylab = "Cumulative Density Function", 
        xlab = "Deviance Residuals", ...)
    lines(x$resid[alpha * nsim, ], (1:n - 0.5)/n, lty = 2)
    lines(x$resid[(1 - alpha) * nsim, ], (1:n - 0.5)/n, lty = 2)
    hist(x$NumRuns, xlab = "Number of Runs", main = "", breaks = breaks)
    abline(v = x$ObsRuns, lwd = 2)
    par(opar)
    invisible()
}


###################################################
### code chunk number 31: MLDS.Stex:1515-1516
###################################################
plot(kk.diag.prob, pch = ".", cex = 3)


###################################################
### code chunk number 32: MLDS.Stex:1550-1551
###################################################
kk[residuals(kk.mlds$obj) < -2.5, ]


###################################################
### code chunk number 33: MLDS.Stex:1585-1588
###################################################
library(psyphy)
kk.mlds2 <- psyfun.2asym(cbind(resp, 1 - resp) ~ . - 1, 
	data = kk.mlds$obj$data, link = probit.2asym)


###################################################
### code chunk number 34: MLDS.Stex:1598-1599
###################################################
pmc(kk.mlds)


###################################################
### code chunk number 35: MLDS.Stex:1643-1654
###################################################
opar <- par(xaxt = "n", yaxt = "n", bty = "n")
 plot(c(0, 14), c(0, 1), type = "n", xlab = "", ylab = "")
 lines(c(1, 12), c(0.5, 0.5), lwd = 5)
 points(c(3.25, 5.5, 7, 9.15, 10.25, 11.25), rep(0.5, 6), pch = 16, cex = 1.55)
 text(3.25, 0.2, "a", cex = 3)
  text(5.5, 0.2, "b", cex = 3)
 text(7, 0.2, "c", cex = 3)
 text(9.15, 0.2, "a'", cex = 3)
  text(10.25, 0.2, "b'", cex = 3)
 text(11.25, 0.2, "c'", cex = 3)
par(opar)


###################################################
### code chunk number 36: MLDS.Stex:1681-1699
###################################################
opar <- par(lwd = 3, pty = "s")
th <- seq(0, pi, len = 200)
plot(cos(th), sin(th), type = "l", axes = FALSE,
	xlab = "", ylab = "", xlim = c(-1.3, 1.3), ylim = c(-1.3, 1.3))
th <- seq(pi, 1.2*pi/2, len = 3)
lines(cos(th), sin(th))
rd <- c(pi,1.2*pi/2)
lines(cos(rd), sin(rd), col = "black", lty = 2)
points(cos(th), sin(th), pch = 16, col = "black")
text(1.2 * cos(th), 1.2 * sin(th), c("a", "b", "c"), cex = 3)
th <- seq(pi/2,0.2*pi/2, len = 3)
points(cos(th), sin(th), pch = 16, col = "black")
lines(cos(th), sin(th))
rd <- c(pi/2,0.2*pi/2)
lines(cos(rd), sin(rd), col = "black", lty = 2)
text(1.2 * cos(th), 1.2 * sin(th), c("a'", "b'", "c'"), cex = 3)
mtext("a", line = 1, at = -1, adj = 0, cex = 4)
par(opar)


###################################################
### code chunk number 37: MLDS.Stex:1702-1739
###################################################
opar <- par(lwd = 3, pty = "s")
asr <- 0.475
dd <- function(th, x, y, d) {
	sqrt((cos(th) - x)^2 + (asr *  sin(th) - y)^2)  - d }

th <- seq(0, pi, len = 200)
plot(cos(th), asr * sin(th), type = "l", axes = FALSE,
	xlab = "", ylab = "", xlim = c(-1.3, 1.3), 
	ylim = c(-1.3, 1.3))
th <- seq(pi, 5 * pi/8, len = 3) #- 1  #points a, b, c
dst <- sqrt(rowSums((apply(cbind(cos(th),  asr * sin(th)), 2, diff))^2))
points(cos(th),  asr * sin(th), pch = 16, col = "black")
lines( cos(th), asr * sin(th))
lines(cos(th[-2]), asr *  sin(th[-2]), col = "black", lty = 2)
text(1.225 * cos(th), 1.225 *  asr * sin(th), c("a", "b", "c"), cex = 3)
th2 <- th[3] - pi/16  # 2 *  pi /2 + 0.3
th2 <- c(th2, 
	uniroot(dd, c(th2, th2 - pi/6), 
	x =  cos(th2), y = asr * sin(th2), d = dst[1])$root
	)
th2 <- c(th2, 
	uniroot(dd, c(th2[2], th2[2] - pi^2/10), 
		x = cos(th2[2]), y =  asr * sin(th2[2]), d = dst[2])$root
		)
points(cos(th2),  asr * sin(th2), pch = 16, col = "black")
lines(cos(th2),  asr * sin(th2))
lines(cos(th2[-2]),  asr * sin(th2[-2]), col = "black", lty = 2)
text(1.25 * cos(th2), 1.25 *  asr * sin(th2), c("a'", "b'", "c'"), cex = 3)
ac0 <- sum(rowSums(diff(cbind(cos(th[-2]), asr * 	sin(th[-2])))^2))
ac1 <- sum(rowSums(diff(cbind(cos(th2[-2]), asr * 
	sin(th2[-2])))^2))
txt <- paste("ac = ", round(ac0, 2), ", a'c' = ", 
	round(ac1, 2), sep = "")
text(0, -0.75, txt, cex = 3)
mtext("b", line = 1, at = -1, adj = 0, cex = 4)
par(opar)



###################################################
### code chunk number 38: MLDS.Stex:1816-1817
###################################################
kk.6pt <- simu.6pt(kk.mlds, nsim = 10000, nrep = 3)


###################################################
### code chunk number 39: MLDS.Stex:1819-1820
###################################################
str(kk.6pt)


###################################################
### code chunk number 40: MLDS.Stex:1833-1835
###################################################
hist(kk.6pt$boot.samp)
abline(v = kk.6pt$lik6pt, lwd = 2)



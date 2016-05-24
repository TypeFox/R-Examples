### R code from vignette source 'PsychFun2/PsychFun2.Stex'

###################################################
### code chunk number 1: PsychFun2.Stex:9-30
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
### code chunk number 2: PsychFun2.Stex:114-116
###################################################
library(psyphy)
names(ecc2)


###################################################
### code chunk number 3: PsychFun2.Stex:126-127
###################################################
sz20.6 <- subset(ecc2, Size == 20.6)


###################################################
### code chunk number 4: PsychFun2.Stex:139-140
###################################################
qnorm(0.25)


###################################################
### code chunk number 5: PsychFun2.Stex:147-151
###################################################
q25 <- rep(qnorm(0.25), nrow(sz20.6))
det1.glm <- glm(cbind(Correct, Incorrect) ~ Contr - 1,
  binomial(probit), sz20.6, subset = task == "DET",
  offset = q25)


###################################################
### code chunk number 6: PsychFun2.Stex:155-156 (eval = FALSE)
###################################################
##  ~ Contr + offset(q25) - 1


###################################################
### code chunk number 7: PsychFun2.Stex:170-171
###################################################
det2.glm <- update(det1.glm, . ~ . + I(Contr^2))


###################################################
### code chunk number 8: PsychFun2.Stex:184-185
###################################################
det3.glm <- update(det2.glm, . ~ . - Contr)


###################################################
### code chunk number 9: PsychFun2.Stex:202-218
###################################################
p <- with(sz20.6, Correct/(Correct + Incorrect))
plot(p ~ Contr, data = sz20.6, subset = task == "DET", 
  xlab = "Contrast", ylab = "Proportion Correct",
  ylim = c(0, 1), xlim = c(0.001, 1),  log = "x", type = "n")
abline(0.25, 0, lty = 2, lwd = 2, col = "grey")  
cc <- seq(0.001, 1, len = 200)
nq25 <- rep(qnorm(0.25), length(cc))
nd <- data.frame(Contr =  cc, q25 = nq25)
lines(cc, predict(det1.glm, newdata = nd, type = "response"),
	lty = 3, lwd = 3)
lines(cc, predict(det2.glm, newdata = nd, type = "response"),
	lty = 2, lwd = 3)
lines(cc, predict(det3.glm, newdata = nd, type = "response"),
	lwd = 3)
points(p ~ Contr, data = sz20.6,  subset = task == "DET", 
  pch = 21, bg = "white", cex = 1.5)


###################################################
### code chunk number 10: PsychFun2.Stex:268-269 (eval = FALSE)
###################################################
## mafc.probit(.m = 3)


###################################################
### code chunk number 11: PsychFun2.Stex:288-294
###################################################
sz20.6$LContr <- log10(sz20.6$Contr)
m1 <- glm(cbind(Correct, Incorrect) ~ LContr, 
  binomial(mafc.probit(.m = 4)), data = sz20.6)
m2 <- update(m1, . ~ . + task)
m3 <- update(m2, . ~ . + task:LContr)
anova(m1, m2, m3,  test = "Chisq")


###################################################
### code chunk number 12: PsychFun2.Stex:332-349
###################################################
p <- with(sz20.6, Correct/(Correct + Incorrect))
pch <- 1:2
plot(p ~ Contr, data = sz20.6, pch = pch[unclass(sz20.6$task)],
  xlab = "Contrast", ylab = "Proportion Correct",
  ylim = c(0, 1), xlim = c(0.02, 0.3),  log = "x", type = "n")
xx <- seq(-1.7, -0.5, len = 100)
nd <- data.frame(LContr = rep(xx , 2),
  task = rep(c("DET", "ID"), each = 100))
pred <- predict(m2, newdata = nd, type =  "response")
abline(0.25, 0, lty = 2)
lines(10^xx, pred[1:100], lwd = 2) 
lines(10^xx, pred[-(1:100)], lwd = 2, pch =20 +  pch[unclass(sz20.6$task)],
  bg = "white")
  points(p ~ Contr, data = sz20.6, pch = 20 + pch[unclass(sz20.6$task)],
    bg = "white")
 legend(0.02, 0.9, legend = c("DET", "ID"), pch = c(21, 22), 
   pt.cex = 1.8, bty = "n")


###################################################
### code chunk number 13: PsychFun2.Stex:426-427
###################################################
set.seed(10041992)


###################################################
### code chunk number 14: PsychFun2.Stex:429-444
###################################################
b <- 3      # steepness parameter
g <- 0.02 # lower asymptote
d <- 0      # 1 - upper asymptote
a <- 0.027 # threshold
num.tr <- 100
cnt <- 10^seq(-2, -1, length = 8) 
truep <- g + (1 - g - d) * pweibull(cnt, b, a)
ny <- rbinom(length(cnt), num.tr, truep)
nn <- num.tr - ny
ny.mod <- ny
ny.mod[8] <- ny[8] - 1 
nn.mod <- num.tr - ny.mod
weib.sim <- data.frame(ny = ny, nn = nn, cnt = cnt)
weib.sim.mod <- data.frame(ny = ny.mod, nn = nn.mod, 
	cnt = cnt)


###################################################
### code chunk number 15: PsychFun2.Stex:448-451
###################################################
weib.glm <- glm(cbind(ny, nn) ~ log10(cnt), 
	binomial("cloglog"), data = weib.sim)
weib2.glm <- update(weib.glm, data = weib.sim.mod)


###################################################
### code chunk number 16: PsychFun2.Stex:466-468
###################################################
weib3.glm <- psyfun.2asym(cbind(ny.mod, nn.mod) ~ log10(cnt), 
	link = weib.2asym, init.g = 0.05, init.lam = 0.025) 


###################################################
### code chunk number 17: PsychFun2.Stex:483-496
###################################################
p.obs <- ny/(ny + nn)
p.obs.mod <- ny.mod/(ny.mod + nn.mod)
plot(cnt, p.obs, log = "x",  ylim = c(0, 1),
	ylab = "Proportion Correct", xlab = "Stimulus Intensity")
pcnt <- seq(0.01, 0.1, len = 100)
lines(pcnt, predict(weib.glm, data.frame(cnt = pcnt),
	type = "response"), lwd = 2)
lines(pcnt, predict(weib2.glm, data.frame(cnt = pcnt),
	type = "response"), lty = 2, lwd = 2)
points(cnt[8], (ny.mod/(ny.mod + nn.mod))[8], pch = 21, 
	bg = rgb(0.4, 0.4, 0.4, 0.7))   
lines(pcnt, predict(weib3.glm, data.frame(cnt = pcnt),
	type = "response"), lwd = 3, lty = 3)


###################################################
### code chunk number 18: PsychFun2.Stex:604-609
###################################################
data(StairCase)
print(
xyplot(Contrast ~ Trial, StairCase, groups = StairCase, 
	type = c("p", "l"), scale = list(y = list(log = TRUE)))
)


###################################################
### code chunk number 19: PsychFun2.Stex:621-622
###################################################
names(StairCase)


###################################################
### code chunk number 20: PsychFun2.Stex:639-641
###################################################
sc.glm <- glm(Response ~ log10(Contrast), 
	binomial(mafc.logit(2)), StairCase)


###################################################
### code chunk number 21: PsychFun2.Stex:648-649
###################################################
sc.mns <- with(StairCase, tapply(Response, Contrast, mean))


###################################################
### code chunk number 22: PsychFun2.Stex:655-665
###################################################
plot(as.numeric(names(sc.mns)), sc.mns, log = "x",
	xlab = "Contrast", ylab = "Proportion Correct")
cnt <- seq(0.03, 0.5, len = 200)
lines(cnt, predict(sc.glm, newdata = data.frame(Contrast = cnt),
	type = "response"))	
sc2.glm <- psyfun.2asym(cbind(Response, 1 - Response) ~ 
	log10(Contrast), data = StairCase, link = logit.2asym, 
	init.g = 0.5)
lines(cnt, predict(sc2.glm, newdata = data.frame(Contrast = cnt),
	type = "response"), lty = 2)	


###################################################
### code chunk number 23: PsychFun2.Stex:713-727
###################################################
opar <- par(xpd = TRUE)
xx <- seq(-3, 3, len = 100)
plot(xx, pnorm(xx, 0.5), type = "l", lwd = 2,
	xlab = "Stimulus Level", ylab = "Proportion Detected")
q75 <- qnorm(0.75, 0.5)
q5 <- qnorm(0.5, 0.5)
segments(c(-3, q75), c(0.75, 0.75), c(q75, q75), c(0.75, 0), lty = 3, lwd = 2)
segments(c(-3, q5), c(0.5, 0.5), c(q5, q5), c(0.5, 0), lty = 3, lwd = 2)
arrows(q5, 0.05, q5, -0.0325)
arrows(q75, 0.05, q75, -0.0325)
text(-2.5, 0.5, "p = 0.5", pos = 3)
text(-2.5, 0.75, "p = 0.75", pos = 3)
arrows(q5, -0.1, q75, -0.1, angle = 90, code = 3, length = 0.1)
par(opar)


###################################################
### code chunk number 24: PsychFun2.Stex:813-824
###################################################
ln10 <- log(10)
palpha <- 1 - exp(-1)
gam <- weib3.glm$gam
lam <- weib3.glm$lambda
cc <- coef(weib3.glm)
P <-  (palpha - gam)/(1 - gam - lam)
thresh_alpha <- qweibull( P,  shape = cc[2]/ln10,  
	scale = exp(-ln10 * cc[1]/cc[2]) )
weibfit.params <- c(thresh_alpha, cc[2]/ln10)
names(weibfit.params) <- c("alpha", "beta")
weibfit.params


###################################################
### code chunk number 25: PsychFun2.Stex:844-845
###################################################
coef(summary(weib.glm))


###################################################
### code chunk number 26: PsychFun2.Stex:851-852
###################################################
sqrt(diag(vcov(weib.glm)))


###################################################
### code chunk number 27: PsychFun2.Stex:860-862
###################################################
sqrt(c(1, log10(0.026)) %*% vcov(weib.glm) %*% 
	c(1, log10(0.026)))


###################################################
### code chunk number 28: PsychFun2.Stex:867-869
###################################################
predict(weib.glm, newdata = list(cnt = 0.026), 
	se.fit = TRUE)$se.fit


###################################################
### code chunk number 29: PsychFun2.Stex:884-896
###################################################
cnt <- 10^seq(-2, -1, length = 8) 
plot(cnt, p.obs, log = "x", type = "n",  ylim = c(0, 1),
	ylab  = "Proportion Correct", xlab = "Stimulus Intensity")
pcnt <- seq(0.01, 0.1, len = 100)
weib.pred <- predict(weib.glm, newdata = list(cnt = pcnt), 
	type = "response", se.fit = TRUE)
polygon(c(pcnt, rev(pcnt)), 
	with(weib.pred, c(fit + 2 * se.fit, rev(fit - 2 * se.fit))),
	col = "grey", border = "white" )
lines(pcnt, weib.pred$fit, lwd = 2)
points(cnt, (ny/(ny + nn)), pch = 21, 
	bg = "white")  


###################################################
### code chunk number 30: PsychFun2.Stex:910-911
###################################################
confint(weib.glm, level = 0.99)


###################################################
### code chunk number 31: PsychFun2.Stex:947-949
###################################################
rbinom(nrow(weib.glm$data), weib.glm$prior.weights,
   fitted(weib.glm))


###################################################
### code chunk number 32: PsychFun2.Stex:959-973
###################################################
psyfun.boot <- function(obj, N = 100){
   n <- obj$prior.weights
   f <- fitted(obj)
   resp.bt <- matrix(rbinom(N * length(n), n, f), ncol = N)
   bt.res <- sapply(seq_len(N), function(x) {
      r <- resp.bt[, x]
      res.bt <- glm(cbind(r, n - r) ~ model.matrix(obj) - 1,
         binomial(obj$family$link))
      cc <- coef(res.bt)
      names(cc) <- names(coef(obj))
   cc
    })
    bt.res
}


###################################################
### code chunk number 33: PsychFun2.Stex:1003-1010
###################################################
psyfun.stat <- function(d){
   nn <- rowSums(d[, 1:2])
   mm <- model.matrix(as.matrix(d[, 1:2]) ~ log10(d[, 3]))
   t.glm <- glm(cbind(d[, 4], nn - d[, 4]) ~ mm - 1,
      binomial("cloglog"))
      as.vector(coef(t.glm))
}


###################################################
### code chunk number 34: PsychFun2.Stex:1015-1020
###################################################
psyfun.gen <- function(d, mle) {
   nn <- with(d, ny + nn)
   d$resp <- rbinom(nrow(d), nn, mle)
   d	
}


###################################################
### code chunk number 35: PsychFun2.Stex:1033-1039
###################################################
library(boot)
weib.d <- cbind(weib.sim, resp = weib.sim$ny)
weib.boot <- boot(weib.d, statistic = psyfun.stat,
   R = 1000, sim = "parametric",
   ran.gen = psyfun.gen,
   mle = fitted(weib.glm))


###################################################
### code chunk number 36: PsychFun2.Stex:1046-1047
###################################################
weib.boot


###################################################
### code chunk number 37: PsychFun2.Stex:1087-1089 (eval = FALSE)
###################################################
## plot(weib.boot, index = 1)
## plot(weib.boot, index = 2)


###################################################
### code chunk number 38: PsychFun2.Stex:1099-1103 (eval = FALSE)
###################################################
## alpha <- with(weib.boot, exp(-log(10) * t[, 1]/t[, 2]))
## cc <- coef(weib.glm)
## plot(weib.boot, index = 1, t0 = exp(-log(10) * cc[1]/cc[2]), 
##    t = alpha)


###################################################
### code chunk number 39: PsychFun2.Stex:1110-1111
###################################################
plot(weib.boot, index = 1)


###################################################
### code chunk number 40: PsychFun2.Stex:1114-1115
###################################################
plot(weib.boot, index = 2)


###################################################
### code chunk number 41: PsychFun2.Stex:1118-1121
###################################################
alpha <- with(weib.boot, exp(-log(10) * t[, 1]/t[, 2]))
cc <- coef(weib.glm)
plot(weib.boot, index = 1, t0 = exp(-log(10) * cc[1]/cc[2]), t = alpha)


###################################################
### code chunk number 42: PsychFun2.Stex:1141-1145
###################################################
boot.ci(weib.boot, type = c("norm", "basic", "perc", "bca"), 
   L = influence(weib.glm)$hat, index = 1)
boot.ci(weib.boot, type = c("norm", "basic", "perc", "bca"), 
   L = influence(weib.glm)$hat, index = 2)


###################################################
### code chunk number 43: PsychFun2.Stex:1157-1160
###################################################
boot.ci(weib.boot, type = c("norm", "basic", "perc", "bca"), 
   L = influence(weib.glm)$hat, 
   t = alpha, t0 = exp(-log(10) * cc[1]/cc[2]))


###################################################
### code chunk number 44: PsychFun2.Stex:1203-1215
###################################################
cxx <- seq(0.01, 0.1, len = 100)
cc.bt <- weib.boot$t
l10 <- log(10)
with(weib.sim, plot(cnt, ny/(ny + nn), type = "n", 
   xlim = c(0.01, 0.1), ylim = c(0, 1), log = "x",
   xlab = "Contrast", ylab = "P(Correct)"))
invisible(apply(cc.bt, 1, function(x) {
   xx <- c(exp(-l10 * x[1]/x[2]), x[2]/l10)
   lines(cxx, pweibull(cxx, xx[2], xx[1]), 
      col = rgb(0, 0, 0, 0.005) )
	}))
points(I(ny/(ny + nn)) ~ cnt, weib.sim, pch = 21, bg = "white", cex = 1.5)


###################################################
### code chunk number 45: PsychFun2.Stex:1250-1252
###################################################
library(modelfree)
library(mgcv)


###################################################
### code chunk number 46: PsychFun2.Stex:1262-1264
###################################################
Xie.df <- example04 
Xie.df$p <- with(Xie.df, r/m)


###################################################
### code chunk number 47: PsychFun2.Stex:1304-1313
###################################################
xx <- seq(0, 10, len = 200)
opt.bw <- with(Xie.df, 
   bandwidth_cross_validation(r, m, x, 10^c(-1, 1)))
Xie.mf <- with(Xie.df, 
   locglmfit(xx, r, m, x, opt.bw$deviance))
Xie.gam <- gam(cbind(r, m - r) ~ s(x), binomial,
   Xie.df)	
Xie.glm <- psyfun.2asym(cbind(r, m - r) ~ x, Xie.df,
   init.g = 0.5, init.lam = 0.1, mxNumAlt = 100)


###################################################
### code chunk number 48: PsychFun2.Stex:1323-1335
###################################################
plot(xx, Xie.mf$pfit, type = "l", lty = 2,
	ylim = c(0.5, 1), lwd = 2,
	xlab = expression(paste(plain(Log)[2], 
	  " Separation (", plain(log)[2], " pixels)")), 
	ylab = "Proportion Correct")
lines(xx, predict(Xie.gam, newdata = data.frame(x = xx), 
	type = "response"), lwd = 2, lty = 1)	
lines(xx, predict(Xie.glm, newdata = data.frame(x = xx), 
	type = "response"), lwd = 3, lty = 3)
points(p ~ x, Xie.df, pch = 21, bg = "white")
legend(0, 0.8, legend = c("gam", "local", "glm"), 
	lty = 1:3, lwd = c(2, 2, 4), bty = "n")


###################################################
### code chunk number 49: PsychFun2.Stex:1353-1365 (eval = FALSE)
###################################################
## plot(xx, Xie.mf$pfit, type = "l", lty = 2,
## 	ylim = c(0.5, 1), lwd = 2,
## 	xlab = expression(paste(plain(Log)[2], 
## 	  " Separation (", plain(log)[2], " pixels)")), 
## 	ylab = "Proportion Correct")
## lines(xx, predict(Xie.gam, newdata = data.frame(x = xx), 
## 	type = "response"), lwd = 2, lty = 1)	
## lines(xx, predict(Xie.glm, newdata = data.frame(x = xx), 
## 	type = "response"), lwd = 3, lty = 3)
## points(p ~ x, Xie.df, pch = 21, bg = "white")
## legend(0, 0.8, legend = c("gam", "local", "glm"), 
## 	lty = 1:3, lwd = c(2, 2, 4), bty = "n")



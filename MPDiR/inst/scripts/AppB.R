### R code from vignette source 'StatBack/StatBack.Stex'

###################################################
### code chunk number 1: StatBack.Stex:37-56
###################################################
options(width=65, digits = 3)
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
### code chunk number 2: StatBack.Stex:105-120
###################################################
opar <- par(mfrow = c(1, 2))
RT <- scan("florRT.txt")
RT.fit <- MASS:::fitdistr(RT, "gamma")
set.seed(121952)
RT.sim <- rgamma(length(RT), shape = RT.fit$estimate[1],
	rate = RT.fit$estimate[2])
rng <- round(range(c(RT, RT.sim)), 1) + c(-0.1, 0.1)
brks <- seq(rng[1], rng[2], len = 9)
hist(RT, breaks = brks, main = "Reaction Times",
	xlab = "Time (sec)") 
mtext("a", adj = 0, line = 1, cex = 1.5)
hist(RT.sim, breaks = brks, main = "Gamma Random Variates",
	xlab = "Time (sec)") 
mtext("b", adj = 0, line = 1, cex = 1.5)
par(opar)


###################################################
### code chunk number 3: StatBack.Stex:261-268
###################################################
opar <- par(xpd = NA)
x <- seq(0, 12, len = 200)
y <-  dgamma(x, shape = 4)
plot(x, y, type = "l", cex.axis = 1.5, cex.lab = 2)
polygon(x[ c(120, 120:145, 145) ], c(0, y[120:145], 0), border = "black", col = "grey")
text(x[c(120, 145)], c(-0.015, -0.015), letters[1:2], cex = 1.3)
par(opar)


###################################################
### code chunk number 4: StatBack.Stex:345-346 (eval = FALSE)
###################################################
## help("rnorm")


###################################################
### code chunk number 5: StatBack.Stex:357-364
###################################################
opar <- par(mfrow = c(1, 2), cex.lab=2, cex.axis = 1.5, mar = c(5, 5, 4, 2  + 0.1))
x <- seq(-3, 3, len = 100)
plot(x, dnorm(x), type = "l")
mtext("a", 3, adj = 0, cex = 2, line = 1)
plot(x, pnorm(x), type = "l")
mtext("b", 3, adj = 0, cex = 2, line = 1)
par(opar)


###################################################
### code chunk number 6: StatBack.Stex:372-373
###################################################
qnorm(pnorm(1))


###################################################
### code chunk number 7: StatBack.Stex:376-377
###################################################
pnorm(qnorm(0.8413447))


###################################################
### code chunk number 8: StatBack.Stex:415-426
###################################################
opar <- par(mfrow = c(1, 2), cex.lab=2, cex.axis = 1.5, mar = c(5, 5, 4, 2  + 0.1))
x <- seq(-3, 3, len = 100)
plot(x, dnorm(x), type = "l")
lines(x, dnorm(x, mean = 1))
lines(x, dnorm(x, mean = 2))
mtext("a", 3, adj = 0, cex = 2, line = 1)
plot(x, dnorm(x), type = "l")
lines(x, dnorm(x, sd = 2))
lines(x, dnorm(x, sd = 3))
mtext("b", 3, adj = 0, cex = 2, line = 1)
par(opar)


###################################################
### code chunk number 9: StatBack.Stex:443-454
###################################################
opar <- par(mfrow = c(1, 2), cex.lab=2, cex.axis = 1.5, mar = c(5, 5, 4, 2  + 0.1))
x <- seq(-3, 3, len = 100)
plot(x, pnorm(x), type = "l")
lines(x, pnorm(x, mean = 1))
lines(x, pnorm(x, mean = 2))
mtext("a", 3, adj = 0, cex = 2, line = 1)
plot(x, pnorm(x), type = "l")
lines(x, pnorm(x, sd = 2))
lines(x, pnorm(x, sd = 3))
mtext("b", 3, adj = 0, cex = 2, line = 1)
par(opar)


###################################################
### code chunk number 10: StatBack.Stex:500-501
###################################################
(X <- rnorm(10,mean=0,sd=1))


###################################################
### code chunk number 11: StatBack.Stex:505-506
###################################################
(Y <- rnorm(10,mean=10,sd=1))


###################################################
### code chunk number 12: StatBack.Stex:646-653
###################################################
opar <- par(mfrow = c(1, 2), cex.lab=2, cex.axis = 1.5, mar = c(5, 5, 4, 2  + 0.1), lwd = 2)
x <- c(seq(-1, 2, len = 1000))
plot(x, dunif(x), type = "l", ylim = c(0, 1.2), xlim = c(-0.5, 1.5))  
mtext("a", 3, adj = 0, cex = 2, line = 1)
plot(x, punif(x), type = "l", ylim = c(0, 1.2), xlim = c(-0.5, 1.5))   
mtext("b", 3, adj = 0, cex = 2, line = 1)
par(opar)  


###################################################
### code chunk number 13: StatBack.Stex:732-743
###################################################
opar <- par(cex.lab = 2, cex.axis = 1.5, mar = c(5, 5, 4, 2) + 0.1)
x <- seq(0, 10, len = 100)
plot(x, dgamma(x, 1), type = "l")
lines(x, dgamma(x,2), type="l")
lines(x, dgamma(x,3), type="l")
lines(x, dgamma(x,4), type="l")
lines(x, dgamma(x,5), type="l")
x <- c(0.5,1.4,2.5,3.5,4.5)
y <- c(0.8,0.38,0.3,0.26, 0.23)
text(x,y,1:5, cex = 1.25)
par(opar)


###################################################
### code chunk number 14: StatBack.Stex:1041-1056
###################################################
opar <- par(lwd = 2)
n <- 3
x <- seq(-3, 3, len = 100)
plot(c(-3,3),c(-5,0), type = "n",
xlab = expression(paste(mu, " (Gaussian mean)")),
ylab = "Log Likelihood")   
s <- c(-0.5, 1, 1.5) ## The 'sample' (let's take no chances :^) % Ha Ha. You don't REALLY believe in statistics do you ...?
ysum <- 0
for (i in 1:n) {
  y <- -(x-s[i])^2
  ysum <- y + ysum
  lines(x,y)
  }
lines(x,ysum,lwd=5)
par(opar)


###################################################
### code chunk number 15: StatBack.Stex:1176-1200
###################################################
samp <- c(549,550,551,551,552,552,552,552,552,553,553,553,554,554,554,554,
	555,555,555,556,556,556,556,556,556,557,557,557,557,557,558,558,559,
	559,560, 561,561,561,561,561,561,562,562,562,563,563,564,564,564,564,
	565,565,565,566,567,568,569,570)

hist(samp, breaks = diff(range(samp)) + 1, #freq = FALSE,
	xlab = expression(paste(lambda[max], " (L cones)")),
	main = "", # "Dartnall et al. (1983)",
	cex.lab = 1.5, cex.axis = 1.3)
#samp.dens <- density(samp, bw = 1.75)
#lines(samp.dens$x, samp.dens$y)
mix.obj <- function(p, x)
{
  e <- p[1] * dnorm((x - p[2])/p[3])/p[3] +
       (1 - p[1]) * dnorm((x - p[4])/p[5])/p[5]
  if(any(e <= 0)) Inf else -sum(log(e))
}
e <- function(p, x) p[1] * dnorm((x - p[2])/p[3])/p[3] +
       (1 - p[1]) * dnorm((x - p[4])/p[5])/p[5]
p0 <- c(0.6, 555, 3, 562, 3)
p <- optim(p0, mix.obj, x = samp, 
		control = list(maxit = 10000))
xx <- seq(545, 575, len = 200)	
lines(xx, sum(table(samp)) * e(p$par, xx), lwd = 3)


###################################################
### code chunk number 16: StatBack.Stex:1309-1320
###################################################
p <- seq(0.01, 0.99, len = 99)
opar <- par(lwd = 2, cex.axis = 1.5, cex.lab = 2)
nH = 3; nT = 7
y <- dbeta(p, nH + 1, nT + 1)
plot(p, y,
	 type = "n",
	 xlab = "P",
	 ylab = "posterior",
	 cex.lab = 1.8, cex.axis = 1.35)
lines(p, y, lwd = 2)
par(opar)


###################################################
### code chunk number 17: StatBack.Stex:1336-1357
###################################################
p <- seq(0.01, 0.99, len = 99)
opar <- par(lwd = 2, cex.axis = 1.5, cex.lab = 2)
nH = 3; nT = 3;
y1 <- dbeta(p,nH+1,nT+1)
nH = 5; nT = 5;
y2 <- dbeta(p,nH+1,nT+1)
nH = 10; nT = 10;
y3 <- dbeta(p,nH+1,nT+1)
nH = 8; nT = 12;
yp <- dbeta(p,nH+1,nT+1);
plot(p,yp,
	 type = "n",
	 xlab = "P",
	 ylab = "probability density",
	 cex.lab = 1.8, cex.axis = 1.35)
lines(p,y1,lwd=2,lty=3)
lines(p,y2,lwd=2,lty=2)
lines(p,y3,lwd=2,lty=3)

lines(p,yp,lwd=2,lty=1)
par(opar)


###################################################
### code chunk number 18: StatBack.Stex:1472-1487
###################################################

opar <- par(lwd = 2, cex.axis = 1.5, cex.lab = 2, mar = c(5, 6, 4, 2) + 0.1,mfrow  = c(1, 1))
n <- 100
smp <-  matrix(rexp(n * 10000), nrow = n)
xbar <- colMeans(smp)
lambdaMLE <- 1/xbar    # MLE estimate
xsum = n*xbar          # for convenience
logL0 <- -xsum                                 # lambda = 1 (constrained)
logL1 <- -lambdaMLE*xsum + n*log( lambdaMLE )  # lambda MLE
w <-  2 * (logL1 - logL0 )
hist(w, breaks = 40, freq = FALSE, main = "Sample Size = 100")  
xx <- seq(0, max(w), len = 200)
lines(xx, dchisq(xx, 1), lwd = 3)

par(opar)



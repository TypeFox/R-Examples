### R code from vignette source 'article.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: article.Rnw:1-7
###################################################
# This should be set to FALSE in case you want to redo the sampling as well.
# Default is TRUE to save CRAN some time.
usePreCalcResults <- TRUE

# Default here is FALSE because the object created is rather large
usePreCalcResultsExtra <- FALSE 


###################################################
### code chunk number 2: preliminaries
###################################################
options(prompt = 'R> ', continue = '+  ')


###################################################
### code chunk number 3: usd1
###################################################
set.seed(123)
library("stochvol")
data("exrates")
ret <- logret(exrates$USD, demean = TRUE)
par(mfrow = c(2, 1), mar = c(1.9, 1.9, 1.9, 0.5), mgp = c(2, 0.6, 0))
plot(exrates$date, exrates$USD, type = "l",
  main = "Price of 1 EUR in USD")
plot(exrates$date[-1], ret, type = "l", main = "Demeaned log returns")


###################################################
### code chunk number 4: usd2
###################################################
sim <- svsim(500, mu = -9, phi = 0.99, sigma = 0.1)
par(mfrow = c(2, 1))
plot(sim)


###################################################
### code chunk number 5: article.Rnw:239-241 (eval = FALSE)
###################################################
## res <- svsample(ret, priormu = c(-10, 1), priorphi = c(20, 1.1),
##   priorsigma = 0.1)


###################################################
### code chunk number 6: article.Rnw:244-250
###################################################
if (usePreCalcResultsExtra) {
 load("vignette_sampling_draws1.RData")
} else {
 res <- svsample(ret, priormu = c(-10, 1), priorphi = c(20, 1.1), priorsigma = .1)
 #save(res, file = "vignette_sampling_draws1.RData")
}


###################################################
### code chunk number 7: article.Rnw:272-273
###################################################
summary(res, showlatent = FALSE)


###################################################
### code chunk number 8: usd3
###################################################
volplot(res, forecast = 100, dates = exrates$date[-1])


###################################################
### code chunk number 9: usd4
###################################################
res <- updatesummary(res, quantiles = c(0.01, 0.1, 0.5, 0.9, 0.99))
volplot(res, forecast = 100, dates = exrates$date[-1])


###################################################
### code chunk number 10: usd5
###################################################
par(mfrow = c(3, 1))
paratraceplot(res)


###################################################
### code chunk number 11: usd6
###################################################
par(mfrow = c(1, 3))
paradensplot(res, showobs = FALSE)


###################################################
### code chunk number 12: usd7
###################################################
plot(res, showobs = FALSE)


###################################################
### code chunk number 13: usd8
###################################################
myresid <- resid(res)
plot(myresid, ret)


###################################################
### code chunk number 14: article.Rnw:404-410
###################################################
set.seed(123456)
n <- 1000
beta.true <- c(0.1, 0.5)
sigma.true <- 0.01
X <- matrix(c(rep(1, n), rnorm(n, sd = sigma.true)), nrow = n)
y <- rnorm(n, X %*% beta.true, sigma.true)


###################################################
### code chunk number 15: article.Rnw:414-420
###################################################
burnin <- 100
draws <- 5000
b0 <- matrix(c(0, 0), nrow = ncol(X))
B0inv <- diag(c(10^-10, 10^-10))
c0 <- 0.001
C0 <- 0.001


###################################################
### code chunk number 16: article.Rnw:423-427
###################################################
p <- ncol(X)
preCov <- solve(crossprod(X) + B0inv)
preMean <- preCov %*% (crossprod(X, y) + B0inv %*% b0)
preDf <- c0 + n/2 + p/2


###################################################
### code chunk number 17: article.Rnw:430-433
###################################################
draws1 <- matrix(NA_real_, nrow = draws, ncol = p + 1)
colnames(draws1) <- c(paste("beta", 0:(p-1), sep = "_"), "sigma")
sigma2draw <- 1


###################################################
### code chunk number 18: article.Rnw:436-444
###################################################
for (i in -(burnin-1):draws) {
  betadraw <- as.numeric(mvtnorm::rmvnorm(1, preMean,
    sigma2draw * preCov))
  tmp <- C0 + 0.5 * (crossprod(y - X %*% betadraw) +
    crossprod((betadraw - b0), B0inv) %*% (betadraw - b0))
  sigma2draw <- 1 / rgamma(1, preDf, rate = tmp)
  if (i > 0) draws1[i, ] <- c(betadraw, sqrt(sigma2draw))
}


###################################################
### code chunk number 19: article.Rnw:447-448
###################################################
colMeans(draws1)


###################################################
### code chunk number 20: homo
###################################################
par(mar = c(3.1, 1.8, 1.9, .5), mgp = c(1.8, .6, 0))
plot(coda::mcmc(draws1), show.obs = FALSE)


###################################################
### code chunk number 21: article.Rnw:458-459 (eval = FALSE)
###################################################
## plot(coda::mcmc(draws1), show.obs = FALSE)


###################################################
### code chunk number 22: article.Rnw:481-486
###################################################
mu.true <- log(sigma.true^2)
phi.true <- 0.97
vv.true <- 0.3
simresid <- svsim(n, mu = mu.true, phi = phi.true, sigma = vv.true)
y <- X %*% beta.true + simresid$y


###################################################
### code chunk number 23: article.Rnw:489-495
###################################################
draws <- 50000
burnin <- 1000
thinning <- 10
priormu <- c(-10, 2)
priorphi <- c(20, 1.5)
priorsigma <- 1


###################################################
### code chunk number 24: article.Rnw:498-505
###################################################
draws2 <- matrix(NA_real_, nrow = floor(draws / thinning),
  ncol = 3 + n + p)
colnames(draws2) <- c("mu", "phi", "sigma",
  paste("beta", 0:(p-1), sep = "_"), paste("h", 1:n, sep = "_"))
betadraw <- c(0, 0)
svdraw <- list(para = c(mu = -10, phi = 0.9, sigma = 0.2),
  latent = rep(-10, n))


###################################################
### code chunk number 25: article.Rnw:512-529 (eval = FALSE)
###################################################
## for (i in -(burnin-1):draws) {
##   ytilde <- y - X %*% betadraw
##   svdraw <- svsample2(ytilde, startpara = para(svdraw),
##     startlatent = latent(svdraw), priormu = priormu,
##     priorphi = priorphi, priorsigma = priorsigma)
##   normalizer <- as.numeric(exp(-latent(svdraw) / 2))
##   Xnew <- X * normalizer
##   ynew <- y * normalizer
##   Sigma <- solve(crossprod(Xnew) + B0inv)
##   mu <- Sigma %*% (crossprod(Xnew, ynew) + B0inv %*% b0)
##   betadraw <- as.numeric(mvtnorm::rmvnorm(1, mu, Sigma))
##   if (i > 0 & i %% thinning == 0) { 
##     draws2[i/thinning, 1:3] <- para(svdraw)
##     draws2[i/thinning, 4:5] <- betadraw
##     draws2[i/thinning, 6:(n+5)] <- latent(svdraw)
##   }
## }


###################################################
### code chunk number 26: article.Rnw:532-561
###################################################
if (usePreCalcResults) {
 load("vignette_sampling_draws2.RData")
} else {
for (i in -(burnin-1):draws) {

 # draw latent volatilities and AR-parameters:
 ytilde <- y - X %*% betadraw
 svdraw <- svsample2(ytilde, startpara = para(svdraw),
                     startlatent = latent(svdraw), priormu = priormu,
                     priorphi = priorphi, priorsigma = priorsigma)

 # draw the betas:
 normalizer <- as.numeric(exp(-latent(svdraw)/2))
 Xnew <- X * normalizer
 ynew <- y * normalizer
 Sigma <- solve(crossprod(Xnew) + B0inv)
 mu <- Sigma %*% (crossprod(Xnew, ynew) + B0inv %*% b0)
 betadraw <- as.numeric(mvtnorm::rmvnorm(1, mu, Sigma))

 # store the results:
 if (i > 0 & i %% thinning == 0) { 
  draws2[i/thinning, 1:3] <- para(svdraw)
  draws2[i/thinning, 4:5] <- betadraw
  draws2[i/thinning, 6:(n+5)] <- latent(svdraw)
 }
}
 draws2selection <- draws2[,4:8]
 #save(draws2selection, file = 'vignette_sampling_draws2.RData')
}


###################################################
### code chunk number 27: article.Rnw:565-566 (eval = FALSE)
###################################################
## plot(coda::mcmc(draws2[, 4:7]), show.obs = FALSE)


###################################################
### code chunk number 28: hetero
###################################################
par(mar = c(3.1, 1.8, 1.9, .5), mgp = c(1.8, .6, 0))
plot(coda::mcmc(draws2selection[, 1:4]), show.obs = FALSE)


###################################################
### code chunk number 29: article.Rnw:584-585 (eval = FALSE)
###################################################
## colMeans(draws2[, 4:8])


###################################################
### code chunk number 30: article.Rnw:588-589
###################################################
colMeans(draws2selection)


###################################################
### code chunk number 31: scatter
###################################################
x <- log(exrates$USD[-length(exrates$USD)])
y <- log(exrates$USD[-1])
X <- matrix(c(rep(1, length(x)), x), nrow = length(x))
par(mfrow=c(1,1), mar = c(2.9, 2.9, 2.7, .5), mgp = c(1.8,.6,0), tcl = -.4)
plot(x,y, xlab=expression(log(p[t])), ylab=expression(log(p[t+1])),
    main="Scatterplot of lagged daily log prices of 1 EUR in USD",
    col="#00000022", pch=16, cex=1)
abline(0,1)


###################################################
### code chunk number 32: article.Rnw:654-875
###################################################
if (usePreCalcResults) {
 load("vignette_sampling_realworld.RData")
} else {
set.seed(34567890)

#configuration parameters
draws <- 100000
thinning <- 100
burnin <- 10000
priormu <- c(-10, 1)
priorphi <- c(20, 1.5)
priorsigma <- .1
p <- 2

n <- length(x)

## design matrix:
X <- matrix(c(rep(1, n), x), nrow=n)

## prior for beta (or beta|sigma in homoskedastic regression):
b0 <- matrix(c(0, 0), nrow=ncol(X))
B0 <- diag(c(10^10, 10^10))

## prior for sigma^2 (only used in homoskedastic regression)
c0 <- 0.001
C0 <- 0.001

B0inv <- solve(B0)

## initialize some space to hold results:
realres <- vector("list", 3)

## AR(1)-SV:
realres[[1]] <- matrix(NA_real_, nrow = floor(draws/thinning), ncol = 3 + n + p)
colnames(realres[[1]]) <- c("mu", "phi", "sigma", paste("beta", 0:(p-1), sep="_"),
			   paste("h", 1:n, sep="_"))

## some indicators:
paras <- 1:3
betas <- 3+(1:p)
latents <- 3+p+(1:n)

## starting values:
betadraw <- rep(.1, p)
svdraw <- list(para = c(mu = -10, phi = .8, sigma = .2),
	       latent = rep(-10, n))

## sampler:
tim <- system.time(
for (i in -(burnin-1):draws) {
 if (i%%1000 == 0) cat("Iteration", i, "done.\n")
 # draw latent volatilities and SV-parameters:
 ytilde <- y - X %*% betadraw
 svdraw <- svsample2(ytilde, startpara=para(svdraw),
		      startlatent=latent(svdraw), priormu=priormu,
		      priorphi=priorphi, priorsigma=priorsigma)

 # draw the betas:
 normalizer <- as.numeric(exp(-latent(svdraw)/2))
 Xnew <- X * normalizer
 ynew <- y * normalizer
 Sigma <- solve(crossprod(Xnew) + B0inv)
 mu <- Sigma %*% (crossprod(Xnew, ynew) + B0inv %*% b0)
 betadraw <- as.numeric(mvtnorm::rmvnorm(1, mu, Sigma))

 # store the results:
 if (i > 0 & i %% thinning == 0) { 
  realres[[1]][i/thinning, paras] <- para(svdraw)
  realres[[1]][i/thinning, latents] <- latent(svdraw)
  realres[[1]][i/thinning, betas] <- betadraw
 }
}
)[["elapsed"]]

## AR(1) homoskedastic:
## pre-calculate some values:
preCov <- solve(crossprod(X) + B0inv)
preMean <- preCov %*% (crossprod(X, y) + B0inv %*% b0)
preDf <- c0 + n/2 + p/2

## allocate space:
realres[[2]] <- matrix(as.numeric(NA), nrow=floor(draws/thinning), ncol=p+1)
colnames(realres[[2]]) <- c(paste("beta", 0:(p-1), sep="_"), "sigma")

## starting values:
betadraw <- rep(0, p)
sigma2draw <- var(y)

## sampler:
tim2 <- system.time(
for (i in -(burnin-1):draws) {
 if (i%%1000 == 0) cat("Iteration", i, "done.\n")
 # draw beta:
 betadraw <- as.numeric(mvtnorm::rmvnorm(1, preMean, sigma2draw*preCov))
 
 # draw sigma^2:
 tmp <- C0 + .5*(crossprod(y - X%*%betadraw) +
		crossprod((betadraw-b0), B0inv) %*% (betadraw-b0))
 sigma2draw <- 1/rgamma(1, preDf, rate=tmp)

 # store the results:
 if (i > 0 & i %% thinning == 0) { 
  realres[[2]][i/thinning, 1:p] <- betadraw
  realres[[2]][i/thinning, p+1] <- sqrt(sigma2draw)
 }
}
)[["elapsed"]]

## AR(1)-GARCH(1,1):
## allocate space:
realres[[3]] <- matrix(NA_real_, nrow = floor(draws/thinning), ncol = 3 + n + p)
colnames(realres[[3]]) <- c("a_0", "a_1", "b_1", paste("beta", 0:(p-1), sep="_"),
			   paste("sigma", 1:n, sep="_"))

## some auxiliary functions:
loglik <- function(y, s2) {
 sum(dnorm(y, 0, sqrt(s2), log=TRUE))
}

s2s <- function(y, para, y20, s20) {
 s2 <- rep(NA_real_, length(y))
 s2[1] <- para[1] + para[2]*y20 + para[3]*s20
 for (i in 2:length(y)) s2[i] <- para[1] + para[2]*y[i-1]^2 + para[3]*s2[i-1]
 s2	
}

## This is a very simple (and slow) random-walk MH sampler
## which can most certainly be improved...

## initial values:

# crude approximation of residual variance:
tmplm <- lm.fit(X, y)
s20 <- var(resid(tmplm))

y20 <- 0

## starting values:
## For the starting values for the parameters,
## use ML results from fGarch::garchFit on the
## residuals

# tmpfit <- fGarch::garchFit(~garch(1,1), resid(tmplm), trace=FALSE)
# startvals <- fGarch::coef(tmpfit)[c("omega", "alpha1", "beta1")]
startvals <- c(1.532868e-07, 3.246417e-02, 9.644149e-01)

## Random walk MH algorithm needs manual tuning for good mixing.
## Gelman et al. (1996) suggest an acceptance rate of 31.6%
## for spherical 3d multivariate normal (which we don't have :))
## These values seem to work OK:
# mhtuning <- tmpfit@fit$se.coef[c("omega", "alpha1", "beta1")] / 2
mhtuning <- c(3.463851e-08, 2.226946e-03, 2.354829e-03)

betadraw <- as.numeric(coef(tmplm))

paradraw <- startvals
paraprop <- rep(NA_real_, 3)
 
s2draw <- rep(.1, length(y))
s2prop <- rep(NA_real_, length(y))

accepts <- 0

## sampler:
tim3 <- system.time(
for (i in -(burnin-1):draws) {
 
 if (i%%1000 == 0) cat("Iteration", i, "done.\n")
 
 # draw volatilities and GARCH-parameters:
 ytilde <- y - X %*% betadraw  # regression residuals

 paraprop <- rnorm(3, paradraw, mhtuning)
 s2prop <- s2s(ytilde, paraprop, y20, s20)

 if (all(s2prop > 0)) { # s2 needs to be positive, otherwise reject
  logR <- loglik(ytilde, s2prop) - loglik(ytilde, s2draw)  # log acceptance rate
  if (log(runif(1)) < logR) {
   paradraw <- paraprop
   s2draw <- s2prop
   if (i > 0) accepts <- accepts + 1
  }
 }

 # draw the betas:
 normalizer <- 1/sqrt(s2draw)
 Xnew <- X * normalizer
 ynew <- y * normalizer
 Sigma <- solve(crossprod(Xnew) + B0inv)
 mu <- Sigma %*% (crossprod(Xnew, ynew) + B0inv %*% b0)
 betadraw <- as.numeric(mvtnorm::rmvnorm(1, mu, Sigma))

 # store the results:
 if (i > 0 & i %% thinning == 0) { 
  realres[[3]][i/thinning, paras] <- paradraw
  realres[[3]][i/thinning, latents] <- sqrt(s2draw)
  realres[[3]][i/thinning, betas] <- betadraw
 }
}
)[["elapsed"]]


#Standardized residuals for model checking:
#homoskedastic
predmeans2 <- tcrossprod(X, realres[[2]][,1:2])
standresidmeans2 <- rowMeans((y - predmeans2)/realres[[2]][,3])

#heteroskedastic SV
predmeans <- tcrossprod(X, realres[[1]][,4:5])
standresidmeans <- rowMeans((y - predmeans)/exp(t(realres[[1]][,6:ncol(realres[[1]])])/2))

#heteroskedastic GARCH
predmeans3 <- tcrossprod(X, realres[[3]][,4:5])
standresidmeans3 <- rowMeans((y - predmeans3)/t(realres[[3]][,6:ncol(realres[[3]])]))

realresselection <- vector("list", 3)
realresselection[[1]] <- realres[[1]][,c('beta_0', 'beta_1')]
realresselection[[2]] <- realres[[2]][,c('beta_0', 'beta_1')]
realresselection[[3]] <- realres[[3]][,c('beta_0', 'beta_1')]
#save(realresselection, standresidmeans, standresidmeans2, standresidmeans3, file = 'vignette_sampling_realworld.RData')
}


###################################################
### code chunk number 33: betapost
###################################################
smootherfactor <- 1.8
par(mar = c(2.9, 2.9, 2.7, .5), mgp = c(1.8,.6,0), tcl = -.4)
layout(matrix(c(1,2,3,3), byrow=TRUE, nrow=2))
plot(density(realresselection[[1]][,"beta_0"], bw="SJ", adjust=smootherfactor),
    main=bquote(paste("p(",beta[0],"|",bold(y),")", sep="")),
    xlab="", ylab="", xlim=c(-.0005, .001))
lines(density(realresselection[[3]][,"beta_0"], bw="SJ", adjust=smootherfactor), lty=2, col=4)
lines(density(realresselection[[2]][,"beta_0"], bw="SJ", adjust=smootherfactor), lty=3, col=2)
#points(ols$coefficients[1],0, col=2)
#lines(confint(ols)[1,], rep(0,2), col=2)
legend("topleft", c("SV", "GARCH", "homosked."), col=c(1,4,2), lty=1:3)

plot(density(realresselection[[1]][,"beta_1"], bw="SJ", adjust=smootherfactor),
    main=bquote(paste("p(",beta[1],"|",bold(y),")", sep="")),
    xlab="", ylab="", xlim=c(0.9965, 1.0022))
lines(density(realresselection[[3]][,"beta_1"], bw="SJ", adjust=smootherfactor), lty=2, col=4)
lines(density(realresselection[[2]][,"beta_1"], bw="SJ", adjust=smootherfactor), lty=3, col=2)
#points(ols$coefficients[2],0, col=2)
#lines(confint(ols)[2,], rep(0,2), col=2)
legend("topright", c("SV", "GARCH", "homosked."), col=c(1,4,2), lty=1:3)

plotorder <- sample.int(3*nrow(realresselection[[1]]))
cols <- rep(c("#000000aa", "#0000ff66", "#ff000077"), each=nrow(realresselection[[1]]))[plotorder]
pchs <- rep(1:3, each=nrow(realresselection[[1]]))[plotorder]

beta0 <- c(realresselection[[1]][,"beta_0"],
	   realresselection[[3]][,"beta_0"],
	   realresselection[[2]][,"beta_0"])[plotorder]

beta1 <- c(realresselection[[1]][,"beta_1"],
	   realresselection[[3]][,"beta_1"],
	   realresselection[[2]][,"beta_1"])[plotorder]

plot(beta0, beta1, col=cols, pch=pchs,
    xlab=bquote(paste("p(",beta[0],"|", bold(y), ")")),
    ylab=bquote(paste("p(",beta[1],"|", bold(y), ")")),
    main="Scatterplot of posterior draws")
legend("topright", c("SV", "GARCH", "homosked."), col=c(1,4,2), pch=1:3)


###################################################
### code chunk number 34: qqplot
###################################################
par(mfrow=c(3,2), mar = c(3.1, 3.3, 2.0, .5), mgp = c(1.7,.5,0),
    tcl = -.4)
plot(standresidmeans2, main="Residual scatterplot (homoskedastic errors)",
     xlab="Time", ylab="Standardized residuals")
qqplot(qnorm(ppoints(length(standresidmeans2)), mean=0, sd=1),
       standresidmeans2, main="Residual Q-Q plot (homoskedastic errors)",
       xlab = "Theoretical N(0,1)-quantiles", ylab = "Empirical Quantiles")
abline(0,1)
plot(standresidmeans3, main="Residual scatterplot (GARCH errors)",
     xlab="Time", ylab="Standardized residuals")
qqplot(qnorm(ppoints(length(standresidmeans3)), mean=0, sd=1),
       standresidmeans3, main="Residual Q-Q plot (GARCH errors)",
       xlab = "Theoretical N(0,1)-quantiles", ylab = "Empirical Quantiles")
abline(0,1)
plot(standresidmeans, main="Residual scatterplot (SV errors)",
     xlab="Time", ylab="Standardized residuals")
qqplot(qnorm(ppoints(length(standresidmeans)), mean=0, sd=1),
       standresidmeans, main="Residual Q-Q plot (SV errors)",
       xlab = "Theoretical N(0,1)-quantiles", ylab = "Empirical Quantiles")
abline(0,1)



### R code from vignette source '7-bayes.rnw'

###################################################
### code chunk number 1: Setup
###################################################
options(repos="http://cran.r-project.org")

if(!require(Hmisc, quietly=TRUE)) install.packages("Hmisc")
if(!require(xtable, quietly=TRUE)) install.packages("xtable")

rm( list=ls() )

options(width = 67)


###################################################
### code chunk number 2: 7-bayes.rnw:230-242
###################################################
# Synthetic Poisson
# =====================
# syn.poisson.r
set.seed(1)
nobs <- 50000
x1 <- runif(nobs)
x2 <- runif(nobs)
xb <- 2 + .75*x1 - 1.25*x2 # linear predictor
exb <- exp(xb)             # mean
py <- rpois(nobs, exb)     # random Poisson variates
poireg <-glm(py ~ x1 + x2, family = poisson)
coef(summary(poireg))


###################################################
### code chunk number 3: 7-bayes.rnw:246-247
###################################################
sum(residuals(poireg, type="pearson")^2) / poireg$df.residual


###################################################
### code chunk number 4: 7-bayes.rnw:251-252
###################################################
with(poireg, deviance/df.residual)


###################################################
### code chunk number 5: 7-bayes.rnw:294-308
###################################################
# Synthetic Grouped Logit
# ===============================
# syn.glogit.r
nobs <- 50000
set.seed(1)
x1 <- runif(nobs); x2 <- runif(nobs)
d <- rep(1:5, each=10000, times=1)*100 # binomial denominator
xb <- 2 + .75*x1 - 1.25*x2             # linear predictor
exb <- 1/(1+exp(-xb))                  # mean
by <- rbinom(nobs, size = d, p = exb)  # logit variates
dby = d - by                           # set up y and not-y
gby <- glm(cbind(by, dby) ~ x1 + x2,
           family = binomial)
coef(summary(gby))


###################################################
### code chunk number 6: 7-bayes.rnw:312-313
###################################################
sum(residuals(gby, type="pearson")^2) / gby$df.residual


###################################################
### code chunk number 7: 7-bayes.rnw:317-318
###################################################
with(gby, deviance/df.residual)


###################################################
### code chunk number 8: 7-bayes.rnw:330-331
###################################################
confint(gby)


###################################################
### code chunk number 9: Wald
###################################################
confint.default(gby)  # Wald


###################################################
### code chunk number 10: MC
###################################################
# Monte Carlo Estimation: Synthetic Poisson
# =================================
# sim.poi.r
set.seed(1)
mysim <- function() {
 nobs <- 50000
 x1 <- runif(nobs)
 x2 <- runif(nobs)
 xb <- 2 + .75*x1 -1.25*x2
 exb <- exp(xb)
 py <- rpois(nobs, exb)
 poisim <- glm(py ~ x1 + x2, family=poisson)
 pr <- sum(residuals(poisim, type="pearson")^2)
 prdisp <- pr/poisim$df.residual
 dvdisp <- with(poisim, deviance/df.residual)
 beta   <- poisim$coef
 list(prdisp, dvdisp, beta)
}


###################################################
### code chunk number 11: 7-bayes.rnw:407-409
###################################################
reps <- 1000
B <- replicate(reps, mysim())


###################################################
### code chunk number 12: 7-bayes.rnw:415-416
###################################################
summary(unlist(B[1,]))


###################################################
### code chunk number 13: 7-bayes.rnw:420-421
###################################################
summary(unlist(B[2,]))


###################################################
### code chunk number 14: fig-pois-sim-histo
###################################################
beta.hat <- do.call(rbind, B[3,])
par(mfrow=c(1,3), mar=c(5,4,1,2), las=1)
for (i in 0:2)
  hist(beta.hat[,i+1], breaks = 50, main = "",
       xlab = bquote(paste("Simulated ", beta[.(i)])))


###################################################
### code chunk number 15: fig:pois-sim-histo
###################################################
beta.hat <- do.call(rbind, B[3,])
par(mfrow=c(1,3), mar=c(5,4,1,2), las=1)
for (i in 0:2)
  hist(beta.hat[,i+1], breaks = 50, main = "",
       xlab = bquote(paste("Simulated ", beta[.(i)])))


###################################################
### code chunk number 16: 7-bayes.rnw:534-557
###################################################
# Monte Carlo Estimation: Synthetic Negative Binomial
# ===================================================
# sim.nb2.r
library(MASS)
mysim <- function() {
  nobs <- 50000
  x1 <-runif(nobs)
  x2 <-runif(nobs)
  xb <- 2 + .75*x1 - 1.25*x2
  a <- .5                       # alpha
  ia <- 1/.5                    # theta
  exb <- exp(xb)                # log link
  xg <- rgamma(nobs, a, a, ia)  # gamma variates
  xbg <-exb*xg                  # means
  nby <- rpois(nobs, xbg)       # Poisson variates
  nbsim <-glm.nb(nby ~ x1 + x2) # model
  alpha <- nbsim$theta
  pr <- sum(residuals(nbsim, type="pearson")^2)
  prdisp <- pr/nbsim$df.residual
  dvdisp <- nbsim$deviance/nbsim$df.residual
  beta <- nbsim$coef
  list(alpha, prdisp, dvdisp, beta)
}


###################################################
### code chunk number 17: NB.sim
###################################################
set.seed(1)
reps <- 1000
B <- replicate(reps, mysim())


###################################################
### code chunk number 18: 7-bayes.rnw:571-572
###################################################
summary(unlist(B[1,]))


###################################################
### code chunk number 19: 7-bayes.rnw:576-577
###################################################
summary(unlist(B[2,]))


###################################################
### code chunk number 20: 7-bayes.rnw:581-582
###################################################
summary(unlist(B[3,]))


###################################################
### code chunk number 21: fig:nb-sim-histo
###################################################
beta.hat <- do.call(rbind, B[4,])
par(mfrow=c(1,3), mar=c(5,4,1,2), las=1)
for (i in 0:2)
  hist(beta.hat[,i+1], breaks = 50, main = "",
       xlab = bquote(paste("Simulated ", beta[.(i)])))


###################################################
### code chunk number 22: Synzip
###################################################
# Synthetic Zero-inflated Poisson with logit component for 0's
# ============================================================
# syn.zip.r
library(MASS); library(pscl)
set.seed(1)
nobs <- 50000
x1 <- runif(nobs); x2 <- runif(nobs)
xb <- 2 + .75*x1 - 1.25*x2             # Poisson lin predictor
exb <- exp(xb)                         # Poisson fitted values
poy <- rpois(nobs, exb)                # Poisson
pdata <- data.frame(poy, x1, x2)       # Create data frame
pi <- 1/(1+exp(-(.9*x1 + .1*x2 + .2))) # Logit probabilities
pdata$bern <- runif(nobs) > pi         # Filter
zy <- pdata$bern * poy                 # Mixture
zip <- zeroinfl(zy ~ x1 + x2 | x1 + x2,
                dist = "poisson",
                data = pdata)
coef(summary(zip))


###################################################
### code chunk number 23: 7-bayes.rnw:640-641
###################################################
confint(zip)


###################################################
### code chunk number 24: medpar
###################################################
load("../../package/msme/data/medpar.rda")


###################################################
### code chunk number 25: 7-bayes.rnw:711-713 (eval = FALSE)
###################################################
## library(msme)
## data(medpar)


###################################################
### code chunk number 26: Setup
###################################################
library(arm)
# fit initial model to get coefficients
fit.1 <- glm(los ~ hmo + white, data = medpar, family=poisson)
fit.1$coef


###################################################
### code chunk number 27: arm.sim
###################################################
# simulation 1000 "random" values of each coefficient
n.sims <- 1000
sim.1 <- sim(fit.1, n.sims)
pcoef <- coef(sim.1)


###################################################
### code chunk number 28: hmo
###################################################
# get stats for hmo
hmo.coef <- pcoef[,2]
summary(hmo.coef)
quantile(hmo.coef, p = c(0.025, 0.975))


###################################################
### code chunk number 29: white
###################################################
# get stats for white
white.coef <- pcoef[,3]
summary(white.coef)
quantile(white.coef, p = c(0.025, 0.975))


###################################################
### code chunk number 30: fig:nb-sim-histo-arm
###################################################
par(mfrow=c(1,2), mar=c(5,4,1,2), las=1)
for (i in 1:2)
  hist(pcoef[,i+1], breaks = 30, main = "",
       xlab = c("Coefficient for HMO","Coefficient for White")[i])


###################################################
### code chunk number 31: 7-bayes.rnw:777-784
###################################################
library(COUNT)
library(nlme)

data(rwm5yr)
rwm5yr <- rwm5yr[,c(1,2,6:9)]
for(i in 1:6) class(rwm5yr[,i]) <- "numeric"
rwm5yr$id <- factor(rwm5yr$id)


###################################################
### code chunk number 32: 7-bayes.rnw:791-802
###################################################
test.lme <- lme(docvis ~ age + outwork + female + married,
                random = ~1 | id,
                method = "ML",
                data = rwm5yr)

summary(test.lme)$tTable

test.lm <- lm(docvis ~ age + outwork + female + married,
              data = rwm5yr)

anova(test.lme, test.lm) 


###################################################
### code chunk number 33: 7-bayes.rnw:811-827
###################################################
reps <- 1000

new.y <- simulate(test.lm, nsim = reps, seed = 100)
new.rwm <- rwm5yr

should.be.chi2 <- 
  sapply(new.y,
         function(x){
           out <- try({
           new.rwm$docvis <- x
           new.lm <- update(test.lm, data = new.rwm)
           new.lme <- update(test.lme, data = new.rwm)
           anova(new.lme, new.lm)$L.Ratio[2]})
           if (class(out) == "numeric") return(out)
           else return (NA)
         })


###################################################
### code chunk number 34: fig-not-chisq-1
###################################################
par(mfrow=c(1,2), mar=c(4,4,2,1), las=1)
qqplot(should.be.chi2, rchisq(10000, df = 1),
    xlab = "Simulated Dist.", ylab = "Theoretical Dist.")
abline(0, 1)
plot(ecdf(pchisq(should.be.chi2, df=1)),
     main = "Empirical Reference CDF")


###################################################
### code chunk number 35: not-chisq-1
###################################################
par(mfrow=c(1,2), mar=c(4,4,2,1), las=1)
qqplot(should.be.chi2, rchisq(10000, df = 1),
    xlab = "Simulated Dist.", ylab = "Theoretical Dist.")
abline(0, 1)
plot(ecdf(pchisq(should.be.chi2, df=1)),
     main = "Empirical Reference CDF")


###################################################
### code chunk number 36: 7-bayes.rnw:904-909
###################################################
nobs <- 100 
x1 <- runif(nobs)
x2 <- runif(nobs)
xb <- 2 + .75*x1 -1.25*x2
exb <- exp(xb)


###################################################
### code chunk number 37: 7-bayes.rnw:1021-1024
###################################################
m <- 0.5; y <- 4
m^y * exp(-m)/factorial(y)
dpois(y,m)


###################################################
### code chunk number 38: 7-bayes.rnw:1029-1031
###################################################
y*log(m) - m - lfactorial(y)
log(dpois(y, m))


###################################################
### code chunk number 39: 7-bayes.rnw:1035-1036
###################################################
dpois(y, m, log=TRUE)


###################################################
### code chunk number 40: medpar
###################################################
load("../../package/msme/data/medpar.rda")


###################################################
### code chunk number 41: 7-bayes.rnw:1059-1061 (eval = FALSE)
###################################################
## library(msme)
## data(medpar)


###################################################
### code chunk number 42: 7-bayes.rnw:1078-1082
###################################################
MLpoi <- glm(los ~ hmo, family = poisson, data = medpar)
summary(MLpoi)
logLik(MLpoi)
confint.default(MLpoi)


###################################################
### code chunk number 43: M-H
###################################################
jll.poisson <- function(y, mu, m, a) {
  dpois(x = y, lambda = mu, log = TRUE)
}

jll <- function(y, y.hat, ...) UseMethod("jll")

predict.expFamily <- function(y, b.hat, X, offset = 0) {
  lin.pred <- as.matrix(X) %*% b.hat + offset
  y.hat <- unlink(y, lin.pred)
  return(y.hat)
}

Sjll <- function(b.hat, X, y, offset = 0, ...) { 
  y.hat <- predict(y, b.hat, X, offset)
  sum(jll(y, y.hat, ...))
}

unlink <- function(y, eta, ...) UseMethod("unlink")
unlink.log <- function(y, eta, a=1, m=1) exp(eta)  


###################################################
### code chunk number 44: 7-bayes.rnw:1120-1130
###################################################
mh.formula <- los ~ hmo

family <- "poisson"
link <- "log"

mf <- model.frame(mh.formula, medpar)
y <- model.response(mf, "numeric")
X <- model.matrix(mh.formula, data = medpar)

class(y) <- c(family, link, "expFamily")


###################################################
### code chunk number 45: 7-bayes.rnw:1134-1135
###################################################
Sjll(c(1,0), X, y)


###################################################
### code chunk number 46: 7-bayes.rnw:1151-1159
###################################################
BLogPrior <- function(theta){
  alpha <- theta[1]
  beta <- theta[2]
  fprior.a <- dunif(alpha, -25, 30)
  fprior.b <- dunif(beta, -25, 30)
  fprior <- fprior.a * fprior.b # a, b independent
  if (fprior > 0) return(log(fprior))
}


###################################################
### code chunk number 47: 7-bayes.rnw:1177-1184
###################################################
  
nT <- 50000
Theta.t <- matrix(nrow = nT+1, ncol = 2)
Theta.star <- vector(length = 2)
current.Theta <- c(1, 0)
Theta.t[1,] <- current.Theta
acc <- 1


###################################################
### code chunk number 48: 7-bayes.rnw:1227-1247
###################################################
  
# Metropolis--Hastings MCMC algorithm
for (i in 1:nT) {
# Normal samples or draw for a proposal distribution
# We keep the proposal distributions distinct for clarity.
   Theta.star[1] <- rnorm(1, Theta.t[acc,1], 0.1)
   Theta.star[2] <- rnorm(1, Theta.t[acc,2], 0.1)
# Calculate log(R)
   logR <-
     (Sjll(Theta.star, X, y) + BLogPrior(Theta.star)) -
       (Sjll(Theta.t[acc,], X, y) + BLogPrior(Theta.t[acc,]))
# Draw new uniform random variate
   u <- runif(1)
# Compare u and r
   if (log(u) < logR) {
      acc <- acc + 1
      Theta.t[acc,] <- Theta.star
   }
}



###################################################
### code chunk number 49: 7-bayes.rnw:1252-1253
###################################################
(acc - 1) / nT


###################################################
### code chunk number 50: fig-met-hast-sim
###################################################
burn.in.length <- 500
burn <- 1:burn.in.length
use <- (burn.in.length+1):acc 
par(mfrow = c(2,2), mar = c(5,4,1,2), las = 1)
plot(Theta.t[burn, 1], xlab = "Index", ylab = "alpha", type= "l")
abline(h = coef(MLpoi)[1], lwd = 2)
plot(Theta.t[burn, 2], xlab = "Index", ylab = "beta", type = "l")
abline(h = coef(MLpoi)[2], lwd = 2)
plot(Theta.t[use, 1], xlab = "Index", ylab = "alpha", type = "l")
abline(h = coef(MLpoi)[1], lwd = 2) 
plot(Theta.t[use, 2], xlab = "Index", ylab = "beta", type = "l")
abline(h = coef(MLpoi)[2], lwd = 2)


###################################################
### code chunk number 51: met-hast-sim
###################################################
burn.in.length <- 500
burn <- 1:burn.in.length
use <- (burn.in.length+1):acc 
par(mfrow = c(2,2), mar = c(5,4,1,2), las = 1)
plot(Theta.t[burn, 1], xlab = "Index", ylab = "alpha", type= "l")
abline(h = coef(MLpoi)[1], lwd = 2)
plot(Theta.t[burn, 2], xlab = "Index", ylab = "beta", type = "l")
abline(h = coef(MLpoi)[2], lwd = 2)
plot(Theta.t[use, 1], xlab = "Index", ylab = "alpha", type = "l")
abline(h = coef(MLpoi)[1], lwd = 2) 
plot(Theta.t[use, 2], xlab = "Index", ylab = "beta", type = "l")
abline(h = coef(MLpoi)[2], lwd = 2)


###################################################
### code chunk number 52: fig-met-hast-out
###################################################
par(mfrow=c(1,2), mar = c(5,4,1,2), las = 1)
hist(Theta.t[use, 1], breaks=50, main = "", xlab = "Alpha")
hist(Theta.t[use, 2], breaks=50, main = "", xlab = "Beta")


###################################################
### code chunk number 53: met-hast-out
###################################################
par(mfrow=c(1,2), mar = c(5,4,1,2), las = 1)
hist(Theta.t[use, 1], breaks=50, main = "", xlab = "Alpha")
hist(Theta.t[use, 2], breaks=50, main = "", xlab = "Beta")


###################################################
### code chunk number 54: 7-bayes.rnw:1316-1320
###################################################
apply(Theta.t[use, ], 2, mean, na.rm = TRUE)
apply(Theta.t[use, ], 2, sd, na.rm = TRUE)
quantile (Theta.t[use, 1], na.rm = TRUE, probs = c(.025, .975))
quantile (Theta.t[use, 2], na.rm = TRUE, probs = c(.025, .975))


###################################################
### code chunk number 55: MCMCpack
###################################################
library(MCMCpack)
p.fit <- MCMCpoisson(los ~ hmo,  
                     burnin = 5000, mcmc = 100000, 
                     data = medpar)
summary(p.fit)



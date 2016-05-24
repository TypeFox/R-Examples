### R code from vignette source '2-estimation.rnw'

###################################################
### code chunk number 1: Setup
###################################################
options(repos="http://cran.r-project.org")

if(!require(Hmisc, quietly=TRUE)) install.packages("Hmisc")
if(!require(xtable, quietly=TRUE)) install.packages("xtable")

rm( list=ls() )

options(width = 67)


###################################################
### code chunk number 2: 2-estimation.rnw:291-294
###################################################
y <- 5
log(factorial(y))
lfactorial(y)


###################################################
### code chunk number 3: 2-estimation.rnw:299-300
###################################################
lfactorial


###################################################
### code chunk number 4: 2-estimation.rnw:402-405
###################################################
jll.watson <- function(theta, x) {
  sum(log(1 + theta) - log(theta) - 2*log(1 + x / theta))
}


###################################################
### code chunk number 5: 2-estimation.rnw:412-420
###################################################
watson.fit <- function(x, ...) {
  optim(0.5,
        jll.watson,
        x = x,
        method = "Brent",
        lower = 0, upper = 1,
        control = list(fnscale= -1), ...)
}


###################################################
### code chunk number 6: 2-estimation.rnw:428-429
###################################################
large.sample <- rep(1:10, 10)/20


###################################################
### code chunk number 7: 2-estimation.rnw:432-433
###################################################
large.sample.fit <- watson.fit(large.sample)


###################################################
### code chunk number 8: 2-estimation.rnw:444-445
###################################################
large.sample.fit


###################################################
### code chunk number 9: 2-estimation.rnw:626-630
###################################################
large.sample <- rep(1:10, 10)/20
large.sample.fit <- watson.fit(large.sample, hessian = TRUE)
large.sample.fit$par
(large.se <- sqrt(diag(solve(-large.sample.fit$hessian))))


###################################################
### code chunk number 10: 2-estimation.rnw:634-635
###################################################
large.sample.fit$par + c(-1,1) * large.se


###################################################
### code chunk number 11: 2-estimation.rnw:642-646
###################################################
small.sample <- rep(1:10, 1)/20
small.sample.fit <- watson.fit(small.sample, hessian = TRUE)
small.sample.fit$par
(small.se <- sqrt(diag(solve(-small.sample.fit$hessian))))


###################################################
### code chunk number 12: 2-estimation.rnw:650-651
###################################################
small.sample.fit$par + c(-1,1) * small.se


###################################################
### code chunk number 13: 2-estimation.rnw:699-701
###################################################
mll.watson <- function (x, data)
  sapply(x, function(y) jll.watson(y, data))


###################################################
### code chunk number 14: 2-estimation.rnw:705-708
###################################################
profiles <- data.frame(thetas = (1:1000)/1000)
profiles$ll.large <- mll.watson(profiles$thetas, large.sample)
profiles$ll.small <- mll.watson(profiles$thetas, small.sample)


###################################################
### code chunk number 15: 2-estimation.rnw:718-719
###################################################
for (i in 2:3) profiles[,i] <- profiles[,i] - max(profiles[,i])


###################################################
### code chunk number 16: fig-profile-compare
###################################################
par(las = 1, mar = c(4,4,2,1))
plot(ll.small ~ thetas, data = profiles,
     type = "l", ylim = c(-4,0),
     ylab = "Log-Likelihood (Scaled)",
     xlab = expression(paste("Candidate Values of ", theta)))
lines(ll.large ~ thetas, data = profiles, lty = 2)
abline(h = -1.92)
abline(v = large.sample.fit$par)


###################################################
### code chunk number 17: profile-compare
###################################################
par(las = 1, mar = c(4,4,2,1))
plot(ll.small ~ thetas, data = profiles,
     type = "l", ylim = c(-4,0),
     ylab = "Log-Likelihood (Scaled)",
     xlab = expression(paste("Candidate Values of ", theta)))
lines(ll.large ~ thetas, data = profiles, lty = 2)
abline(h = -1.92)
abline(v = large.sample.fit$par)


###################################################
### code chunk number 18: 2-estimation.rnw:807-827
###################################################
set.seed(1234)
gamma.sample <- rgamma(1000, scale = 4, shape = 2)

jll.gamma <- function(params, data) {
  sum(dgamma(data,
             scale = params[2],
             shape = params[1],
             log = TRUE))
}

gamma.fit <- function(data, ...) {
  optim(c(2,2),
        jll.gamma, 
        data = data,
        control = list(fnscale = -1), ...)
}

test <- gamma.fit(gamma.sample, hessian = TRUE)

test


###################################################
### code chunk number 19: 2-estimation.rnw:838-851
###################################################
jll.gamma.shape <- function(params, alpha, data) {
  sum(dgamma(data, scale = params, shape = alpha, log = TRUE))
}

gamma.fit.shape <- function(alpha, data, ...) {
  optim(1,
        jll.gamma.shape, 
        data = data,
        alpha = alpha,
        method = "Brent",
        lower = 0, upper = 10,
        control = list(fnscale = -1), ...)$value
}


###################################################
### code chunk number 20: 2-estimation.rnw:857-864
###################################################
gamma.profile <-
  data.frame(candidates = seq(1.9, 2.3, length.out = 100))

gamma.profile$profile.right <-
  sapply(gamma.profile$candidates,
         gamma.fit.shape,
         data = gamma.sample)


###################################################
### code chunk number 21: 2-estimation.rnw:870-882
###################################################
gamma.fit.shape.wrong <- function(alpha, data, mle) {
  sum(dgamma(data, scale = mle, shape = alpha, log = TRUE))
}
  
gamma.profile$profile.wrong <-
  sapply(gamma.profile$candidates,
         gamma.fit.shape.wrong,
         data = gamma.sample,
         mle = test$par[2])

for (i in 2:3)
  gamma.profile[,i] <- gamma.profile[,i] - max(gamma.profile[,i])


###################################################
### code chunk number 22: fig-gamma-profile
###################################################
par(las = 1, mar = c(4,4,2,1))
plot(profile.right ~ candidates, type="l",
     data = gamma.profile, ylim = c(-2, 0), 
     ylab = "Log-Likelihood (Scaled)",
     xlab = "Candidate Values of a")
lines(profile.wrong ~ candidates, type="l",
      data = gamma.profile,
      lty = 2)
abline(h = -1.92)
abline(v = test$par[1])


###################################################
### code chunk number 23: gamma-profile
###################################################
par(las = 1, mar = c(4,4,2,1))
plot(profile.right ~ candidates, type="l",
     data = gamma.profile, ylim = c(-2, 0), 
     ylab = "Log-Likelihood (Scaled)",
     xlab = "Candidate Values of a")
lines(profile.wrong ~ candidates, type="l",
      data = gamma.profile,
      lty = 2)
abline(h = -1.92)
abline(v = test$par[1])


###################################################
### code chunk number 24: 2-estimation.rnw:967-975
###################################################
testnum <- 2
set.seed(2468)
r1 <- runif(testnum)
r2 <- runif(testnum)
c(r1, r2)
testnum <- 4
set.seed(2468)
(r12 <- runif(testnum))


###################################################
### code chunk number 25: 2-estimation.rnw:1008-1012
###################################################
rndexp <- function(obs = 10000, shape = 3) {
  xe <- -(shape)*log(runif(obs))
  return(xe)
}


###################################################
### code chunk number 26: 2-estimation.rnw:1015-1017
###################################################
set.seed(1)
rndexp(10, 3)


###################################################
### code chunk number 27: 2-estimation.rnw:1028-1032
###################################################
rndchi2 <- function(obs = 10000, dof = 3) {
  z2 <- matrix(qnorm(runif(obs*dof))^2, nrow = obs)
  return(rowSums(z2))
}


###################################################
### code chunk number 28: 2-estimation.rnw:1035-1037
###################################################
set.seed(1)
rndchi2(10, 2)


###################################################
### code chunk number 29: 2-estimation.rnw:1044-1054
###################################################
rndpoi  <- function(mu) {
  g <- exp(-mu) 
  em <- -1 
  t <- 1   
  while(t > g) {
    em <- em + 1
    t <- t * runif(1)
  }
  return(floor(em + 0.5))
}


###################################################
### code chunk number 30: 2-estimation.rnw:1057-1059
###################################################
set.seed(1)
rndpoi(4)


###################################################
### code chunk number 31: 2-estimation.rnw:1067-1077
###################################################
rndpoi  <- function(obs = 50000, mu = 4) {
  g <- exp(-mu) 
  em <- rep(-1, obs) 
  t <- rep(1, obs)   
  while(any(t > g)) {
    em <- em + (t > g)
    t[t > g] <- t[t > g] * runif(sum(t > g))
  }
  return(floor(em + 0.5))
}


###################################################
### code chunk number 32: 2-estimation.rnw:1081-1086
###################################################
set.seed(1)
xp <- rndpoi(10000, 4)
str(xp)
mean(xp)
var(xp)


###################################################
### code chunk number 33: 2-estimation.rnw:1162-1279 (eval = FALSE)
###################################################
## 
## jll.watson <- function(theta, x) {
##   sum(log(1 + theta) - log(theta) - 2*log(1 + x / theta))
## }
## 
## dwatson <- function(x, theta) {
##   (1 + theta) / theta / (1 + x / theta)^2
## }
## 
## pwatson <- function(q, theta) {
##   integrate(function(x) dwatson(x, theta), 
##             lower = 0, upper = q)$value
##   }
## 
## qwatson <- function(p, theta) {
##   uniroot(function(x) pwatson(x, theta) - p, 
##           lower = .Machine$double.eps, upper = 1)$root
## }
## 
## qwatson <- function(p, theta) {
##   (1 / ( 1 - p / (1 + theta) ) - 1) * theta
## }
## 
## rwatson <- function(n, theta) {
##   sapply(runif(n), function(x) qwatson(x, theta))
## }
## 
## # library(boot)
## 
## watson.fit <- function(x, ...) {
##   optim(0.5,
##         jll.watson,
##         x = x,
##         method = "Brent",
##         lower = 0, upper = 1,
##         control = list(fnscale= -1), ...)
## }
## 
## small.sample <- rep(1:10, 1)/20
## 
## small.sample.fit <- watson.fit(small.sample, hessian = TRUE)
## 
## reps <- 999
## 
## ## PB doesn't include variability of MLE??
## 
## watson.gen <- function(data, mle) {
##   rwatson(length(data), mle)
## }
## 
## watson.boot.small <-
##  sapply(1:reps,
##         function (x) {
##           new.sample <- watson.gen(small.sample,
##                                    small.sample.fit$par)
##           fit <- watson.fit(new.sample, hessian=TRUE)
##           pivot <- (fit$par - small.sample.fit$par) /
##                          sqrt(diag(solve(-fit$hessian)))
##           return(pivot)
##         })
## 
## #watson.boot.small <-    ### Convert this to a loop
## #   boot(small.sample,
## #        function(x) {
## #          fit <- watson.fit(x, hessian=TRUE)
## #          pivot <- (fit$par - small.sample.fit$par) /
## #            sqrt(diag(solve(-fit$hessian)))
## #        },
## #        R = 999,
## #        sim = "parametric",
## #        ran.gen = watson.gen,
## #        mle = small.sample.fit$par)
## 
## # watson.boot.small
## 
## small.sample.fit$par - sqrt(diag(solve(-small.sample.fit$hessian))) *
##    quantile(watson.boot.small, p = c(0.975, 0.025))
## 
## par(mfrow = c(1,3), mar=c(4,4,2,1), las=1)
## plot(density(watson.boot.small))
## plot(ecdf(watson.boot.small))
## qqnorm(watson.boot.small)
## qqline(watson.boot.small)
## 
## # plot(watson.boot.small)
## 
## ## Non-parametric
## 
## watson.gen <- function(data, mle) {
##   sample(data, replace=TRUE)
## }
## 
## watson.boot.small <-
##  sapply(1:reps,
##         function (x) {
##           new.sample <- watson.gen(small.sample,
##                                    small.sample.fit$par)
##           fit <- watson.fit(new.sample, hessian=TRUE)
##           pivot <- (fit$par - small.sample.fit$par) /
##                          sqrt(diag(solve(-fit$hessian)))
##           return(pivot)
##         })
## 
## small.sample.fit$par - sqrt(diag(solve(-small.sample.fit$hessian))) *
##    quantile(watson.boot.small, p = c(0.975, 0.025))
## 
## par(mfrow = c(1,3), mar=c(4,4,2,1), las=1)
## plot(density(watson.boot.small))
## plot(ecdf(watson.boot.small))
## qqnorm(watson.boot.small)
## qqline(watson.boot.small)
## 
## 
## 
## 
## 
## 



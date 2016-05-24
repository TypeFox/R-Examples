### R code from vignette source '4-logistic.rnw'

###################################################
### code chunk number 1: Setup
###################################################
options(repos="http://cran.r-project.org")
options(useFancyQuotes = FALSE)

if(!require(Hmisc, quietly=TRUE)) install.packages("Hmisc")
if(!require(xtable, quietly=TRUE)) install.packages("xtable")

rm(list=ls())

options(width = 62)


###################################################
### code chunk number 2: 4-logistic.rnw:699-700
###################################################
load("../../package/msme/data/titanic.rda")


###################################################
### code chunk number 3: 4-logistic.rnw:703-705 (eval = FALSE)
###################################################
## library(msme)
## data(titanic)


###################################################
### code chunk number 4: 4-logistic.rnw:708-721
###################################################
y <- titanic$survived
x <- titanic$age

mu <- rep(mean(y), nrow(titanic))    # initialize mu
eta <- log(mu/(1-mu))                # initialize eta
for (i in 1:4) {                     # loop for 4 iterations
    w <-  mu*(1-mu)                  # weight = variance
    z <- eta + (y - mu)/(mu*(1-mu))  # working response
    mod <- lm(z ~ x, weights = w)    # weighted regression
    eta <- mod$fit                   # linear predictor
    mu <- 1/(1+exp(-eta))            # fitted value
    cat(i, coef(mod), "\n")          # displayed iteration log
}


###################################################
### code chunk number 5: 4-logistic.rnw:730-731
###################################################
coef(summary(mod))


###################################################
### code chunk number 6: 4-logistic.rnw:735-736
###################################################
confint(mod)


###################################################
### code chunk number 7: 4-logistic.rnw:743-747
###################################################
glm.test <- glm(survived ~ age,
                family = binomial,
                data = titanic)
coef(summary(glm.test))


###################################################
### code chunk number 8: 4-logistic.rnw:757-758
###################################################
getAnywhere(vcov.lm)


###################################################
### code chunk number 9: 4-logistic.rnw:764-765
###################################################
(se <- sqrt(diag(summary(mod)$cov.unscaled)))


###################################################
### code chunk number 10: 4-logistic.rnw:801-806
###################################################
glm.qb <- glm(survived ~ age,
              family = quasibinomial,
              data = titanic)
summary(glm.qb)$dispersion
coef(summary(glm.qb))


###################################################
### code chunk number 11: 4-logistic.rnw:832-833
###################################################
confint(mod)


###################################################
### code chunk number 12: 4-logistic.rnw:843-844
###################################################
confint(glm.test)


###################################################
### code chunk number 13: 4-logistic.rnw:858-859
###################################################
confint.default(glm.test)


###################################################
### code chunk number 14: 4-logistic.rnw:886-888
###################################################
mod$coef - 1.96*se
mod$coef + 1.96*se


###################################################
### code chunk number 15: 4-logistic.rnw:905-947
###################################################
irls_logit <- function(formula, data, tol = 0.000001) {# arguments

## Set up the model components
    mf <- model.frame(formula, data)        # define model frame as mf
    y <- model.response(mf, "numeric")      # set model response as y
    X <- model.matrix(formula, data = data) # predictors in matrix X

## Check for missing values; stop if any.
    if (any(is.na(cbind(y, X)))) stop("Some data are missing.")

## Initialize mu, eta, the deviance, etc.
    mu <- rep(mean(y), length(y))
    eta <- log(mu/(1-mu))                             
    dev <- 2 * sum(y*log(1/mu) +
           (1 - y) * log(1/(1-mu)))
    deltad <- 1                                       
    i <- 1                                            

## Loop through the IRLS algorithm
    while (abs(deltad) > tol ) {                  # IRLS loop begin
        w <-  mu * (1-mu)                         # weight
        z <- eta + (y - mu)/w                     # working response
        mod <- lm(z ~ X-1, weights = w)           # weighted regression
        eta <- mod$fit                            # linear predictor
        mu <- 1/(1+exp(-eta))                     # fitted value
        dev.old <- dev                            
        dev <- 2 * sum(y * log(1/mu) +
                       (1 - y) * log(1/(1 - mu))) # deviance
        deltad <- dev - dev.old                   # change
        cat(i, coef(mod), deltad, "\n")           # iteration log
        i <- i + 1                                # iterate
    }

## Build some post-estimation statistics
    df.residual <- summary(mod)$df[2]
    pearson.chi2 <- sum((y - mu)^2 / (mu * (1 - mu))) / df.residual

## Return a compact result
    result <- list(coefficients = coef(mod),
                   se.beta.hat = sqrt(diag(summary(mod)$cov.unscaled)))
  return(result)
}


###################################################
### code chunk number 16: 4-logistic.rnw:972-973
###################################################
load("../../package/msme/data/medpar.rda")


###################################################
### code chunk number 17: 4-logistic.rnw:978-980 (eval = FALSE)
###################################################
## library(msme)
## data(medpar)


###################################################
### code chunk number 18: 4-logistic.rnw:985-987
###################################################
i.logit <- irls_logit(died ~ hmo + white,
                      data = medpar)


###################################################
### code chunk number 19: 4-logistic.rnw:998-999
###################################################
i.logit


###################################################
### code chunk number 20: 4-logistic.rnw:1004-1006
###################################################
with(i.logit, coefficients - 1.96 * se.beta.hat)
with(i.logit, coefficients + 1.96 * se.beta.hat)


###################################################
### code chunk number 21: 4-logistic.rnw:1010-1012
###################################################
(Z <- with(i.logit, coefficients / se.beta.hat))
(pvalues <- 2*pnorm(abs(Z), lower.tail = FALSE))


###################################################
### code chunk number 22: 4-logistic.rnw:1016-1020
###################################################
glm.logit <- glm(died ~ hmo + white,
                   family = binomial,
                   data = medpar)
coef(summary(glm.logit))


###################################################
### code chunk number 23: 4-logistic.rnw:1030-1031
###################################################
confint.default(glm.logit)


###################################################
### code chunk number 24: 4-logistic.rnw:1038-1046
###################################################
estim <-coef(summary(glm.logit))                                           
beta <- estim[,1]                                                          
se <- estim[,2]                                                            
confint.model <- (cbind(beta, 
                        beta - se*qnorm(.975), 
                        beta + se*qnorm(.975)))                                                             
colnames(confint.model) <- c("Beta", "low_CI", "Hi_CI")                    
confint.model                                                              


###################################################
### code chunk number 25: glogit
###################################################
irls_glogit <- function(formula, data,
                        tol = 0.000001, offset = 0, m = 1 ) {

## Set up the model components
   mf <- model.frame(formula, data)
   y <- model.response(mf, "numeric")
   X <- model.matrix(formula, data = data)

## Check for missing values; stop if any.
   if (any(is.na(cbind(y, X)))) stop("Some data are missing.")

## Check for nonsensical values; stop if any.
   if (any(y < 0 | y > m)) stop("Some data are absurd.")

## Initialize mu, eta, the deviance, etc.
   mu <- rep(mean(y), length(y))
   eta <- log(mu/(m - mu))
   dev <- 2 * sum(y * log(pmax(y,1)/mu) +
              (m - y)* log(pmax(y,1)/(m - mu)) )
   deltad <- 1
   i <- 1

## Loop through the IRLS algorithm
   while (abs(deltad) > tol ) {
     w <- mu*(m-mu)/m
     z <- eta + (y-mu)/w - offset
     mod <- lm(z ~ X-1, weights=w)
     eta <- mod$fit + offset
     mu <- m/(1+exp(-eta))
     dev.old <- dev
     dev <- 2 * sum(y*log(pmax(y,1)/mu) +
                    (m - y)* log(pmax(y,1)/(m - mu)))
     deltad <- dev - dev.old
     i <- i + 1
   }

## Build some post-estimation statistics
   df.residual <- summary(mod)$df[2]
   pearson.chi2 <- sum((y - mu)^2 /
                       (mu * (1 - mu/m))) / df.residual

## Return a brief result
   return(list(coef = coef(mod),
               se = sqrt(diag(summary(mod)$cov.unscaled)),
               pearson.chi2 = pearson.chi2))
 }


###################################################
### code chunk number 26: 4-logistic.rnw:1240-1241
###################################################
pmax(c(0, 1, 1), c(1, -1, 0))


###################################################
### code chunk number 27: 4-logistic.rnw:1273-1274
###################################################
load("../../package/msme/data/doll.rda")


###################################################
### code chunk number 28: 4-logistic.rnw:1277-1279 (eval = FALSE)
###################################################
## library(msme)
## data(doll)


###################################################
### code chunk number 29: irls - doll
###################################################
i.glog <- irls_glogit(deaths ~ smokes + ordered(age),
                      data = doll,
                      m = doll$pyears)
i.glog


###################################################
### code chunk number 30: 4-logistic.rnw:1307-1308
###################################################
with(i.glog, coef / se)


###################################################
### code chunk number 31: 4-logistic.rnw:1316-1321
###################################################
glm.glog <- glm(cbind(deaths, pyears - deaths) ~ 
                smokes + ordered(age),
                data = doll,
                family = binomial)
coef(summary(glm.glog))


###################################################
### code chunk number 32: 4-logistic.rnw:1379-1384
###################################################
jll <- function(y, mu, m, a) UseMethod("jll")
linkFn <- function(mu, m, a) UseMethod("linkFn")
lPrime <- function(mu, m, a) UseMethod("lPrime")
unlink <- function(y, eta, m, a) UseMethod("unlink")
variance <- function(mu, m, a) UseMethod("variance")


###################################################
### code chunk number 33: 4-logistic.rnw:1404-1407
###################################################
devianceResids <- function(y, mu, m, a)
  sign(y - mu) * sqrt(2 * abs(jll(y, mu, m, a) -
                              jll(y, y,  m, a)))


###################################################
### code chunk number 34: 4-logistic.rnw:1410-1412
###################################################
devIRLS <- function(object, ...)
   sum(devianceResids(object, ...)^2)


###################################################
### code chunk number 35: 4-logistic.rnw:1435-1440
###################################################
initialize <- function(y, m) {
   ret.y <- rep(mean(y), length(y))
   class(ret.y) <- class(y)
   ret.y
}


###################################################
### code chunk number 36: 4-logistic.rnw:1447-1448
###################################################
variance.binomial <- function(mu, m, a) mu * (1 - mu/m)


###################################################
### code chunk number 37: 4-logistic.rnw:1456-1458
###################################################
jll.binomial <- function(y, mu, m, a) 
  dbinom(x = y, size = m, prob = mu / m, log = TRUE)


###################################################
### code chunk number 38: logit bits
###################################################
linkFn.logit <- function(mu, m, a) log(mu / (m - mu))
lPrime.logit <- function(mu, m, a) m / (mu * (m - mu))
unlink.logit <- function(y, eta, m, a) m / (1 + exp(-eta))


###################################################
### code chunk number 39: irls-def
###################################################
irls <- function(formula, data, family, link,
                 tol = 1e-6,
                 offset = 0,
                 m = 1,
                 a = 1,
                 verbose = 0) {

### Prepare the model components as previously
  mf <- model.frame(formula, data)
  y <- model.response(mf, "numeric")
  X <- model.matrix(formula, data = data)

### Check for missing values 
  if (any(is.na(cbind(y, X)))) stop("Some data are missing.")

### Arrange the class information
  class(y) <- c(family, link, "expFamily")

### Establish a start point
  mu <- initialize(y, m)
  eta <- linkFn(mu, m, a)
  dev <- devIRLS(y, mu, m, a)
  deltad <- 1
  i <- 1

### IRLS loop (as before)
  while (abs(deltad) > tol ) {
    w <-  1 / (variance(mu, m, a) * lPrime(mu, m, a)^2)
    z <- eta + (y - mu) * lPrime(mu, m, a) - offset
    mod <- lm(z ~ X - 1, weights = w)
    eta <- mod$fit + offset
    mu <- unlink(y, eta, m, a)
    dev.old <- dev
    dev <- devIRLS(y, mu, m, a)
    deltad <- dev - dev.old
    if(verbose > 0) cat(i, coef(mod), deltad, "\n")
    i <- i + 1
  }

### Post-estimation statistics
  df.residual <- summary(mod)$df[2]
  pearson.chi2 <- sum((y - mu)^2 /
                  variance(mu, m, a)) / df.residual
  ll <- sum(jll(y, mu, m, a))

### Return a rich object --- allows use of print.glm
  result <-
    list(coefficients = coef(mod),
         se.beta.hat = sqrt(diag(summary(mod)$cov.unscaled)),
         model = mod,
         call = match.call(),
         nobs = length(y),
         eta = eta,
         mu = mu,
         df.residual = df.residual,
         df.null = length(y) - 1,
         deviance = dev,
         null.deviance = NA,
         p.dispersion = pearson.chi2,
         pearson = pearson.chi2 * df.residual,
         loglik = ll,
         family = list(family = family),
         X = X,
         i = i - 1,
         residuals = devianceResids(y, mu, m, a),
         aic = -2 * ll + 2 * summary(mod)$df[1])
  class(result) <- c("msme","glm")
  return(result)
}


###################################################
### code chunk number 40: binomial
###################################################
irls.logit <- irls(died ~ hmo + white,
                  family = "binomial", link = "logit",
                  data = medpar)
irls.logit


###################################################
### code chunk number 41: 4-logistic.rnw:1646-1649
###################################################
with(irls(died ~ 1,
          family = "binomial", link = "logit",
          data = medpar), c(deviance, df.residual))


###################################################
### code chunk number 42: 4-logistic.rnw:1665-1708
###################################################
summary.msme <- function(object, ...) {

### Create a coefficient table  
  z <- with(object, coefficients / se.beta.hat)
  zTable <-
      with(object, 
           data.frame(Estimate = coefficients,
                      SE = se.beta.hat,
                      Z = z,
                      p = 2 * pnorm(-abs(z)),
                      LCL = coefficients - 1.96 * se.beta.hat,
                      UCL = coefficients + 1.96 * se.beta.hat))
  rownames(zTable) <- colnames(object$X)

### Prepare part of the coefficient table for printing  
  z.print <- zTable
  z.print$p <- formatC(z.print$p, digits = 3, format="g")

### Build a list of output objects
  summ <- list(call = object$call,
               coefficients = zTable,
               deviance = object$deviance,
               null.deviance = object$null.deviance,
               df.residual = object$df.residual,
               df.null = object$df.null) 

### Write out a set of results
  cat("\nCall:\n")
  print(object$call)
  cat("\nDeviance Residuals:\n")
  print(summary(as.numeric(object$residuals)))
  cat("\nCoefficients:\n")
  print(z.print, digits = 3, ...)
  cat("\nNull deviance:", summ$null.deviance,
      " on ", summ$df.null, "d.f.") 
  cat("\nResidual deviance:", summ$deviance,
      " on ", summ$df.residual, "d.f.") 
  cat("\nAIC: ", object$aic)
  cat("\n\nNumber of optimizer iterations: ", object$i, "\n\n")

### Return the list but do not print it.  
  return(invisible(summ))
}


###################################################
### code chunk number 43: 4-logistic.rnw:1806-1809
###################################################
linkFn.probit <- function(mu, m, a) qnorm(mu / m)
lPrime.probit <- function(mu, m, a) 1 / (m * dnorm(qnorm(mu/m)))
unlink.probit <- function(y, eta, m, a) m * pnorm(eta)


###################################################
### code chunk number 44: 4-logistic.rnw:1812-1815
###################################################
i.probit <- irls(died ~ hmo + white, 
                  family = "binomial", link = "probit",
                  data = medpar)


###################################################
### code chunk number 45: 4-logistic.rnw:1847-1850
###################################################
coef(summary(i.glm <- glm(died ~ hmo + white, 
                          family = binomial(probit),
                          data = medpar)))


###################################################
### code chunk number 46: 4-logistic.rnw:1945-1948
###################################################
linkFn.log <- function(mu, m, a) log(mu)
lPrime.log <- function(mu, m, a) 1/mu
unlink.log <- function(y, eta, m, a) exp(eta)  


###################################################
### code chunk number 47: 4-logistic.rnw:1951-1956
###################################################
variance.negBinomial <- function(mu, m, a) mu + a*mu^2

jll.negBinomial <- function(y, mu, m, a) {
  dnbinom(y, mu = mu, size = 1 / a, log = TRUE)
}


###################################################
### code chunk number 48: 4-logistic.rnw:1965-1968
###################################################
library(MASS)
ml.nb <- glm.nb(los ~ hmo + white, 
                data = medpar)


###################################################
### code chunk number 49: 4-logistic.rnw:2010-2011
###################################################
1 / ml.nb$theta


###################################################
### code chunk number 50: 4-logistic.rnw:2016-2020
###################################################
irls.nb <- irls(los ~ hmo + white, a = 0.4846026,
                family = "negBinomial", link = "log",
                data = medpar)
nb.summ <- summary(irls.nb)


###################################################
### code chunk number 51: 4-logistic.rnw:2052-2053
###################################################
ml.nb$theta


###################################################
### code chunk number 52: 4-logistic.rnw:2057-2060
###################################################
glm.nb <- glm(los ~ hmo + white,
              data = medpar,
              family = negative.binomial(2.063547))


###################################################
### code chunk number 53: 4-logistic.rnw:2098-2099
###################################################
sqrt(summary(glm.nb)$dispersion)


###################################################
### code chunk number 54: 4-logistic.rnw:2105-2106
###################################################
nb.summ$coefficients[,"SE"] * sqrt(irls.nb$p.dispersion)


###################################################
### code chunk number 55: 4-logistic.rnw:2210-2212 (eval = FALSE)
###################################################
## library(msme)
## data(heart)


###################################################
### code chunk number 56: 4-logistic.rnw:2214-2215
###################################################
load("../../package/msme/data/heart.rda")


###################################################
### code chunk number 57: heart
###################################################
heart


###################################################
### code chunk number 58: 4-logistic.rnw:2228-2233
###################################################
heart.nb <- irls(death ~ anterior + hcabg + factor(killip),
                 a = 0.0001,
                 offset = log(heart$cases),
                 family = "negBinomial", link = "log",
                 data = heart)


###################################################
### code chunk number 59: 4-logistic.rnw:2268-2273
###################################################
h.glm <- glm(death ~ anterior + hcabg + factor(killip),
             family = negative.binomial(10000),             
             offset = log(cases),
             data = heart)
coef(summary(h.glm))


###################################################
### code chunk number 60: 4-logistic.rnw:2279-2280
###################################################
sqrt(heart.nb$p.dispersion)


###################################################
### code chunk number 61: 4-logistic.rnw:2283-2284
###################################################
with(heart.nb, se.beta.hat * sqrt(p.dispersion))


###################################################
### code chunk number 62: 4-logistic.rnw:2288-2321 (eval = FALSE)
###################################################
## irls_nb <- function(formula, data, offset=0, a=1) {
##   mf <- model.frame(formula, data)
##   y <- model.response(mf, "numeric")
##   X <- model.matrix(formula, data = data)
##   if (any(is.na(cbind(y, X)))) stop("Some data are missing.")
##   mu <- (y + mean(y))/2
##   eta <- log(mu)
##   dev <- 2 * sum(y*log(y/mu) - (y + 1/a) * log((1+a*y)/(1+a*mu)))
##   deltad <- 1
##   i <- 1
##   while (abs(deltad) > 0.000001 ) {
##     w <- 1/((mu + a*mu*mu) * (1/mu)^2)
##     z <- eta + (y - mu)*(1/mu) - offset
##     mod <- lm(z ~ X-1, weights=w)
##     eta <- mod$fit + offset
##     mu <- exp(eta)
##     dev.old <- dev
##     dev <- 2 * sum(y*log(y/mu) - (y + 1/a) * log((1+a*y)/(1+a*mu)))
##     deltad <- dev - dev.old
##     cat(i, coef(mod), deltad, "\n")
##     i <- i + 1
##   }
##   df.residual <- summary(mod)$df[2]
##   pearson.chi2 <- sum((y - mu)^2 / (mu + a*mu^2)) / df.residual
##   return(list(coef = coef(mod),
##               se = sqrt(diag(summary(mod)$cov.unscaled)) 
##               ))
## }  
## heart.test <- irls_nb(death ~ anterior + hcabg + factor(killip),
##                 a = 0.0001,
##                 offset = log(heart$cases),
##                 data = heart)
## heart.test


###################################################
### code chunk number 63: 4-logistic.rnw:2835-2836
###################################################
irls.nb$p.dispersion


###################################################
### code chunk number 64: 4-logistic.rnw:2871-2872
###################################################
anova(glm.glog)


###################################################
### code chunk number 65: 4-logistic.rnw:2876-2877
###################################################
12.10 / 4


###################################################
### code chunk number 66: 4-logistic.rnw:2882-2884
###################################################
sum(residuals(glm.glog, type="pearson")^2) / 
    glm.glog$df.residual


###################################################
### code chunk number 67: 4-logistic.rnw:2926-2931
###################################################
P__disp <- function(x) {
   pr <- sum(residuals(x, type="pearson")^2)
   dispersion <- pr / x$df.residual
   return(c(pearson.chi2 = pr, dispersion = dispersion))
}


###################################################
### code chunk number 68: 4-logistic.rnw:2938-2939
###################################################
P__disp(glm.glog)


###################################################
### code chunk number 69: 4-logistic.rnw:3118-3122
###################################################
jll.poisson <- function(y, mu, m, a) {
  dpois(x = y, lambda = mu, log = TRUE)
}
variance.poisson <- function(mu, m, a) mu


###################################################
### code chunk number 70: heart.p
###################################################
heart.p <- irls(death ~ anterior + hcabg + factor(killip),
                offset = log(heart$cases),
                family = "poisson", link = "log",
                data = heart)


###################################################
### code chunk number 71: 4-logistic.rnw:3137-3139
###################################################
LR <- -2*(heart.p$loglik - heart.nb$loglik)
pchisq(LR, 1, lower.tail = FALSE)/2


###################################################
### code chunk number 72: 4-logistic.rnw:3142-3172 (eval = FALSE)
###################################################
## logLik.msme <- function(x, ...) {
##   val <- x$loglik
##   attr(val, "nall") <- nrow(x$X)
##   attr(val, "nobs") <- nrow(x$X)
##   attr(val, "df") <- length(x$coefficients)
##   class(val) <- "logLik"  
##   val
## }
## 
## alrt <- function(x1, x2, ...) {
##   jll1 <- logLik(x1)
##   jll2 <- logLik(x2)
##   df1 <- attr(jll1, "df")
##   df2 <- attr(jll2, "df")
##   jll.diff <- abs(c(jll1) - c(jll2))
##   df.diff <- abs(df1 - df2)
##   p.value <- 1 - pchisq(2 * jll.diff, df = df.diff) 
##   results <- list(out.tab = data.frame(model = c(1,2),
##                                        jll = c(jll1, jll2),
##                                        df = c(df1, df2)),
##                   jll.diff = jll.diff,
##                   df.diff = df.diff,
##                   p = p.value)
##   cat("\nLL of model 1: ", jll1, " df: ", df1, 
##       "\nLL of model 2: ", jll2, " df: ", df2, 
##       "\nDifference: ", jll.diff, " df: ", df.diff, 
##       "\np-value against H_0: no difference between models ", 
##       p.value, "\n")
##   return(invisible(results))
## }



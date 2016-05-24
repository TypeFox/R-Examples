### R code from vignette source '5-advanced-glm.rnw'

###################################################
### code chunk number 1: 5-advanced-glm.rnw:7-9
###################################################
rm(list = ls())
options(useFancyQuotes = FALSE)


###################################################
### code chunk number 2: Setup
###################################################
options(repos="http://cran.r-project.org")

if(!require(Hmisc, quietly=TRUE)) install.packages("Hmisc")
if(!require(xtable, quietly=TRUE)) install.packages("xtable")

options(width = 65)


###################################################
### code chunk number 3: 5-advanced-glm.rnw:104-107
###################################################
jll.bernoulli <- function(y, y.hat, ...) {
  dbinom(x = y, size = 1, prob = y.hat, log = TRUE)
}


###################################################
### code chunk number 4: 5-advanced-glm.rnw:112-113
###################################################
jll <- function(y, y.hat, ...) UseMethod("jll")


###################################################
### code chunk number 5: 5-advanced-glm.rnw:129-133
###################################################
Sjll <- function(b.hat, X, y, offset = 0, ...) {
  y.hat <- predict(y, b.hat, X, offset)
  sum(jll(y, y.hat, ...))
}


###################################################
### code chunk number 6: 5-advanced-glm.rnw:144-149
###################################################
predict.expFamily <- function(y, b.hat, X, offset = 0) {
  lin.pred <- as.matrix(X) %*% b.hat + offset
  y.hat <- unlink(y, lin.pred)
  return(y.hat)
}


###################################################
### code chunk number 7: 5-advanced-glm.rnw:155-157
###################################################
unlink <- function(y, eta) UseMethod("unlink")
unlink.logit1 <- function(y, eta) 1 / (1 + exp(-eta))


###################################################
### code chunk number 8: 5-advanced-glm.rnw:166-169
###################################################
y <- c(1,0,0,1,1,1,0,1)
X <- as.matrix(cbind(1, 1:8))
beta.hat <- c(0,1)


###################################################
### code chunk number 9: 5-advanced-glm.rnw:177-178
###################################################
class(y) <- c("bernoulli","logit1","expFamily")


###################################################
### code chunk number 10: 5-advanced-glm.rnw:183-184
###################################################
Sjll(beta.hat, X, y)


###################################################
### code chunk number 11: 5-advanced-glm.rnw:206-221
###################################################
maximize <- function(start, f, X, y, offset = 0, ...) {
  optim(par = start,           
        fn = f,
        X = X,
        y = y,
        offset = offset,
        method = "BFGS",
        control = list(
          fnscale = -1,
          reltol = 1e-16,
          maxit = 10000),
        hessian = TRUE,
        ...
        )
}


###################################################
### code chunk number 12: 5-advanced-glm.rnw:229-231
###################################################
test <- maximize(beta.hat, Sjll, X, y)
test$par


###################################################
### code chunk number 13: 5-advanced-glm.rnw:297-299
###################################################
sqrt(2 * dbinom(x = 1, size = 1, prob = 1, log = TRUE) -
         dbinom(x = 1, size = 1, prob = 0.5, log = TRUE))


###################################################
### code chunk number 14: 5-advanced-glm.rnw:307-311
###################################################
devianceResiduals <- function(y, b.hat, X, offset = 0) {
  y.hat <- predict(y, b.hat, X, offset)
  sign(y - y.hat) * sqrt(2 * (jll(y, y) - jll(y, y.hat)))
}


###################################################
### code chunk number 15: 5-advanced-glm.rnw:328-329
###################################################
sum(devianceResiduals(y, test$par, X)^2)


###################################################
### code chunk number 16: 5-advanced-glm.rnw:336-338
###################################################
fit.null <- maximize(0.5, Sjll, 1, y)
sum(devianceResiduals(y, fit.null$par, 1)^2)


###################################################
### code chunk number 17: 5-advanced-glm.rnw:346-347
###################################################
glm(y ~ X[,-1], family=binomial)


###################################################
### code chunk number 18: 5-advanced-glm.rnw:368-373
###################################################
kickStart <- function(y, X, offset)
                             UseMethod("kickStart")
kickStart.default <- function(y, X, offset = 0) {
  coef(lm(I(y - offset) ~ X - 1))
}


###################################################
### code chunk number 19: 5-advanced-glm.rnw:445-446
###################################################
getAnywhere(coef.default)


###################################################
### code chunk number 20: 5-advanced-glm.rnw:464-523
###################################################
ml_glm <- function(formula,
                   data,
                   family,
                   link,
                   offset = 0,
                   start = NULL,
                   verbose = FALSE,
                   ...) {

### Handle the input
  mf <- model.frame(formula, data)
  y <- model.response(mf, "numeric")

### Prepare model infrastructure
  class(y) <- c(family, link, "expFamily")
  X <- model.matrix(formula, data = data)

### Check for missing data.  Stop if any.
  if (any(is.na(cbind(y, X)))) stop("Some data are missing!")

### Initialize the search, if needed
  if (is.null(start))  start <- kickStart(y, X, offset)
  
### Maximize the joint log likelihood
  fit <- maximize(start, Sjll, X, y, offset, ...)

### Check for optim failure and report and stop
  if (verbose | fit$convergence > 0)  print(fit)

### Extract and compute quantities of interest
  beta.hat <- fit$par
  se.beta.hat <- sqrt(diag(solve(-fit$hessian)))
  residuals <- devianceResiduals(y, beta.hat, X, offset, ...)

### Fit null model and determine null deviance 
  fit.null <- maximize(mean(y), Sjll, 1, y, offset, ...)
  null.deviance <-
    sum(devianceResiduals(y, fit.null$par, 1, offset, ...)^2)

### Report the results, with the needs of print.glm in mind
  results <- list(fit = fit,
                  X = X,
                  y = y,
                  call = match.call(),
                  obs = length(y),
                  df.null = length(y) - 1,
                  df.residual = length(y) - length(beta.hat),
                  deviance = sum(residuals^2),    
                  null.deviance = null.deviance, 
                  residuals = residuals,
                  coefficients = beta.hat,
                  se.beta.hat = se.beta.hat,
                  aic = - 2 * fit$val + 2 * length(beta.hat),
                  i = fit$counts[1])

### Use (new) msme class and glm class
  class(results) <- c("msme","glm")
  return(results)
}


###################################################
### code chunk number 21: 5-advanced-glm.rnw:529-531 (eval = FALSE)
###################################################
## library(msme)
## data(medpar)


###################################################
### code chunk number 22: 5-advanced-glm.rnw:535-536
###################################################
load("../../package/msme/data/medpar.rda") 


###################################################
### code chunk number 23: mort.glm
###################################################
mort.glm <- ml_glm(died ~ hmo + white,
                   data = medpar,
                   family = "bernoulli",
                   link = "logit1")


###################################################
### code chunk number 24: 5-advanced-glm.rnw:573-574
###################################################
mort.glm$fit


###################################################
### code chunk number 25: 5-advanced-glm.rnw:586-596
###################################################
residuals.msme <- function(object,
                           type = c("deviance","standard"),
                           ...) {
  type <- match.arg(type)
  if (type == "standard") {
    object$residuals / sqrt(1 - diag(hatvalues(object)))
  } else {
    object$residuals
  }
}


###################################################
### code chunk number 26: 5-advanced-glm.rnw:606-607
###################################################
source("../../package/msme/R/summary.msme.R")


###################################################
### code chunk number 27: 5-advanced-glm.rnw:661-664
###################################################
jll.poisson <- function(y, y.hat, ...) {
  dpois(y, lambda = y.hat, log = TRUE)
}


###################################################
### code chunk number 28: 5-advanced-glm.rnw:670-671
###################################################
unlink.log <- function(y, eta, a=1, m=1) exp(eta)  


###################################################
### code chunk number 29: 5-advanced-glm.rnw:677-680
###################################################
kickStart.log <- function(y, X, offset = 0) {
  coef(lm(I((log(y + 0.1) - offset) ~ X - 1)))
}


###################################################
### code chunk number 30: 5-advanced-glm.rnw:687-691
###################################################
ml.poi <- ml_glm(los ~ hmo + white,
                 family = "poisson",
                 link = "log",
                 data = medpar)


###################################################
### code chunk number 31: 5-advanced-glm.rnw:778-780
###################################################
jll.ztp <- function(y, y.hat, ...) 
  dpois(y, lambda = y.hat, log = TRUE) - log(1 - exp(-y.hat))


###################################################
### code chunk number 32: ml.ztp
###################################################
ml.ztp <- ml_glm(los ~ hmo + white,
                 family = "ztp",
                 link = "log",
                 data = medpar)


###################################################
### code chunk number 33: 5-advanced-glm.rnw:793-794
###################################################
summary(ml.ztp)


###################################################
### code chunk number 34: 5-advanced-glm.rnw:828-834
###################################################
jll.negBinomial <- function(y, y.hat, scale, ...) {
  dnbinom(y,
          mu = y.hat,
          size = 1 / scale,
          log = TRUE)
}


###################################################
### code chunk number 35: 5-advanced-glm.rnw:848-849
###################################################
jll <- function(y, y.hat, scale, ...) UseMethod("jll")


###################################################
### code chunk number 36: 5-advanced-glm.rnw:865-870
###################################################
Sjll2 <- function(b.hat, X, y, p, offset = 0, ...) {
  y.hat <- predict(y, b.hat[1:p], X[,1:p], offset)
  scale.hat <- predict_s(y, b.hat[-(1:p)], X[,-(1:p)])
  sum(jll(y, y.hat, scale.hat, ...))
}


###################################################
### code chunk number 37: 5-advanced-glm.rnw:889-894
###################################################
predict_s <- function(y, b.hat, X) {
  lin.pred <- as.matrix(X) %*% b.hat
  scale <- unlink_s(y, lin.pred)
  return(scale)
}


###################################################
### code chunk number 38: 5-advanced-glm.rnw:908-910
###################################################
unlink_s.log_s <- function(y, eta) exp(eta)
unlink_s <- function(y, eta) UseMethod("unlink_s")


###################################################
### code chunk number 39: 5-advanced-glm.rnw:918-921
###################################################
y <- c(1,0,0,1,1,3,0,9)
X <- as.matrix(cbind(1, 1:8, 1, 1:8))
beta.hat <- c(0,1,0,1)


###################################################
### code chunk number 40: 5-advanced-glm.rnw:926-927
###################################################
class(y) <- c("negBinomial","log","log_s","expFamily")


###################################################
### code chunk number 41: nb-eval
###################################################
Sjll2(beta.hat, X, y, 2)


###################################################
### code chunk number 42: test
###################################################
test <- maximize(beta.hat, Sjll2, X, y, p = 2)
test$par


###################################################
### code chunk number 43: 5-advanced-glm.rnw:987-989
###################################################
getDispersion <- function(y, scale) UseMethod("getDispersion")
getDispersion.negBinomial <- function(y, scale) 1


###################################################
### code chunk number 44: 5-advanced-glm.rnw:992-1000
###################################################
devianceResiduals2 <- function(y, b.hat, X, p, offset = 0) {
  y.hat <- predict(y, b.hat[1:p], X[,1:p], offset)
  scale <- predict_s(y, b.hat[-(1:p)], X[,-(1:p)])
  sign(y - y.hat) *  
    sqrt(2 * getDispersion(y, scale) * 
         (   jll(y, y,     scale) -
             jll(y, y.hat, scale)    )) 
}


###################################################
### code chunk number 45: 5-advanced-glm.rnw:1011-1012
###################################################
sum(devianceResiduals2(y, test$par, X, 2)^2)


###################################################
### code chunk number 46: 5-advanced-glm.rnw:1020-1021
###################################################
library(MASS)


###################################################
### code chunk number 47: our-nb
###################################################
test <- maximize(beta.hat[1:3], Sjll2, X[,1:3], y, p = 2)
test$par


###################################################
### code chunk number 48: 5-advanced-glm.rnw:1053-1054
###################################################
exp(-test$par[3])


###################################################
### code chunk number 49: 5-advanced-glm.rnw:1059-1060
###################################################
sum(devianceResiduals2(y, test$par, X[,1:3], 2)^2)


###################################################
### code chunk number 50: null.eg (eval = FALSE)
###################################################
## #null <- maximize2(beta.hat[c(1,3)], Sjll2, X[,c(1,3)], y, 1)[1]
## null <- maximize(beta.hat[c(1,3)], Sjll2, X[,c(1,3)], y, p = 1)[1]
## sum(devianceResiduals2(y, c(null, test$par[3]), X[,c(1,3)], 1)^2)


###################################################
### code chunk number 51: 5-advanced-glm.rnw:1118-1121
###################################################
kickStart.log <- function(y, X, family, offset = NULL) {
  coef(lm(log(y + 0.01) ~ X - 1, offset = offset))
}


###################################################
### code chunk number 52: 5-advanced-glm.rnw:1141-1216
###################################################
ml_glm2 <- function(formula1,
                    formula2 = ~1, data,
                    family,
                    mean.link,
                    scale.link,
                    m = 1,
                    offset = 0,
                    start = NULL,
                    verbose = FALSE) {

### Handle the input
  mf <- model.frame(formula1, data)
  y <- model.response(mf, "numeric")

### Prepare model infrastructure
  class(y) <- c(family, mean.link, scale.link, "expFamily")
  X1 <- model.matrix(formula1, data = data)
  X2 <- model.matrix(formula2, data = data)
  colnames(X2) <- paste(colnames(X2), "_s", sep="")
  p <- ncol(X1)
  X <- cbind(X1, X2)

### Check for missing data.  Stop if any.
  if (any(is.na(cbind(y, X)))) stop("Some data are missing!")

### Initialize the search 
  if (is.null(start)) {
    start <- c(kickStart(y, X1, offset),
               1, # Shameless hack
               rep(0, ncol(X) - p - 1))
    names(start) <- c(colnames(X1), colnames(X2))
  }

### Maximize the joint log likelihood
  fit <- maximize(start, Sjll2, X, y, offset, p = p)

### Check for optim failure and report and stop
  if (verbose | fit$convergence > 0)  print(fit)

### Extract and compute quantities of interest
  beta.hat <- fit$par
  se.beta.hat <- sqrt(diag(solve(-fit$hessian)))
  residuals <- devianceResiduals2(y, beta.hat, X, p, offset)

#### Deviance residuals for null 
  fit.null <- maximize(c(mean(y), 1), 
                        Sjll2,
                        X[,c(1,p+1)],  y, offset, p = 1)
  null.deviance <-
    sum(devianceResiduals2(y,
                          c(fit.null$par[1], fit$par[p+1]),
                          X[, c(1,p+1)],
                          1,
                          offset)^2)

### Report the results, with the needs of print.glm in mind
  results <- list(fit = fit,
                  loglike = fit$val,
                  X = X,
                  y = y,
                  p = p,
                  call = match.call(),
                  obs = length(y),
                  df.null = length(y) - 2,
                  df.residual = length(y) - length(beta.hat),
                  deviance = sum(residuals^2),    
                  null.deviance = null.deviance, 
                  residuals = residuals,
                  coefficients = beta.hat,
                  se.beta.hat = se.beta.hat,
                  aic = - 2 * fit$val + 2 * length(beta.hat),
                  i = fit$counts[1])
  class(results) <- c("msme","glm")
  return(results)
}


###################################################
### code chunk number 53: Test.3.g
###################################################
test.3.g <- ml_glm2(los ~ hmo + white,
                    formula2 = ~1,
                    data = medpar,
                    family = "negBinomial",
                    mean.link = "log",
                    scale.link = "log_s")


###################################################
### code chunk number 54: 5-advanced-glm.rnw:1292-1294
###################################################
test.3.g$coefficients
exp(0.7244261)


###################################################
### code chunk number 55: 5-advanced-glm.rnw:1300-1301
###################################################
summary(with(test.3.g, predict(y, coefficients[1:p], X[,1:p])))


###################################################
### code chunk number 56: 5-advanced-glm.rnw:1309-1310
###################################################
test.3.g$loglike


###################################################
### code chunk number 57: 5-advanced-glm.rnw:1365-1366
###################################################
coef(test.3.g)[4]


###################################################
### code chunk number 58: 5-advanced-glm.rnw:1370-1371
###################################################
exp(-coef(test.3.g)[4])


###################################################
### code chunk number 59: 5-advanced-glm.rnw:1397-1400
###################################################
unlink_s.log_s <- function(y, eta) exp(eta)
unlink_s.identity_s <- function(y, eta) eta
unlink_s.inverse_s <- function(y, eta) 1 / eta


###################################################
### code chunk number 60: nb-h
###################################################
test.3a.g <- ml_glm2(los ~ hmo + white,
                     formula2 = ~ white + hmo,
                     data = medpar,
                     family = "negBinomial",
                     mean.link = "log",
                     scale.link = "log_s")
summary(test.3a.g, dig = 2)


###################################################
### code chunk number 61: 5-advanced-glm.rnw:1465-1466
###################################################
exp(coef(test.3a.g))


###################################################
### code chunk number 62: AIC
###################################################
test.3.g$aic
test.3a.g$aic


###################################################
### code chunk number 63: 5-advanced-glm.rnw:1484-1485
###################################################
test.3.g$se.beta.hat


###################################################
### code chunk number 64: 5-advanced-glm.rnw:1496-1497
###################################################
with(test.3.g, se.beta.hat * exp(coefficients))


###################################################
### code chunk number 65: 5-advanced-glm.rnw:1508-1515
###################################################
jll.normal <- function(y, y.hat, scale, ...) {
  dnorm(y,
        mean = y.hat,
        sd = scale, log = TRUE)
}
getDispersion.normal <- function(y, scale) scale^2
unlink.identity <- function(y, eta) eta


###################################################
### code chunk number 66: OLS
###################################################
load("../../package/msme/data/ufc.rda")
ufc <- na.omit(ufc)
test.1.g <- ml_glm2(height.m ~ dbh.cm,
                    formula2 = ~1,
                    data = ufc,
                    family = "normal",
                    mean.link = "identity",
                    scale.link = "log_s")


###################################################
### code chunk number 67: 5-advanced-glm.rnw:1560-1561
###################################################
test.1.lm <- lm(height.m ~ dbh.cm, data = ufc)


###################################################
### code chunk number 68: 5-advanced-glm.rnw:1589-1591
###################################################
exp(test.1.g$coefficients[3])
summary(test.1.lm)$sigma


###################################################
### code chunk number 69: GLS
###################################################
test.2.g <- ml_glm2(height.m ~ dbh.cm,
                    formula2 = ~ dbh.cm,
                    data = ufc,
                    family = "normal",
                    mean.link = "identity",
                    scale.link = "log_s")


###################################################
### code chunk number 70: 5-advanced-glm.rnw:1609-1617
###################################################
logLik.msme <- function(object, ...) {
  val <- object$fit$value
  attr(val, "nall") <- nrow(object$X)
  attr(val, "nobs") <- nrow(object$X)
  attr(val, "df") <- length(object$fit$par)
  class(val) <- "logLik"  
  val
}


###################################################
### code chunk number 71: 5-advanced-glm.rnw:1620-1641
###################################################
alrt <- function(x1, x2, ...) {
  jll1 <- logLik(x1)
  jll2 <- logLik(x2)
  df1 <- attr(jll1, "df")
  df2 <- attr(jll2, "df")
  jll.diff <- abs(c(jll1) - c(jll2))
  df.diff <- abs(df1 - df2)
  p.value <- 1 - pchisq(2 * jll.diff, df = df.diff) 
  results <- list(out.tab = data.frame(model = c(1,2),
                                       jll = c(jll1, jll2),
                                       df = c(df1, df2)),
                  jll.diff = jll.diff,
                  df.diff = df.diff,
                  p = p.value)
  cat("\nLL of model 1: ", jll1, " df: ", df1, 
      "\nLL of model 2: ", jll2, " df: ", df2, 
      "\nDifference: ", jll.diff, " df: ", df.diff, 
      "\np-value against H_0: no difference between models ", 
      p.value, "\n")
  return(invisible(results))
}


###################################################
### code chunk number 72: 5-advanced-glm.rnw:1644-1645
###################################################
alrt(test.1.g, test.2.g)


###################################################
### code chunk number 73: 5-advanced-glm.rnw:1698-1707 (eval = FALSE)
###################################################
## 
## jll2.gamma <- function(y, y.hat, scale) {
##   dgamma(y,
##          shape = 1 / scale,
##          scale = y.hat * scale, log = TRUE)
## }
## 
## getDispersion.gamma <- function(y, scale) scale 
## 


###################################################
### code chunk number 74: 5-advanced-glm.rnw:1710-1713 (eval = FALSE)
###################################################
## kickStart.inverse <- function(y, X, family, offset = NULL) {
##   coef(lm(I(1/(y)) ~ X - 1), offset = offset)
## }


###################################################
### code chunk number 75: Test Other Models (eval = FALSE)
###################################################
## 
## library(COUNT)
## 
## data(rwm5yr)
## 
## test.3a.g <- ml_glm2(docvis ~ outwork + female + age + factor(edlevel),
##                      formula2 = ~ outwork + female + age + factor(edlevel),
##                      data = rwm5yr,
##                      family = "negBinomial",
##                      mean.link = "log",
##                      scale.link = "log_s")
## 
## load("~/Dropbox/projects/modelling in R/Chapter 4/rwm1984.RData")
## 
## nb.1 <- ml_glm2(docvis ~ outwork + female + edlevel4,
##                 formula2 = ~ outwork + female + edlevel4,
##                 data = rwm1984,
##                 family = "negBinomial",
##                 mean.link = "log",
##                 scale.link = "log_s")
## summary(nb.1)
## 
## 
## 



### R code from vignette source '6-panel-data.rnw'
### Encoding: ASCII

###################################################
### code chunk number 1: Setup
###################################################
options(repos="http://cran.r-project.org")

if(!require(Hmisc, quietly=TRUE)) install.packages("Hmisc")

rm(list=ls())

options(width = 70)


###################################################
### code chunk number 2: 6-panel-data.rnw:185-187 (eval = FALSE)
###################################################
## library(msme)
## data(medpar)


###################################################
### code chunk number 3: 6-panel-data.rnw:189-190
###################################################
load("../../package/msme/data/medpar.rda")


###################################################
### code chunk number 4: 6-panel-data.rnw:193-199
###################################################
# P__disp.r Calculate Pearson and Pearson dispersion
P__disp <- function(x) {
   pr <- sum(residuals(x, type="pearson")^2)
   dispersion <- pr/x$df.residual
   c(pearson.chi2 = pr, dispersion = dispersion)
}


###################################################
### code chunk number 5: 6-panel-data.rnw:202-207
###################################################
library(MASS)

nbtest1 <- glm.nb(los ~ hmo + white, data = medpar)
P__disp(nbtest1)           # dispersion
1/nbtest1$theta            # alpha     


###################################################
### code chunk number 6: 6-panel-data.rnw:226-231
###################################################

source("../../package/msme/R/ml_glm_2.R")
source("../../package/msme/R/summary.msme.R")
source("../../package/msme/R/residuals.msme.R")



###################################################
### code chunk number 7: nbtest2
###################################################

nbtest2 <- ml_glm2(los ~ hmo + white, 
                    formula2 = ~1,
                    data = medpar,
                    family = "negBinomial",
                    mean.link = "log",
                    scale.link = "inverse_s")


###################################################
### code chunk number 8: 6-panel-data.rnw:284-296
###################################################
pearson2 <- function(object, ...) {
  dispersion <-
    with(object, {
      y.hat <- predict(y, coefficients[1:p], X[,1:p], offset)
      a.hat <- predict_s(y, coefficients[-(1:p)], X[,-(1:p)])
      pearson <- sum(((y - y.hat)^2) /
                     variance(y, y.hat, a.hat))
      pearson / df.residual})
  return(dispersion)
}

pearson2(nbtest2)


###################################################
### code chunk number 9: ufenb.fit
###################################################
medpar$pr <- factor(substr(medpar$provnum, 3, 6)) 

ufenb <- ml_glm2(los ~ hmo + white + pr, 
                 formula2 = ~1,
                 data = medpar,
                 family = "negBinomial",
                 mean.link = "log",
                 scale.link = "identity_s")       


###################################################
### code chunk number 10: disp-ufenb
###################################################
pearson2(ufenb)


###################################################
### code chunk number 11: 6-panel-data.rnw:569-581
###################################################
jll.gNegBinomial <- function(y, y.hat, groups, ...) {
  y.hat.sum.g <- tapply(y.hat, groups, sum)
  y.sum.g <- tapply(y, groups, sum)
  both.sum.g <- y.hat.sum.g + y.sum.g
  lg.y.hat.sum.g <- tapply(lgamma(y.hat), groups, sum)
  lg.y.sum.g <- tapply(lgamma(y + 1), groups, sum)
  lg.both.sum.g <- tapply(lgamma(y + y.hat), groups, sum)
  ll <- sum(lgamma(y.hat.sum.g) + lgamma(y.sum.g + 1) - 
            lgamma(both.sum.g) + lg.both.sum.g -
            lg.y.sum.g - lg.y.hat.sum.g)

}


###################################################
### code chunk number 12: 6-panel-data.rnw:594-595 (eval = FALSE)
###################################################
## library(msme)


###################################################
### code chunk number 13: 6-panel-data.rnw:598-646
###################################################
jll <- function(y, y.hat, ...) UseMethod("jll")

Sjll <- function(b.hat, X, y, offset = 0, ...) {
  y.hat <- predict(y, b.hat, X, offset)
  sum(jll(y, y.hat, ...))
}

predict.expFamily <- function(y, b.hat, X, offset = 0) {
  lin.pred <- as.matrix(X) %*% b.hat + offset
  y.hat <- unlink(y, lin.pred)
  return(y.hat)
}

unlink <- function(y, eta) UseMethod("unlink")
unlink.log <- function(y, eta) exp(eta)

kickStart <- function(y, X, offset)
                             UseMethod("kickStart")
kickStart.default <- function(y, X, offset = 0) {
  coef(lm(I(y - offset) ~ X - 1))
}
kickStart.log <- function(y, X, offset = 0) {
  coef(lm(I(log(y - offset)) ~ X - 1))
}

devianceResiduals <- function(y, b.hat, X, offset = 0, ...) {
  y.hat <- predict(y, b.hat, X, offset)
  sign(y - y.hat) * sqrt(2 * (jll(y, y, ...) - jll(y, y.hat, ...)))
}

maximize <- function(start, f, X, y, offset = 0, ...) {
  optim(start,           
        f,
        X = X,
        y = y,
        offset = offset,
        control = list(
          fnscale = -1,
          reltol = 1e-16,
          maxit = 10000),
        hessian = TRUE,
        ...
        )
}

load("../../package/msme/data/medpar.rda")
source("../../package/msme/R/residuals.msme.R")
source("../../package/msme/R/summary.msme.R")


###################################################
### code chunk number 14: 6-panel-data.rnw:656-692
###################################################
ml_glm3 <- function(formula,
                   data,
                   family,
                   link,
                   offset = 0,
                   start = NULL,
                   verbose = FALSE,
                   ...) {
  mf <- model.frame(formula, data)
  y <- model.response(mf, "numeric")
  class(y) <- c(family, link, "expFamily")
  X <- model.matrix(formula, data = data)
  if (any(is.na(cbind(y, X)))) stop("Some data are missing!")
  if (is.null(start))  start <- kickStart(y, X, offset)
  fit <- maximize(start, Sjll, X, y, offset, ...)
  if (verbose | fit$convergence > 0)  print(fit)
  beta.hat <- fit$par
  se.beta.hat <- sqrt(diag(solve(-fit$hessian)))
  residuals <- devianceResiduals(y, beta.hat, X, offset, ...)
  results <- list(fit = fit,
                  X = X,
                  y = y,
                  call = match.call(),
                  obs = length(y),
                  df.null = length(y) - 1,
                  df.residual = length(y) - length(beta.hat),
                  deviance = sum(residuals^2),    
                  null.deviance = NA,
                  residuals = residuals,
                  coefficients = beta.hat,
                  se.beta.hat = se.beta.hat,
                  aic = - 2 * fit$val + 2 * length(beta.hat),
                  i = fit$counts[1])
  class(results) <- c("msme","glm")
  return(results)
}


###################################################
### code chunk number 15: 6-panel-data.rnw:695-696 (eval = FALSE)
###################################################
## data(medpar)


###################################################
### code chunk number 16: med.nb.g
###################################################
med.nb.g <- ml_glm3(los ~ hmo + white,
                   family = "gNegBinomial", link = "log",
                   group = medpar$provnum, 
                   data = medpar)


###################################################
### code chunk number 17: 6-panel-data.rnw:895-897 (eval = FALSE)
###################################################
## library(msme)
## data(ufc)


###################################################
### code chunk number 18: 6-panel-data.rnw:900-901
###################################################
load("../../package/msme/data/ufc.rda")


###################################################
### code chunk number 19: 6-panel-data.rnw:904-906
###################################################
ufc <- na.omit(ufc)
head(ufc)


###################################################
### code chunk number 20: 6-panel-data.rnw:912-914
###################################################
summary(ufc$height.m)
summary(ufc$dbh.cm)


###################################################
### code chunk number 21: 6-panel-data.rnw:924-929
###################################################
library(nlme)
renorm.lme <- lme(height.m ~ dbh.cm, 
                  random = ~ 1 | plot,
                  data = ufc, method = "ML")
summary(renorm.lme)


###################################################
### code chunk number 22: 6-panel-data.rnw:987-1000
###################################################
jll.gnormal <- function(params, y, X, group, ...) {
   p <- ncol(X)
   N_i <- tapply(y, group, length)
   Su <- exp(params[p+1])
   Se <- exp(params[p+2])
   z <- y - X %*% params[1:p]
   gamma_i <- Su^2 / (N_i * Su^2 + Se^2)
   c1 <- (tapply(z^2, group, sum) -
          gamma_i * tapply(z, group, sum)^2) / Se^2
   c2 <- log(N_i * Su^2 / Se^2 + 1)
   c3 <- N_i * log(2 * pi * Se^2)
   return(sum(-0.5 * (c1 + c2 + c3)))
}


###################################################
### code chunk number 23: 6-panel-data.rnw:1009-1013
###################################################
ufc.model <- with(ufc,
                  list(y = height.m,
                       X = model.matrix(~ dbh.cm),
                       group = plot))


###################################################
### code chunk number 24: 6-panel-data.rnw:1019-1021
###################################################
start <- c(1,1, 1, 1)
with(ufc.model, jll.gnormal(start, y, X, group))


###################################################
### code chunk number 25: 6-panel-data.rnw:1027-1038
###################################################
ufc.fit <-
   with(ufc.model,
        optim(start,
              jll.gnormal,
              y = y,
              X = X,
              group = group,
              method = "BFGS",
              control = list(fnscale = -1,
                             reltol = .Machine$double.eps,
                             maxit = 10000)))


###################################################
### code chunk number 26: 6-panel-data.rnw:1051-1053
###################################################
fixed.effects(renorm.lme)
ufc.fit$par[1:2]


###################################################
### code chunk number 27: 6-panel-data.rnw:1056-1058
###################################################
VarCorr(renorm.lme)[,2]
exp(ufc.fit$par[3:4])


###################################################
### code chunk number 28: 6-panel-data.rnw:1061-1101 (eval = FALSE)
###################################################
## reml.gnormal <- function(params, z, group, ...) {
##    N_i <- tapply(z, group, length)
##    Su <- exp(params[1])
##    Se <- exp(params[2])
##    gamma_i <- Su^2 / (N_i * Su^2 + Se^2)
##    c1 <- (tapply(z^2, group, sum) -
##           gamma_i * tapply(z, group, sum)^2) / Se^2
##    c2 <- log(N_i * Su^2 / Se^2 + 1)
##    c3 <- N_i * log(2 * pi * Se^2)
##    return(sum(-0.5 * (c1 + c2 + c3)))
## }
## 
## ufc.model <- with(ufc,
##                   list(z = residuals(lm(height.m ~ dbh.cm)),
##                        group = plot))
## 
## start <- c(0.5, 2)
## 
## with(ufc.model,
##      reml.gnormal(start,
##                  z, group))
## 
## ufc.fit <-
##    with(ufc.model,
##         optim(start,
##               reml.gnormal,
##               z = z,
##               group = group,
##               control = list(fnscale = -1,
##                              maxit = 10000)))
## 
## ufc.fit
## 
## exp(ufc.fit$par)
## 
## renorm.lme <- lme(height.m ~ dbh.cm, 
##                   random = ~ 1 | plot,
##                   data = ufc)
## 
## VarCorr(summary(renorm.lme))


###################################################
### code chunk number 29: 6-panel-data.rnw:1222-1224 (eval = FALSE)
###################################################
## library(msme)
## data(ufc)


###################################################
### code chunk number 30: 6-panel-data.rnw:1227-1228
###################################################
load("../../package/msme/data/ufc.rda")


###################################################
### code chunk number 31: 6-panel-data.rnw:1236-1237
###################################################
sapply(ufc, function(x) sum(is.na(x)))


###################################################
### code chunk number 32: 6-panel-data.rnw:1242-1243
###################################################
ufc <- ufc[!is.na(ufc$dbh.cm),]


###################################################
### code chunk number 33: 6-panel-data.rnw:1249-1251
###################################################
ufc$na <- is.na(ufc$height.m)
ufc$height.m[ufc$na] <- mean(ufc$height.m, na.rm = TRUE) 


###################################################
### code chunk number 34: 6-panel-data.rnw:1262-1273
###################################################
reps <- 20
trace <- vector(mode = "list", length = reps)

for (i in 1:reps)  {
  refit.lm <- lm(height.m ~ dbh.cm, 
                 data = ufc, 
                 na.action = na.exclude)          # M step
  ufc$height.m[ufc$na] <-
    predict(refit.lm, newdata = ufc)[ufc$na]      # E step
  trace[[i]] <- coef(refit.lm)
}


###################################################
### code chunk number 35: 6-panel-data.rnw:1279-1280
###################################################
trace <- do.call(rbind, trace)


###################################################
### code chunk number 36: fig-em-trace
###################################################
par(las = 1, mar=c(4,4,2,1))
plot(dbh.cm ~ `(Intercept)`, data = trace, type = "b")
base.model <- lm(height.m ~ dbh.cm, data = ufc)
points(coef(base.model)[1], coef(base.model)[2],
       cex = 2, pch = 3)


###################################################
### code chunk number 37: em-trace
###################################################
par(las = 1, mar=c(4,4,2,1))
plot(dbh.cm ~ `(Intercept)`, data = trace, type = "b")
base.model <- lm(height.m ~ dbh.cm, data = ufc)
points(coef(base.model)[1], coef(base.model)[2],
       cex = 2, pch = 3)


###################################################
### code chunk number 38: 6-panel-data.rnw:1358-1360 (eval = FALSE)
###################################################
## data(ufc)
## ufc <- na.omit(ufc)


###################################################
### code chunk number 39: 6-panel-data.rnw:1363-1365
###################################################
load("../../package/msme/data/ufc.rda")
ufc <- na.omit(ufc)


###################################################
### code chunk number 40: EM
###################################################
(N_T <- nrow(ufc))
N <- length(unique(ufc$plot))
y_i <- with(ufc, split(height.m, plot))
x_i <- with(ufc, split(dbh.cm, plot))
X_i <- lapply(x_i, function(x) cbind(1, x))
Z_i <- lapply(x_i, function(x) matrix(1, nrow = length(x)))


###################################################
### code chunk number 41: 6-panel-data.rnw:1382-1390
###################################################
lm.start <- with(ufc, lm(height.m ~ dbh.cm))
(s2 <- summary(lm.start)$sigma^2)
(D <- with(ufc,
          var(tapply(residuals(lm.start),
                     plot,
                     mean))))
beta <- coef(lm.start)
dim(beta) <- c(length(beta), 1)


###################################################
### code chunk number 42: 6-panel-data.rnw:1395-1399
###################################################
V_inv <- lapply(Z_i,
                function(Z, D)
                solve(diag(nrow(Z)) + Z %*% D %*% t(Z)),
                D = D) 


###################################################
### code chunk number 43: 6-panel-data.rnw:1405-1468
###################################################
for (i in 1:2000) {

## First, compute the list of E_i as within-group residuals 
## from the fixed effects. 

  E_i <- mapply(function(y, X, beta)
                y - X %*% beta,
                y = y_i, X = X_i, beta = list(beta))

## Then obtain the within-group variance contribution by the 
## weighted inner product of the residuals for each group.

  s2.i <- mapply(function(E, V) 
                 sum(t(E) %*% V %*% E),
                 E = E_i, V = V_inv)
  
## Sum these across the groups, scale them, and use them to 
## update the variance estimate.

  s2 <- s2 - 1 + sum(s2.i) / (s2 * N_T) # This is 6.10

## D_i is the groupwise calculated correction to the matrix D.
## Refer to Equation 6.11; the following is the quantity within
## the squared brackets.

  D.i <-
    mapply(function(D, Z, V, E, s2) 
           D %*% t(Z) %*% V %*% Z %*% D -
           D %*% t(Z) %*% V %*% E %*% t(E) %*% V %*% Z %*% D / 
             s2,
           D = list(D), Z = Z_i, V = V_inv, E = E_i, 
           s2 = list(s2))
  
  D <- D - sum(D.i) / N                  # This is 6.11
  
## We can now update V, given the new D.  Note our use of 
## diag(nrow(Z)) to create a suitably sized identity matrix
## on the fly. 

  V_inv <- lapply(Z_i,
                  function(Z, D)
                  solve(diag(nrow(Z)) + Z %*% D %*% t(Z)),
                  D = D)                 # This is 6.9

## Finally we update Beta

  beta.denom <-
    mapply(function(X, D, Z)
           t(X) %*% solve(diag(nrow(Z)) + 
              Z %*% D %*% t(Z)) %*% X, 
           X = X_i, D = list(D), Z = Z_i, SIMPLIFY = FALSE)

  beta.num <-
    mapply(function(X, D, Z, y)
           t(X) %*% solve(diag(nrow(Z)) + 
              Z %*% D %*% t(Z)) %*% y,
           X = X_i, D = list(D), Z = Z_i, y = y_i, 
           SIMPLIFY = FALSE)

  beta <- solve(Reduce(`+`, beta.denom)) %*% Reduce(`+`, beta.num)
          # 6.12

}


###################################################
### code chunk number 44: 6-panel-data.rnw:1478-1480
###################################################
fixed.effects(renorm.lme)
beta[,1]


###################################################
### code chunk number 45: 6-panel-data.rnw:1484-1486
###################################################
VarCorr(renorm.lme)[,1]
c(D*s2, s2)



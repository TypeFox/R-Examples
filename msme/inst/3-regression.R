### R code from vignette source '3-regression.rnw'

###################################################
### code chunk number 1: Setup
###################################################
options(repos="http://cran.r-project.org")

if(!require(Hmisc, quietly=TRUE)) install.packages("Hmisc")
if(!require(xtable, quietly=TRUE)) install.packages("xtable")

rm( list=ls() )

require(ggplot2)
theme_set(theme_bw())

require(gridExtra)

none <- theme_blank()

theme_update(axis.line = theme_segment(colour = "black", size = 1.05),
             panel.grid.major = none,
             panel.border = none,
             panel.grid.minor = none)

vplayout <- function(x,y) {
  viewport(layout.pos.row = x, layout.pos.col = y)
}

options(width = 67)


###################################################
### code chunk number 2: 3-regression.rnw:121-123
###################################################
y <- c(3, 5, 7)
x <- c(1, 2, 3)


###################################################
### code chunk number 3: 3-regression.rnw:132-135
###################################################
least.squares <- function(p, x, y) {
  sum((y - (p[1] + p[2] * x))^2)
}


###################################################
### code chunk number 4: 3-regression.rnw:147-149
###################################################
start.searching.here <- c(intercept = 0, slope = 0)
least.squares(start.searching.here, x, y)


###################################################
### code chunk number 5: 3-regression.rnw:159-162
###################################################
optim(par = start.searching.here, 
      fn = least.squares, 
      x = x, y = y)$par


###################################################
### code chunk number 6: 3-regression.rnw:191-193
###################################################
mean(y) - cov(x,y) / var(x) * mean(x)   # Beta 0
cov(x,y) / var(x)                       # Beta 1


###################################################
### code chunk number 7: 3-regression.rnw:214-219
###################################################
L1.obj <- function(p, x, y) {
  sum(abs(y - (p[1] + p[2] * x)))
}

optim(c(0,0), L1.obj, x=x, y=y)$par


###################################################
### code chunk number 8: 3-regression.rnw:224-229
###################################################
Linf.obj <- function(p, x, y) {
  max(abs(y - (p[1] + p[2] * x)))
}

optim(start.searching.here, Linf.obj, x=x, y=y)$par


###################################################
### code chunk number 9: 3-regression.rnw:367-369
###################################################
X <- as.matrix(cbind(1, x))
(beta.hat <- solve(t(X) %*% X) %*% t(X) %*% y)


###################################################
### code chunk number 10: 3-regression.rnw:412-415
###################################################
(y.hat <- X %*% beta.hat)
(sigma.2 <- as.numeric(var(y - y.hat)))
(vcov.beta.hat <- sigma.2 * solve(crossprod(X)))


###################################################
### code chunk number 11: 3-regression.rnw:451-452
###################################################
X.hat <- X %*% solve(t(X) %*% X) %*% t(X)


###################################################
### code chunk number 12: 3-regression.rnw:457-459
###################################################
(var.y.hat <- sigma.2 * diag(X.hat))
(var.e.hat <- sigma.2 * (1 - diag(X.hat)))


###################################################
### code chunk number 13: QR intro
###################################################
(xR <- qr.R(qr(X)))
(xQ <- qr.Q(qr(X)))
xQ %*% xR


###################################################
### code chunk number 14: 3-regression.rnw:556-558
###################################################
Y <- matrix(c(3,5,7), nrow=3)
xR


###################################################
### code chunk number 15: 3-regression.rnw:562-563
###################################################
t(xQ) %*% Y


###################################################
### code chunk number 16: 3-regression.rnw:571-572
###################################################
backsolve(xR, t(xQ) %*% Y)


###################################################
### code chunk number 17: 3-regression.rnw:585-586
###################################################
(hat.values <- diag(crossprod(t(xQ))))


###################################################
### code chunk number 18: 3-regression.rnw:628-630 (eval = FALSE)
###################################################
## library(msme)
## data(ufc)


###################################################
### code chunk number 19: 3-regression.rnw:633-634
###################################################
load("../../package/msme/data/ufc.rda")


###################################################
### code chunk number 20: 3-regression.rnw:640-641
###################################################
ufc <- na.omit(ufc)


###################################################
### code chunk number 21: 3-regression.rnw:650-652
###################################################
X <- cbind(1, ufc$dbh.cm)
Y <- ufc$height.m


###################################################
### code chunk number 22: More QR
###################################################
xR <- qr.R(qr(X))
xQ <- qr.Q(qr(X))


###################################################
### code chunk number 23: 3-regression.rnw:664-668
###################################################
(beta.hat <- backsolve(xR, t(xQ) %*% Y))
y.hat <- X %*% beta.hat
(sigma.2 <- as.numeric(var(Y - y.hat)))
(vcov.beta.hat <- sigma.2 * solve(crossprod(xR)))


###################################################
### code chunk number 24: 3-regression.rnw:671-672
###################################################
sqrt(diag(vcov.beta.hat))


###################################################
### code chunk number 25: 3-regression.rnw:686-689
###################################################
e.hat <- Y - y.hat
e.hat.ss <- e.hat / 
   sqrt(sigma.2 * (1 - diag(crossprod(t(xQ)))))


###################################################
### code chunk number 26: 3-regression.rnw:693-695
###################################################
var(e.hat)
var(e.hat.ss)


###################################################
### code chunk number 27: fig-lm-diag-1
###################################################
qplot(y.hat, e.hat, alpha = I(0.2),
      ylab = "Residuals", xlab = "Fitted Values") + 
  geom_abline(intercept = 0, slope = 0) +
  geom_smooth(aes(x = y.hat, y = e.hat), se = FALSE)  


###################################################
### code chunk number 28: fig-lm-diag-2
###################################################
qplot(y.hat, abs(sqrt(e.hat.ss)), alpha = I(0.2),
      ylab = "Standardized Studentized Residuals",
      xlab = "Fitted Values") + 
  geom_smooth(aes(x = y.hat, y = abs(sqrt(e.hat.ss))),
      se = FALSE)


###################################################
### code chunk number 29: lm-diag-2
###################################################
grid.newpage()
pushViewport(viewport(layout = grid.layout(1,2)))
print(
qplot(y.hat, e.hat, alpha = I(0.2),
      ylab = "Residuals", xlab = "Fitted Values") + 
  geom_abline(intercept = 0, slope = 0) +
  geom_smooth(aes(x = y.hat, y = e.hat), se = FALSE)  
, vp = vplayout(1,1))
print(
qplot(y.hat, abs(sqrt(e.hat.ss)), alpha = I(0.2),
      ylab = "Standardized Studentized Residuals",
      xlab = "Fitted Values") + 
  geom_smooth(aes(x = y.hat, y = abs(sqrt(e.hat.ss))),
      se = FALSE)
, vp = vplayout(1,2))


###################################################
### code chunk number 30: 3-regression.rnw:816-819
###################################################
jll.normal <- function(p, x, y) {
  sum(dnorm(y, p[1] + p[2] * x, p[3], log = TRUE))
}


###################################################
### code chunk number 31: 3-regression.rnw:831-835
###################################################
optim(par = c(intercept = 0, slope = 0, sigma = 1),
      fn = jll.normal,
      control = list(fnscale = -1),
      x = x, y = y)$par


###################################################
### code chunk number 32: 3-regression.rnw:901-902
###################################################
options()$contrasts


###################################################
### code chunk number 33: 3-regression.rnw:924-931
###################################################
jll.normal <- function(params, X, y) {
  p <- length(params)
  beta <- params[-p]
  sigma <- params[p]
  linpred <- X %*% beta  
  sum(dnorm(y, mean = linpred, sd = sigma, log = TRUE))
}


###################################################
### code chunk number 34: 3-regression.rnw:981-988
###################################################
jll.normal <- function(params, X, y) {
  p <- length(params)
  beta <- params[-p]
  sigma <- exp(params[p])
  linpred <- X %*% beta  
  sum(dnorm(y, mean = linpred, sd = sigma, log = TRUE))
}


###################################################
### code chunk number 35: ml_g
###################################################
ml_g <- function(formula, data, ...) {

### Prepare the data, relying on the formula class for 
### handling model specification
  mf <- model.frame(formula, data)
  y <- model.response(mf, "numeric")
  X <- model.matrix(formula, data = data)

### Check for missing data.  Stop if any.
  if (any(is.na(cbind(y, X)))) stop("Some data are missing.")

### Declare the joint log likelihood function 
  jll.normal <- function(params, X, y) {
     p <- length(params)
     beta <- params[-p]
     sigma <- exp(params[p])
     linpred <- X %*% beta  
     sum(dnorm(y, mean = linpred, sd = sigma, log = TRUE))
     }

### Initialize the search
  ls.reg <- lm(y ~ X - 1)
  beta.hat.ls <- coef(ls.reg)
  sigma.hat.ls <- sd(residuals(ls.reg))
  start <- c(beta.hat.ls, sigma.hat.ls)

### Maximize the joint log likelihood
  fit <- optim(start,           
               jll.normal,
               X = X,
               y = y,
               control = list(
                 fnscale = -1,
                 maxit = 10000),
               hessian = TRUE
               )

### Check for optim failure and report and stop
  if (fit$convergence > 0) {
    print(fit)
    stop("optim failed to converge!")
  }

### Post-processing 
  beta.hat <- fit$par
  se.beta.hat <- sqrt(diag(solve(-fit$hessian)))

### Reporting
  results <- list(fit = fit,
                  X = X,
                  y = y,
                  call = match.call(),
                  beta.hat = beta.hat,
                  se.beta.hat = se.beta.hat,
                  sigma.hat = exp(beta.hat[length(beta.hat)]))

### Prepare for S3 deployment (see next Section!)
  class(results) <- c("ml_g_fit","lm")
  return(results)
}


###################################################
### code chunk number 36: 3-regression.rnw:1189-1190
###################################################
ufc.g.reg <- ml_g(height.m ~ dbh.cm, data = ufc)


###################################################
### code chunk number 37: 3-regression.rnw:1258-1261
###################################################
coef.ml_g_fit <- function(x, ...) {
  x$beta.hat[-length(x$beta.hat)]
}


###################################################
### code chunk number 38: 3-regression.rnw:1266-1267
###################################################
coef(ufc.g.reg)


###################################################
### code chunk number 39: 3-regression.rnw:1275-1276
###################################################
ufc.g.reg


###################################################
### code chunk number 40: 3-regression.rnw:1289-1292
###################################################
fitted.ml_g_fit <- function(x, ...) {
   as.numeric(x$X %*% coef(x))
}


###################################################
### code chunk number 41: 3-regression.rnw:1303-1304
###################################################
str(fitted(ufc.g.reg))


###################################################
### code chunk number 42: hatvalues
###################################################
hatvalues.ml_g_fit <- function(x, ...) {
  tcrossprod(qr.Q(qr(x$X)))
}


###################################################
### code chunk number 43: 3-regression.rnw:1334-1335
###################################################
str(hatvalues(ufc.g.reg))


###################################################
### code chunk number 44: 3-regression.rnw:1349-1359
###################################################
residuals.ml_g_fit <-
   function(x, type = c("raw","ss"), ...) {
   type <- match.arg(type)
   e.hat <- x$y - fitted(x)
   if (type == "ss") {
      e.hat <- e.hat /
              (x$sigma.hat * sqrt(1 - diag(hatvalues(x)))) 
    }
   return(e.hat)
}


###################################################
### code chunk number 45: 3-regression.rnw:1371-1374
###################################################
ufc.g.res.s <- residuals(ufc.g.reg, type = "ss")
str(ufc.g.res.s)
var(ufc.g.res.s)


###################################################
### code chunk number 46: 3-regression.rnw:1449-1476
###################################################
plot.ml_g_fit <- function(x, ...) {
  require(ggplot2)
  require(gridExtra)
  e.hat <- residuals(x)
  e.hat.ss <- residuals(x, type="ss")
  y.hat <- fitted(x)
  n <- nrow(x$X)
  pp1 <- qplot(y.hat, x$y, alpha = I(0.2),        ### Plot 1
            ylab = "Observations", xlab = "Fitted Values") + 
         geom_abline(intercept = 0, slope = 1) +
         geom_smooth(aes(x = y.hat, y = x$y), se = FALSE)
  pp2 <- qplot(y.hat, e.hat, alpha = I(0.2),      ### Plot 2
               ylab = "Residuals", xlab = "Fitted Values") + 
         geom_abline(intercept = 0, slope = 0) +
         geom_smooth(aes(x = y.hat, y = e.hat), se = FALSE)
  pp3 <- qplot(y.hat, abs(sqrt(e.hat.ss)),         ### Plot 3
               alpha = I(0.2),
               ylab = "Sqrt (Abs( Stand. Res.))",
               xlab = "Fitted Values") + 
         geom_smooth(aes(x = y.hat, y = abs(sqrt(e.hat.ss))),
               se = FALSE)
  pp4 <- qplot(sort(e.hat.ss), qnorm((1:n)/(n+1)), ### Plot 4
               alpha = I(0.2),
               xlab = "Standardized Studentized Residuals",
               ylab = "Normal Quantiles")
  grid.arrange(pp1, pp2, pp3, pp4, ncol=2)
}


###################################################
### code chunk number 47: fig-ufc-lm-diag
###################################################
plot(ufc.g.reg)


###################################################
### code chunk number 48: ufc-lm-diag
###################################################
plot(ufc.g.reg)


###################################################
### code chunk number 49: 3-regression.rnw:1528-1536
###################################################
logLik.ml_g_fit <- function(x, ...) {
  val <- x$fit$value
  attr(val, "nall") <- nrow(x$X)
  attr(val, "nobs") <- nrow(x$X)
  attr(val, "df") <- length(x$fit$par)
  class(val) <- "logLik"  
  val
}


###################################################
### code chunk number 50: 3-regression.rnw:1540-1541
###################################################
logLik(ufc.g.reg)


###################################################
### code chunk number 51: 3-regression.rnw:1558-1559
###################################################
AIC(ufc.g.reg)


###################################################
### code chunk number 52: 3-regression.rnw:1575-1595
###################################################
summary.ml_g_fit <- function(object, dig = 3, ...) {
  zTable <- with(object, 
              data.frame(Estimate = beta.hat,
                         SE = se.beta.hat,
                         Z = beta.hat / se.beta.hat,
                         LCL = beta.hat - 1.96 * se.beta.hat,
                         UCL = beta.hat + 1.96 * se.beta.hat))
  rownames(zTable) <- c(colnames(object$X), "Log Sigma")
  p <- length(object$fit$par)
  n <- nrow(object$X)
  df <- c(p, n, p)
  summ <- list(call = object$call,
               coefficients = zTable,
               df = df, 
               residuals = residuals(object),
               aliased = rep(FALSE, p),
               sigma = object$sigma.hat)
  class(summ) <- c("summary.ml_g_fit", "summary.lm")
  return(summ)
}


###################################################
### code chunk number 53: 3-regression.rnw:1599-1600
###################################################
summary(ufc.g.reg)


###################################################
### code chunk number 54: 3-regression.rnw:1606-1607
###################################################
methods(class = "ml_g_fit")


###################################################
### code chunk number 55: 3-regression.rnw:1625-1626
###################################################
length(methods(class = "lm"))


###################################################
### code chunk number 56: 3-regression.rnw:1647-1648
###################################################
print(coef(summary(lm(height.m ~ dbh.cm, ufc))), digits = 3)


###################################################
### code chunk number 57: 3-regression.rnw:1659-1660
###################################################
library(splines)


###################################################
### code chunk number 58: 3-regression.rnw:1666-1667
###################################################
ufc.g.spline <- ml_g(height.m ~ bs(dbh.cm),  data = ufc)


###################################################
### code chunk number 59: fig-ufc-sp-diag
###################################################
plot(ufc.g.spline)


###################################################
### code chunk number 60: ufc-sp-diag
###################################################
plot(ufc.g.spline)


###################################################
### code chunk number 61: 3-regression.rnw:1695-1697
###################################################
logLik(ufc.g.reg)
logLik(ufc.g.spline)


###################################################
### code chunk number 62: 3-regression.rnw:1704-1725
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
### code chunk number 63: 3-regression.rnw:1747-1748
###################################################
alrt(ufc.g.reg, ufc.g.spline)


###################################################
### code chunk number 64: 3-regression.rnw:1754-1756
###################################################
ufc.g.quad <- ml_g(height.m ~ dbh.cm + I(dbh.cm^2),  data = ufc)
alrt(ufc.g.quad, ufc.g.spline)


###################################################
### code chunk number 65: 3-regression.rnw:1809-1812
###################################################
y <- matrix(c(10,9,7,4,8,12,11,7,3,5,3,12,9,10), ncol=1)
X <- matrix(c(1,4,6,0,5,4,1,4,7,1,2,2,1,4,
              6,0,7,5,1,5,1,1,2,5,0,3,4,4),ncol=2)


###################################################
### code chunk number 66: 3-regression.rnw:1833-1859 (eval = FALSE)
###################################################
## ## Works
##   
## ufc.trans <- ml_g(I(log(height.m)) ~ I(1 / (dbh.cm + 1)), 
##                   data = na.omit(ufc))
## plot(ufc.trans)
## 
## ## Works
## 
## ufc.poly <- ml_g(height.m ~ dbh.cm + I(dbh.cm^2), 
##                  data = na.omit(ufc))
## plot(ufc.poly)
## 
## library(splines)
## 
## ## Doesn't work yet
##             
##   jll.normal <- function(params, X, y) {
##     p <- length(params)
##     beta <- params[-p]
##     log.sigma <- params[p]
##     linpred <- X %*% beta  # Form the linear predictor
##     sum(dnorm((y - linpred) / exp(log.sigma),
##                mean = 0,
##                sd = 1,
##                log = TRUE) - log.sigma)
##   }



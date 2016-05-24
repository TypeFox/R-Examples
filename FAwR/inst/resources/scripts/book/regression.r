### R code from vignette source 'regression.rnw'

###################################################
### code chunk number 1: Setup
###################################################
options(width=60)
options(repos="http://cran.r-project.org")
#if(!require(systemfit, quietly=TRUE)) 
#  install.packages("systemfit", dependencies=TRUE)
#if(!require(nlreg, quietly=TRUE)) 
#  install.packages("nlreg", dependencies=TRUE)
if(!require(ggplot2, quietly=TRUE)) {
  install.packages("ggplot2", dependencies=TRUE)
  require(ggplot2, quietly=TRUE)
}
if(!require(nlreg, quietly=TRUE)) {
  install.packages("nlreg", dependencies=TRUE)
  require(nlreg, quietly=TRUE)
}
require(MASS, quietly=TRUE)
require(lattice, quietly=TRUE)
require(nlme, quietly=TRUE)
require(mgcv, quietly=TRUE)


###################################################
### code chunk number 2: sweetgum
###################################################
Stangle("../../ch2/sweave/sweetgum.rnw")
source("sweetgum.R")


###################################################
### code chunk number 3: fig-sg-plot
###################################################
par(las = 1)
plot(vol.m3 ~ dbh.cm, 
     data = sweetgum, 
     xlab = "Diameter (cm)", 
     ylab = expression(paste("Volume (", m^3, ")")))


###################################################
### code chunk number 4: sg-plot
###################################################
par(las = 1)
plot(vol.m3 ~ dbh.cm, 
     data = sweetgum, 
     xlab = "Diameter (cm)", 
     ylab = expression(paste("Volume (", m^3, ")")))


###################################################
### code chunk number 5: linear regression 1
###################################################
sweetgum.lm.d <- lm(vol.m3 ~ dbh.cm, data=sweetgum)


###################################################
### code chunk number 6: fig-diag-sg-lm1
###################################################
par(mfrow = c(2, 2), mar=c(4, 4, 3, 1), las = 1)
plot(sweetgum.lm.d)


###################################################
### code chunk number 7: diag-sg-lm1
###################################################
par(mfrow = c(2, 2), mar=c(4, 4, 3, 1), las = 1)
plot(sweetgum.lm.d)


###################################################
### code chunk number 8: log transform
###################################################
sweetgum$log.vol.m3 <- log(sweetgum$vol.m3)
sweetgum$log.dbh.cm <- log(sweetgum$dbh.cm)
sweetgum.lm.ld <- lm(log.vol.m3 ~ log.dbh.cm, 
                     data = sweetgum)


###################################################
### code chunk number 9: fig-diag-sg-lm2
###################################################
par(mfrow = c(2, 2), mar=c(4, 4, 3, 1), las=1)
plot(sweetgum.lm.ld)


###################################################
### code chunk number 10: diag-sg-lm2
###################################################
par(mfrow = c(2, 2), mar=c(4, 4, 3, 1), las=1)
plot(sweetgum.lm.ld)


###################################################
### code chunk number 11: capture
###################################################
summ.ld <- capture.output(summary(sweetgum.lm.ld))


###################################################
### code chunk number 12: regression.rnw:535-536 (eval = FALSE)
###################################################
## summary(sweetgum.lm.ld)


###################################################
### code chunk number 13: capture-a
###################################################
cat(summ.ld[1:3], sep="\n")


###################################################
### code chunk number 14: capture-b
###################################################
cat(summ.ld[5:7], sep="\n")


###################################################
### code chunk number 15: capture-c
###################################################
cat(summ.ld[9:14], sep="\n")


###################################################
### code chunk number 16: capture-d
###################################################
cat(summ.ld[16:18], sep="\n")


###################################################
### code chunk number 17: regression.rnw:592-596
###################################################
reg.table <- as.data.frame(coef(summary(sweetgum.lm.ld)))
reg.table$Estimate[2] + 
  reg.table$`Std. Error`[2] * 
  qt(c(0.025, 0.975), summary(sweetgum.lm.ld)$df[2])


###################################################
### code chunk number 18: regression.rnw:602-603
###################################################
confint(sweetgum.lm.ld)


###################################################
### code chunk number 19: regression.rnw:632-637
###################################################
library(gmodels)
estimable(sweetgum.lm.ld, 
          rbind(Intercept = c(1,0),  
                Slope     = c(0,1)), 
          conf.int = 0.95)[,c(1,2,4,6,7)]


###################################################
### code chunk number 20: regression.rnw:648-649
###################################################
methods(estimable)


###################################################
### code chunk number 21: regression.rnw:721-722
###################################################
exp(summary(sweetgum.lm.ld)$sigma^2 / 2)


###################################################
### code chunk number 22: regression.rnw:739-752
###################################################
sweetgum.vol.hat <- function(dbh.cm, 
                             ht.dbh.lm = sweetgum.lm.ld, 
                             correct = TRUE) {
   sweetgum.hat <- data.frame(dbh.cm = dbh.cm)
   sweetgum.hat$log.dbh.cm <- log(sweetgum.hat$dbh.cm)
   correction.2 <- ifelse(correct, 
                          exp(summary(ht.dbh.lm)$sigma^2 / 2),
                          1)
   sweetgum.hat$vol.m3 <- 
     exp(predict(ht.dbh.lm, newdata = sweetgum.hat)) * 
       correction.2
   return(sweetgum.hat$vol.m3)
}


###################################################
### code chunk number 23: regression.rnw:756-757
###################################################
sweetgum.vol.hat(10, sweetgum.lm.ld)


###################################################
### code chunk number 24: regression.rnw:764-765
###################################################
sweetgum.vol.hat(c(10,12), sweetgum.lm.ld)


###################################################
### code chunk number 25: regression.rnw:769-770
###################################################
sweetgum.vol.hat(c(10,12), sweetgum.lm.ld, correct = FALSE)


###################################################
### code chunk number 26: fig-lm-use
###################################################
par(las = 1)
plot(vol.m3 ~ dbh.cm, 
     data = sweetgum, 
     xlab = "Diameter (cm)", 
     ylab = expression(paste("Volume (", m^3, ")")))
curve(sweetgum.vol.hat, from=0, to=100, add=TRUE)


###################################################
### code chunk number 27: lm-use
###################################################
par(mar=c(5,5,0,5))
par(las = 1)
plot(vol.m3 ~ dbh.cm, 
     data = sweetgum, 
     xlab = "Diameter (cm)", 
     ylab = expression(paste("Volume (", m^3, ")")))
curve(sweetgum.vol.hat, from=0, to=100, add=TRUE)


###################################################
### code chunk number 28: regression.rnw:817-820
###################################################
exp(predict(sweetgum.lm.ld, 
            newdata = data.frame(log.dbh.cm = log(c(10, 20))),
            interval = "prediction"))


###################################################
### code chunk number 29: regression.rnw:837-839
###################################################
names(sweetgum.lm.ld)
names(summary(sweetgum.lm.ld))


###################################################
### code chunk number 30: regression.rnw:855-857
###################################################
class(sweetgum.lm.ld)
head(methods(class = class(sweetgum.lm.ld)))


###################################################
### code chunk number 31: regression.rnw:864-866
###################################################
sweetgum.lm.ld$call
summary(sweetgum.lm.ld)$sigma


###################################################
### code chunk number 32: regression.rnw:873-880
###################################################
sigma <- function(x) {
      if (class(x) == "lm") {
        return(summary(x)$sigma)
      } else {
        stop("Object is not a linear model (class 'lm').")
      }
} 


###################################################
### code chunk number 33: regression.rnw:883-884
###################################################
sigma(sweetgum.lm.ld)


###################################################
### code chunk number 34: regression.rnw:919-920
###################################################
summary(sweetgum.lm.ld)$coef


###################################################
### code chunk number 35: ufc
###################################################
Stangle("../../ch2/sweave/ufc.rnw")
source("ufc.R")


###################################################
### code chunk number 36: regression.rnw:942-946
###################################################
hd.lm.1 <- lm(I(log(height.m)) ~ dbh.cm * species, 
              data = ufc.tree, 
              subset = height.m > 0)
anova(hd.lm.1)


###################################################
### code chunk number 37: regression.rnw:966-970
###################################################
hd.lm.2 <- lm(I(log(height.m)) ~ dbh.cm, 
              data = ufc.tree, 
              subset = height.m > 0)
anova(hd.lm.1, hd.lm.2)


###################################################
### code chunk number 38: regression.rnw:985-988
###################################################
require(MASS)
data(npk)
anova(lm(yield ~ block + N * P + K, npk))


###################################################
### code chunk number 39: regression.rnw:996-998
###################################################
anova(lm(terms(yield ~ block + N * P + K, 
               keep.order = TRUE), npk))


###################################################
### code chunk number 40: logging
###################################################
ufc.tree$log.height.m <- log(ufc.tree$height.m)
hd.lm.4a <- lm(log.height.m ~ dbh.cm * species, 
               data = ufc.tree, 
               subset = height.m > 0)
hd.lm.4b <- lm(I(log(height.m)) ~ dbh.cm * species, 
               data = ufc.tree, 
               subset=height.m > 0)


###################################################
### code chunk number 41: regression.rnw:1052-1055
###################################################
hd.lm.5 <- lm(height.m ~ dbh.cm * species, 
              weights = dbh.cm^-2, 
              data = ufc.tree) 


###################################################
### code chunk number 42: regression.rnw:1067-1073
###################################################
unweighted <- 
  coef(lm(height.m ~ dbh.cm * species - 1 - dbh.cm, 
          data = ufc.tree))
weighted <- 
  coef(lm(height.m ~ dbh.cm * species - 1 - dbh.cm, 
          weights = dbh.cm^-2, data = ufc.tree))


###################################################
### code chunk number 43: fig-weighted-comp
###################################################
intercepts <- 1:10; slopes <- 11:20
par(mfrow = c(1,2), las = 1, mar = c(4,4,3,1))
plot(unweighted[intercepts], weighted[intercepts], 
     main = "Intercepts", type = "n", 
     xlab = "Unweighted", ylab = "Weighted")
abline(0, 1, col="darkgrey")
text(unweighted[intercepts], weighted[intercepts], 
     levels(ufc.tree$species))
plot(unweighted[slopes], weighted[slopes],
     main = "Slopes", type = "n", 
     xlab = "Unweighted", ylab = "Weighted")
abline(0, 1, col = "darkgrey")
text(unweighted[slopes], weighted[slopes], 
     levels(ufc.tree$species))


###################################################
### code chunk number 44: weighted-comp
###################################################
intercepts <- 1:10; slopes <- 11:20
par(mfrow = c(1,2), las = 1, mar = c(4,4,3,1))
plot(unweighted[intercepts], weighted[intercepts], 
     main = "Intercepts", type = "n", 
     xlab = "Unweighted", ylab = "Weighted")
abline(0, 1, col="darkgrey")
text(unweighted[intercepts], weighted[intercepts], 
     levels(ufc.tree$species))
plot(unweighted[slopes], weighted[slopes],
     main = "Slopes", type = "n", 
     xlab = "Unweighted", ylab = "Weighted")
abline(0, 1, col = "darkgrey")
text(unweighted[slopes], weighted[slopes], 
     levels(ufc.tree$species))


###################################################
### code chunk number 45: gls
###################################################
library(nlme)
sweetgum.gls.ld <- gls(log.vol.m3 ~ log.dbh.cm, 
                       weights = varPower(form = ~ dbh.cm),
                       data = sweetgum)


###################################################
### code chunk number 46: regression.rnw:1161-1162
###################################################
summary(sweetgum.gls.ld)


###################################################
### code chunk number 47: regression.rnw:1176-1177
###################################################
anova(sweetgum.gls.ld, sweetgum.lm.ld)


###################################################
### code chunk number 48: regression.rnw:1186-1187 (eval = FALSE)
###################################################
## methods(class = "gls")


###################################################
### code chunk number 49: get.gutten
###################################################
Stangle("../../ch2/sweave/gutten.rnw")
source("gutten.R")


###################################################
### code chunk number 50: regression.rnw:1285-1286
###################################################
library(ggplot2)


###################################################
### code chunk number 51: fig-gut-plot
###################################################
qplot(x = age.bh, y = dbh.cm, group = tree,
      linetype = site, facets = ~ location,
      xlab = "Age (y)", ylab = "Dbh (cm)",
      geom = "line", 
      data = gutten) +
  scale_x_continuous(breaks = 40 * (0:3)) +
  scale_y_continuous(breaks = 10 * (0:5))


###################################################
### code chunk number 52: gut-plot
###################################################
print(
qplot(x = age.bh, y = dbh.cm, group = tree,
      linetype = site, facets = ~ location,
      xlab = "Age (y)", ylab = "Dbh (cm)",
      geom = "line", 
      data = gutten) +
  scale_x_continuous(breaks = 40 * (0:3)) +
  scale_y_continuous(breaks = 10 * (0:5))
      )


###################################################
### code chunk number 53: model
###################################################
dbh.growth <- 
  deriv(~ asymptote * (1 - exp(-log(2)/scale * x)),
        c("asymptote","scale"),        
        function(x, asymptote, scale){},
        hessian = TRUE)


###################################################
### code chunk number 54: regression.rnw:1391-1392
###################################################
handy.tree <- subset(gutten, tree.ID == "1.1")


###################################################
### code chunk number 55: regression.rnw:1401-1402
###################################################
max(handy.tree$dbh.cm, na.rm=TRUE)


###################################################
### code chunk number 56: regression.rnw:1414-1419
###################################################
handy.nls <- 
  nls(dbh.cm ~ dbh.growth(age.bh, asymptote, scale),
      start = list(asymptote = 29, scale = 10),
      na.action = na.exclude,
      data = handy.tree)


###################################################
### code chunk number 57: regression.rnw:1476-1477
###################################################
rms.curv(handy.nls) 


###################################################
### code chunk number 58: regression.rnw:1499-1507
###################################################
dbh.growth.PB <-deriv(~ asymptote * (1 - exp(-exp(scale) * x)),
                    c("asymptote","scale"),
                    function(x, asymptote, scale){},
                    hessian = TRUE)
handy.nls.PB <- nls(dbh.cm ~ dbh.growth.PB(age.bh, asymptote, scale),
                   start = list(asymptote=29, scale=-3.5),
                   data = handy.tree)
rms.curv(handy.nls.PB) 


###################################################
### code chunk number 59: regression.rnw:1525-1526
###################################################
mean(residuals(handy.nls), na.rm=TRUE)


###################################################
### code chunk number 60: fig-nls-resid
###################################################
par(mfrow=c(1,3), mar=c(4,4,2,1), las=1)
plot(fitted(handy.nls), residuals(handy.nls, type="pearson"), 
     xlab = "Fitted Values", ylab = "Standardized Residuals")
abline(h=0, col="red")
qqnorm(residuals(handy.nls, type="pearson"))
qqline(residuals(handy.nls, type="pearson"))
plot(fitted(handy.nls), handy.tree$dbh.cm, 
     xlab = "Fitted Values", ylab = "Observed Values")
abline(0, 1, col="red")


###################################################
### code chunk number 61: nls-resid
###################################################
par(mfrow=c(1,3), mar=c(4,4,2,1), las=1)
plot(fitted(handy.nls), residuals(handy.nls, type="pearson"), 
     xlab = "Fitted Values", ylab = "Standardized Residuals")
abline(h=0, col="red")
qqnorm(residuals(handy.nls, type="pearson"))
qqline(residuals(handy.nls, type="pearson"))
plot(fitted(handy.nls), handy.tree$dbh.cm, 
     xlab = "Fitted Values", ylab = "Observed Values")
abline(0, 1, col="red")


###################################################
### code chunk number 62: fig-profile
###################################################
handy.prof <- profile(handy.nls)
opar <- par(mfrow=c(1,2), mar=c(4,4,2,1), las=1)
plot(handy.prof, conf = c(0.95))


###################################################
### code chunk number 63: profile
###################################################
handy.prof <- profile(handy.nls)
opar <- par(mfrow=c(1,2), mar=c(4,4,2,1), las=1)
plot(handy.prof, conf = c(0.95))


###################################################
### code chunk number 64: fig-fig-contours
###################################################
pairs(handy.prof)


###################################################
### code chunk number 65: fig-contours
###################################################
pairs(handy.prof)


###################################################
### code chunk number 66: regression.rnw:1640-1643
###################################################
(cov.hat.mat <- 
   matrix(summary(handy.nls)$cov.unscaled, nrow=2) * 
   summary(handy.nls)$sigma^2)


###################################################
### code chunk number 67: regression.rnw:1647-1648
###################################################
cov.hat <- cov.hat.mat[1,2]


###################################################
### code chunk number 68: regression.rnw:1653-1654
###################################################
cov.hat / prod(summary(handy.nls)$coefficients[,2])


###################################################
### code chunk number 69: regression.rnw:1668-1669
###################################################
summary(handy.nls)


###################################################
### code chunk number 70: regression.rnw:1679-1682
###################################################
(my.t <- qt(0.975, summary(handy.nls)$df[2]))
coef(summary(handy.nls))[,1:2] %*% 
  matrix(c(1,-my.t,1,my.t), nrow=2)


###################################################
### code chunk number 71: confint
###################################################
confint(handy.prof)


###################################################
### code chunk number 72: regression.rnw:1702-1704
###################################################
sqrt(sum(residuals(handy.nls)^2, na.rm=TRUE) / 
     summary(handy.nls)$df[2])


###################################################
### code chunk number 73: regression.rnw:1708-1709
###################################################
sd(handy.tree$dbh.cm, na.rm=TRUE)


###################################################
### code chunk number 74: handy.dbh.hat
###################################################
handy.dbh.hat <- function(age.bh) 
  predict(handy.nls, newdata = data.frame(age.bh = age.bh))


###################################################
### code chunk number 75: fig-handy-use
###################################################
par(las = 1)
plot(dbh.cm ~ age.bh, data = handy.tree, 
     xlim = c(0, max(handy.tree$age.bh, na.rm=TRUE)),
     ylim = c(0, max(handy.tree$dbh.cm, na.rm=TRUE)),
     ylab = "Diameter (cm)", xlab = "Age (y)")
curve(handy.dbh.hat, add = TRUE)


###################################################
### code chunk number 76: handy-use
###################################################
par(las = 1)
plot(dbh.cm ~ age.bh, data = handy.tree, 
     xlim = c(0, max(handy.tree$age.bh, na.rm=TRUE)),
     ylim = c(0, max(handy.tree$dbh.cm, na.rm=TRUE)),
     ylab = "Diameter (cm)", xlab = "Age (y)")
curve(handy.dbh.hat, add = TRUE)


###################################################
### code chunk number 77: handier.nls
###################################################
handier.nls <- nls(dbh.cm ~ SSasymp(age.bh, Asym, R0, lrc),
                 na.action = na.exclude,
                 data = handy.tree)


###################################################
### code chunk number 78: regression.rnw:1781-1782
###################################################
anova(handy.nls, handier.nls)


###################################################
### code chunk number 79: regression.rnw:1794-1795
###################################################
confint(handier.nls)


###################################################
### code chunk number 80: apropos
###################################################
apropos("^SS")


###################################################
### code chunk number 81: regression.rnw:1891-1895
###################################################
apropos("^SS")[sapply(apropos("^SS"), 
                      function(x) {
                        "selfStart" %in% 
                        class(eval(parse(text = x)))}) ]


###################################################
### code chunk number 82: regression.rnw:1899-1900
###################################################
length(apropos("."))


###################################################
### code chunk number 83: regression.rnw:1914-1917
###################################################
allometric <- deriv(~ alpha * x^beta, c("alpha", "beta"), 
                    function(x, alpha, beta){}, 
                    hessian = TRUE)


###################################################
### code chunk number 84: regression.rnw:1921-1931
###################################################
allometric.init <- function (mCall, data, LHS) {
  xy <- data.frame(sortedXyData(mCall[["x"]], LHS, data))
  if (nrow(xy) < 3) 
    stop("Too few observations to fit allometric function")
  pars <- as.vector(coef(lm(I(log(y)) ~ I(log(x)), 
                            data = xy)))
  pars[1] <- exp(pars[1])
  names(pars) <- mCall[c("alpha", "beta")]
  return(pars)
}


###################################################
### code chunk number 85: regression.rnw:1935-1938
###################################################
SSallometric <- selfStart(allometric,
                          allometric.init,
                          c("alpha", "beta"))


###################################################
### code chunk number 86: SSallometric
###################################################
nls(vol.m3 ~ SSallometric(dbh.cm, alpha, beta), 
    data = sweetgum)


###################################################
### code chunk number 87: regression.rnw:1977-1983
###################################################
normal.ll <- function(parameters, x, y) {
  sum(dnorm(y, 
            parameters[1] + parameters[2] * x, 
            parameters[3], 
            log = TRUE))
}


###################################################
### code chunk number 88: regression.rnw:1988-1994
###################################################
good.fit <- optim(c(intercept = 1, slope = 1, sigma = 1), 
                  normal.ll, 
                  hessian = TRUE, 
                  control = list(fnscale = -1),
                  x = sweetgum$log.dbh.cm, 
                  y = sweetgum$log.vol.m3)


###################################################
### code chunk number 89: regression.rnw:1999-2000
###################################################
good.fit$par


###################################################
### code chunk number 90: regression.rnw:2004-2005
###################################################
sqrt(diag(solve(-good.fit$hessian)))


###################################################
### code chunk number 91: regression.rnw:2026-2032
###################################################
normal.ll.nl <- function(parameters, x, y) {
  sum( dnorm(y, 
             parameters[1] * x ^ parameters[2], 
             parameters[3], 
             log = TRUE ))
}


###################################################
### code chunk number 92: regression.rnw:2036-2042
###################################################
good.fit <- optim(c(intercept = 1, slope = 1, sigma = 1), 
                  normal.ll.nl, 
                  hessian = TRUE, 
                  control = list(fnscale = -1),
                  x = sweetgum$dbh.cm, 
                  y = sweetgum$vol.m3)


###################################################
### code chunk number 93: regression.rnw:2046-2047
###################################################
good.fit$par


###################################################
### code chunk number 94: regression.rnw:2051-2052
###################################################
sqrt(diag(solve(-good.fit$hessian)))


###################################################
### code chunk number 95: regression.rnw:2086-2092
###################################################
t3.ll <- function(parameters, x, y) {
  sum(dt((y - parameters[1] - x * parameters[2]) / 
         exp(parameters[3]), 
         df = 10, 
         log = TRUE) - parameters[3])
}


###################################################
### code chunk number 96: regression.rnw:2096-2102
###################################################
good.fit.t <- optim(c(intercept = 1, slope = 1, sigma = 1), 
                    t3.ll, 
                    hessian = TRUE, 
                    control = list(fnscale = -1),
                    x = sweetgum$log.dbh.cm, 
                    y = sweetgum$log.vol.m3)


###################################################
### code chunk number 97: regression.rnw:2106-2107
###################################################
good.fit.t$par


###################################################
### code chunk number 98: regression.rnw:2111-2112
###################################################
sqrt(diag(solve(-good.fit.t$hessian)))


###################################################
### code chunk number 99: regression.rnw:2122-2123
###################################################
exp(good.fit.t$par[3]) * sqrt(10 / (10 - 2))


###################################################
### code chunk number 100: regression.rnw:2209-2211
###################################################
system("rm -fr package-ChA")
package.skeleton(name = "package-ChA")



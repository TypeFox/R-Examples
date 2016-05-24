## ----knitup, include = FALSE---------------------------------------------
library(knitr)
require(printr, quietly = TRUE)
set.seed(225452)

## ----eval = FALSE--------------------------------------------------------
#  devtools::install_github("sachsmc/pseval")

## ----setup---------------------------------------------------------------
library(pseval)
library(survival)

## ------------------------------------------------------------------------
fakedata <- generate_example_data(n = 500)
head(fakedata)

## ----psdesign------------------------------------------------------------
binary.ps <- psdesign(data = fakedata, Z = Z, Y = Y.obs, S = S.obs, BIP = BIP)
binary.ps

## ----casecont------------------------------------------------------------
fakedata.cc <- fakedata
missdex <- sample((1:nrow(fakedata.cc))[fakedata.cc$Y.obs == 0], 
       size = floor(sum(fakedata.cc$Y.obs == 0) * .8))
fakedata.cc[missdex, ]$S.obs <- NA
fakedata.cc$weights <- ifelse(fakedata.cc$Y.obs == 1, 1, .2)

## ----casecont2-----------------------------------------------------------
binary.cc <- psdesign(data = fakedata.cc, Z = Z, Y = Y.obs, S = S.obs, BIP = BIP, weights = weights)
binary.cc

## ----helpimp-------------------------------------------------------------
?add_integration

## ----impp----------------------------------------------------------------
binary.ps <- binary.ps + integrate_parametric(S.1 ~ BIP)
binary.ps

## ----imp2----------------------------------------------------------------
binary.ps + integrate_parametric(S.0 ~ BIP)

## ----imp3----------------------------------------------------------------
library(splines)
binary.ps + integrate_parametric(S.1 ~ BIP^2)
binary.ps + integrate_parametric(S.1 ~ bs(BIP, df = 3))

## ----riskhelp------------------------------------------------------------
?add_riskmodel

## ----riskbin-------------------------------------------------------------
binary.ps <- binary.ps + risk_binary(model = Y ~ S.1 * Z, D = 50, risk = risk.logit)
binary.ps

## ----est-----------------------------------------------------------------
binary.est <- binary.ps + ps_estimate(method = "BFGS")
binary.boot <- binary.est + ps_bootstrap(n.boots = 50, progress.bar = FALSE, 
                            start = binary.est$estimates$par, method = "BFGS")
binary.boot

## ----alltog, eval = FALSE------------------------------------------------
#  binary.est <- psdesign(data = fakedata, Z = Z, Y = Y.obs, S = S.obs, BIP = BIP) +
#    integrate_parametric(S.1 ~ BIP) +
#    risk_binary(model = Y ~ S.1 * Z, D = 50, risk = risk.logit) +
#    ps_estimate(method = "BFGS")

## ----summary-------------------------------------------------------------
smary <- summary(binary.boot)

## ----calcrisk------------------------------------------------------------
head(calc_risk(binary.boot, contrast = "VE", n.samps = 20))
head(calc_risk(binary.boot, contrast = function(R0, R1) 1 - R1/R0, n.samps = 20))

## ----plot1---------------------------------------------------------------
plot(binary.boot, contrast = "VE", lwd = 2)
abline(h = smary$VE.estimates[2], lty = 2)

expit <- function(x) exp(x)/(1 + exp(x))
trueVE <- function(s){
  
  r0 <- expit(-1 - 0 * s)
  r1 <- expit(-1 - .75 * s)
  1 - r1/r0
  
}

rug(binary.boot$augdata$S.1)
curve(trueVE(x), add = TRUE, col = "red")
legend("bottomright", legend = c("estimated VE", "marginal VE", "true VE"), 
       col = c("black", "black", "red"), lty = c(1, 3, 1), lwd = c(2, 1, 1))

## ----logrr---------------------------------------------------------------
plot(binary.boot, contrast = "logRR", lwd = 2)
plot(binary.boot, contrast = "RR", log = "y", lwd = 2)

## ----ve------------------------------------------------------------------
ve.est <- calc_risk(binary.boot, CI.type = "pointwise", n.samps = 200)
head(ve.est)

## ----plotciag------------------------------------------------------------
plot(binary.boot, contrast = "VE", lwd = 2, CI.type = "band")
sbs <- calc_risk(binary.boot, CI.type = "pointwise", n.samps = 200)
lines(Y.lower.CL.2.5 ~ S.1, data = sbs, lty = 3, lwd = 2)
lines(Y.upper.CL.97.5 ~ S.1, data = sbs, lty = 3, lwd = 2)
legend("bottomright", lwd = 2, lty = 1:3, legend = c("estimate", "simultaneous CI", "pointwise CI"))

## ----ggpt----------------------------------------------------------------
library(ggplot2)
VE.est <- calc_risk(binary.boot, n.samps = 200)
ggplot(VE.est, 
       aes(x = S.1, y = Y, ymin = Y.lower.CL.0.95, ymax = Y.upper.CL.0.95)) + 
  geom_line() + geom_ribbon(alpha = .2) + ylab(attr(VE.est, "Y.function"))

## ----ccest---------------------------------------------------------------
cc.fit <- binary.cc + integrate_parametric(S.1 ~ BIP) + 
  risk_binary(D = 10) + ps_estimate()
cc.fit

## ----surv1---------------------------------------------------------------
surv.fit <- psdesign(fakedata, Z = Z, Y = Surv(time.obs, event.obs), 
                     S = S.obs, BIP = BIP, CPV = CPV) + 
  integrate_semiparametric(formula.location = S.1 ~ BIP, formula.scale = S.1 ~ 1) + 
  risk_exponential(D = 10) + ps_estimate(method = "BFGS") + ps_bootstrap(n.boots = 20)
surv.fit
plot(surv.fit)

## ------------------------------------------------------------------------
with(fakedata, table(S.obs.cat, BIP.cat))

## ----catfit--------------------------------------------------------------
cat.fit <- psdesign(fakedata, Z = Z, Y = Y.obs, 
                     S = S.obs.cat, BIP = BIP.cat) + 
  integrate_nonparametric(formula = S.1 ~ BIP) + 
  risk_binary(Y ~ S.1 * Z, D = 10, risk = risk.probit) + ps_estimate(method = "BFGS")
cat.fit
plot(cat.fit)

## ----catfitps------------------------------------------------------------
cat.fit.ps <- psdesign(fakedata, Z = Z, Y = Y.obs, 
                     S = S.obs, BIP = BIP.cat) + 
  integrate_nonparametric(formula = S.1 ~ BIP) + 
  risk_binary(Y ~ S.1 * Z, D = 10, risk = risk.logit) + ps_estimate(method = "pseudo-score") + 
  ps_bootstrap(n.boots = 20, method = "pseudo-score")
summary(cat.fit.ps)
plot(cat.fit.ps)


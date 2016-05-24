###################################################
### chunk number 1: setup
###################################################
options(prompt = "R> ", continue = "+  ", width = 64,
  digits = 4, show.signif.stars = FALSE, useFancyQuotes = FALSE)

options(SweaveHooks = list(onefig =   function() {par(mfrow = c(1,1))},
                           twofig =   function() {par(mfrow = c(1,2))},                           
                           threefig = function() {par(mfrow = c(1,3))},
                           fourfig =  function() {par(mfrow = c(2,2))},
			   sixfig =   function() {par(mfrow = c(3,2))}))

library("AER")

set.seed(1071)


###################################################
### chunk number 2: data-journals
###################################################
data("Journals")
journals <- Journals[, c("subs", "price")]
journals$citeprice <- Journals$price/Journals$citations
summary(journals)


###################################################
### chunk number 3: linreg-plot eval=FALSE
###################################################
## plot(log(subs) ~ log(citeprice), data = journals)
## jour_lm <- lm(log(subs) ~ log(citeprice), data = journals)
## abline(jour_lm)


###################################################
### chunk number 4: linreg-plot1
###################################################
plot(log(subs) ~ log(citeprice), data = journals)
jour_lm <- lm(log(subs) ~ log(citeprice), data = journals)
abline(jour_lm)


###################################################
### chunk number 5: linreg-class
###################################################
class(jour_lm)


###################################################
### chunk number 6: linreg-names
###################################################
names(jour_lm)


###################################################
### chunk number 7: linreg-summary
###################################################
summary(jour_lm)


###################################################
### chunk number 8: linreg-summary
###################################################
jour_slm <- summary(jour_lm)
class(jour_slm)
names(jour_slm)


###################################################
### chunk number 9: linreg-coef
###################################################
jour_slm$coefficients


###################################################
### chunk number 10: linreg-anova
###################################################
anova(jour_lm)


###################################################
### chunk number 11: journals-coef
###################################################
coef(jour_lm)


###################################################
### chunk number 12: journals-confint
###################################################
confint(jour_lm, level = 0.95)


###################################################
### chunk number 13: journals-predict
###################################################
predict(jour_lm, newdata = data.frame(citeprice = 2.11),
  interval = "confidence")
predict(jour_lm, newdata = data.frame(citeprice = 2.11), 
  interval = "prediction")


###################################################
### chunk number 14: predict-plot eval=FALSE
###################################################
## lciteprice <- seq(from = -6, to = 4, by = 0.25)
## jour_pred <- predict(jour_lm, interval = "prediction",
##   newdata = data.frame(citeprice = exp(lciteprice)))  
## plot(log(subs) ~ log(citeprice), data = journals)
## lines(jour_pred[, 1] ~ lciteprice, col = 1)    
## lines(jour_pred[, 2] ~ lciteprice, col = 1, lty = 2)
## lines(jour_pred[, 3] ~ lciteprice, col = 1, lty = 2)


###################################################
### chunk number 15: predict-plot1
###################################################
lciteprice <- seq(from = -6, to = 4, by = 0.25)
jour_pred <- predict(jour_lm, interval = "prediction",
  newdata = data.frame(citeprice = exp(lciteprice)))  
plot(log(subs) ~ log(citeprice), data = journals)
lines(jour_pred[, 1] ~ lciteprice, col = 1)    
lines(jour_pred[, 2] ~ lciteprice, col = 1, lty = 2)
lines(jour_pred[, 3] ~ lciteprice, col = 1, lty = 2)


###################################################
### chunk number 16: journals-plot eval=FALSE
###################################################
## par(mfrow = c(2, 2))
## plot(jour_lm)
## par(mfrow = c(1, 1))


###################################################
### chunk number 17: journals-plot1
###################################################
par(mfrow = c(2, 2))
plot(jour_lm)
par(mfrow = c(1, 1))


###################################################
### chunk number 18: journal-lht
###################################################
linearHypothesis(jour_lm, "log(citeprice) = -0.5")


###################################################
### chunk number 19: CPS-data
###################################################
data("CPS1988")
summary(CPS1988)


###################################################
### chunk number 20: CPS-base
###################################################
cps_lm <- lm(log(wage) ~ experience + I(experience^2) +
  education + ethnicity, data = CPS1988)


###################################################
### chunk number 21: CPS-visualization-unused eval=FALSE
###################################################
## ex <- 0:56
## ed <- with(CPS1988, tapply(education, 
##   list(ethnicity, experience), mean))[, as.character(ex)]
## fm <- cps_lm
## wago <- predict(fm, newdata = data.frame(experience = ex, 
##   ethnicity = "cauc", education = as.numeric(ed["cauc",])))
## wagb <- predict(fm, newdata = data.frame(experience = ex, 
##   ethnicity = "afam", education = as.numeric(ed["afam",])))
## plot(log(wage) ~ experience, data = CPS1988, pch = ".", 
##   col = as.numeric(ethnicity))
## lines(ex, wago)
## lines(ex, wagb, col = 2)


###################################################
### chunk number 22: CPS-summary
###################################################
summary(cps_lm)


###################################################
### chunk number 23: CPS-noeth
###################################################
cps_noeth <- lm(log(wage) ~ experience + I(experience^2) +
  education, data = CPS1988)
anova(cps_noeth, cps_lm)


###################################################
### chunk number 24: CPS-anova
###################################################
anova(cps_lm)


###################################################
### chunk number 25: CPS-noeth2 eval=FALSE
###################################################
## cps_noeth <- update(cps_lm, formula = . ~ . - ethnicity)


###################################################
### chunk number 26: CPS-waldtest
###################################################
waldtest(cps_lm, . ~ . - ethnicity)


###################################################
### chunk number 27: CPS-spline
###################################################
cps_plm <- lm(log(wage) ~ bs(experience, df = 5) +
  education + ethnicity, data = CPS1988)


###################################################
### chunk number 28: CPS-spline-summary eval=FALSE
###################################################
## summary(cps_plm)


###################################################
### chunk number 29: CPS-BIC
###################################################
cps_bs <- lapply(3:10, function(i) lm(log(wage) ~
  bs(experience, df = i) + education + ethnicity,
  data = CPS1988))
structure(sapply(cps_bs, AIC, k = log(nrow(CPS1988))),
  .Names = 3:10)


###################################################
### chunk number 30: plm-plot eval=FALSE
###################################################
## cps <- data.frame(experience = -2:60, education =
##   with(CPS1988, mean(education[ethnicity == "cauc"])),
##   ethnicity = "cauc")
## cps$yhat1 <- predict(cps_lm, newdata = cps)
## cps$yhat2 <- predict(cps_plm, newdata = cps)
## 
## plot(log(wage) ~ jitter(experience, factor = 3), pch = 19,
##   col = rgb(0.5, 0.5, 0.5, alpha = 0.02), data = CPS1988)
## lines(yhat1 ~ experience, data = cps, lty = 2)
## lines(yhat2 ~ experience, data = cps)
## legend("topleft", c("quadratic", "spline"), lty = c(2,1),
##   bty = "n")


###################################################
### chunk number 31: plm-plot1
###################################################
cps <- data.frame(experience = -2:60, education =
  with(CPS1988, mean(education[ethnicity == "cauc"])),
  ethnicity = "cauc")
cps$yhat1 <- predict(cps_lm, newdata = cps)
cps$yhat2 <- predict(cps_plm, newdata = cps)

plot(log(wage) ~ jitter(experience, factor = 3), pch = 19,
  col = rgb(0.5, 0.5, 0.5, alpha = 0.02), data = CPS1988)
lines(yhat1 ~ experience, data = cps, lty = 2)
lines(yhat2 ~ experience, data = cps)
legend("topleft", c("quadratic", "spline"), lty = c(2,1),
  bty = "n")


###################################################
### chunk number 32: CPS-int
###################################################
cps_int <- lm(log(wage) ~ experience + I(experience^2) +
  education * ethnicity, data = CPS1988)
coeftest(cps_int)


###################################################
### chunk number 33: CPS-int2 eval=FALSE
###################################################
## cps_int <- lm(log(wage) ~ experience + I(experience^2) +
##   education + ethnicity + education:ethnicity,
##   data = CPS1988)


###################################################
### chunk number 34: CPS-sep
###################################################
cps_sep <- lm(log(wage) ~ ethnicity /
  (experience + I(experience^2) + education) - 1,
  data = CPS1988)


###################################################
### chunk number 35: CPS-sep-coef
###################################################
cps_sep_cf <- matrix(coef(cps_sep), nrow = 2)
rownames(cps_sep_cf) <- levels(CPS1988$ethnicity)
colnames(cps_sep_cf) <- names(coef(cps_lm))[1:4]
cps_sep_cf


###################################################
### chunk number 36: CPS-sep-anova
###################################################
anova(cps_sep, cps_lm)


###################################################
### chunk number 37: CPS-sep-visualization-unused eval=FALSE
###################################################
## ex <- 0:56
## ed <- with(CPS1988, tapply(education, list(ethnicity, 
##   experience), mean))[, as.character(ex)]
## fm <- cps_lm
## wago <- predict(fm, newdata = data.frame(experience = ex, 
##   ethnicity = "cauc", education = as.numeric(ed["cauc",])))
## wagb <- predict(fm, newdata = data.frame(experience = ex, 
##   ethnicity = "afam", education = as.numeric(ed["afam",])))
## plot(log(wage) ~ jitter(experience, factor = 2), 
##   data = CPS1988, pch = ".", col = as.numeric(ethnicity))
## 
## 
## plot(log(wage) ~ as.factor(experience), data = CPS1988, 
##   pch = ".")
## lines(ex, wago, lwd = 2)
## lines(ex, wagb, col = 2, lwd = 2)
## fm <- cps_sep
## wago <- predict(fm, newdata = data.frame(experience = ex, 
##   ethnicity = "cauc", education = as.numeric(ed["cauc",])))
## wagb <- predict(fm, newdata = data.frame(experience = ex, 
##   ethnicity = "afam", education = as.numeric(ed["afam",])))
## lines(ex, wago, lty = 2, lwd = 2)
## lines(ex, wagb, col = 2, lty = 2, lwd = 2)


###################################################
### chunk number 38: CPS-region
###################################################
CPS1988$region <- relevel(CPS1988$region, ref = "south")
cps_region <- lm(log(wage) ~ ethnicity + education +
  experience + I(experience^2) + region, data = CPS1988)
coef(cps_region)


###################################################
### chunk number 39: wls1
###################################################
jour_wls1 <- lm(log(subs) ~ log(citeprice), data = journals,
  weights = 1/citeprice^2)


###################################################
### chunk number 40: wls2
###################################################
jour_wls2 <- lm(log(subs) ~ log(citeprice), data = journals,
  weights = 1/citeprice)


###################################################
### chunk number 41: journals-wls1 eval=FALSE
###################################################
## plot(log(subs) ~ log(citeprice), data = journals)
## abline(jour_lm)
## abline(jour_wls1, lwd = 2, lty = 2)
## abline(jour_wls2, lwd = 2, lty = 3)
## legend("bottomleft", c("OLS", "WLS1", "WLS2"),
##   lty = 1:3, lwd = 2, bty = "n")


###################################################
### chunk number 42: journals-wls11
###################################################
plot(log(subs) ~ log(citeprice), data = journals)
abline(jour_lm)
abline(jour_wls1, lwd = 2, lty = 2)
abline(jour_wls2, lwd = 2, lty = 3)
legend("bottomleft", c("OLS", "WLS1", "WLS2"),
  lty = 1:3, lwd = 2, bty = "n")


###################################################
### chunk number 43: fgls1
###################################################
auxreg <- lm(log(residuals(jour_lm)^2) ~ log(citeprice),
  data = journals)
jour_fgls1 <- lm(log(subs) ~ log(citeprice), 
  weights = 1/exp(fitted(auxreg)), data = journals)


###################################################
### chunk number 44: fgls2
###################################################
gamma2i <- coef(auxreg)[2]
gamma2 <- 0
while(abs((gamma2i - gamma2)/gamma2) > 1e-7) {
  gamma2 <- gamma2i
  fglsi <- lm(log(subs) ~ log(citeprice), data = journals, 
    weights = 1/citeprice^gamma2)
  gamma2i <- coef(lm(log(residuals(fglsi)^2) ~
    log(citeprice), data = journals))[2]
}
jour_fgls2 <- lm(log(subs) ~ log(citeprice), data = journals,
  weights = 1/citeprice^gamma2)


###################################################
### chunk number 45: fgls2-coef
###################################################
coef(jour_fgls2)


###################################################
### chunk number 46: journals-fgls
###################################################
plot(log(subs) ~ log(citeprice), data = journals)
abline(jour_lm)
abline(jour_fgls2, lty = 2, lwd = 2)


###################################################
### chunk number 47: usmacro-plot eval=FALSE
###################################################
## data("USMacroG")
## plot(USMacroG[, c("dpi", "consumption")], lty = c(3, 1),
##   plot.type = "single", ylab = "")
## legend("topleft", legend = c("income", "consumption"),
##   lty = c(3, 1), bty = "n")


###################################################
### chunk number 48: usmacro-plot1
###################################################
data("USMacroG")
plot(USMacroG[, c("dpi", "consumption")], lty = c(3, 1),
  plot.type = "single", ylab = "")
legend("topleft", legend = c("income", "consumption"),
  lty = c(3, 1), bty = "n")


###################################################
### chunk number 49: usmacro-fit
###################################################
library("dynlm")
cons_lm1 <- dynlm(consumption ~ dpi + L(dpi), data = USMacroG)
cons_lm2 <- dynlm(consumption ~ dpi + L(consumption), 
  data = USMacroG)


###################################################
### chunk number 50: usmacro-summary1
###################################################
summary(cons_lm1)


###################################################
### chunk number 51: usmacro-summary2
###################################################
summary(cons_lm2)


###################################################
### chunk number 52: dynlm-plot eval=FALSE
###################################################
## plot(merge(as.zoo(USMacroG[,"consumption"]), fitted(cons_lm1),
##   fitted(cons_lm2), 0, residuals(cons_lm1),
##   residuals(cons_lm2)), screens = rep(1:2, c(3, 3)),
##   lty = rep(1:3, 2), ylab = c("Fitted values", "Residuals"),
##   xlab = "Time", main = "")
## legend(0.05, 0.95, c("observed", "cons_lm1", "cons_lm2"), 
##   lty = 1:3, bty = "n")


###################################################
### chunk number 53: dynlm-plot1
###################################################
plot(merge(as.zoo(USMacroG[,"consumption"]), fitted(cons_lm1),
  fitted(cons_lm2), 0, residuals(cons_lm1),
  residuals(cons_lm2)), screens = rep(1:2, c(3, 3)),
  lty = rep(1:3, 2), ylab = c("Fitted values", "Residuals"),
  xlab = "Time", main = "")
legend(0.05, 0.95, c("observed", "cons_lm1", "cons_lm2"), 
  lty = 1:3, bty = "n")


###################################################
### chunk number 54: encompassing1
###################################################
cons_lmE <- dynlm(consumption ~ dpi + L(dpi) +
  L(consumption), data = USMacroG)


###################################################
### chunk number 55: encompassing2
###################################################
anova(cons_lm1, cons_lmE, cons_lm2)


###################################################
### chunk number 56: encompassing3
###################################################
encomptest(cons_lm1, cons_lm2)


###################################################
### chunk number 57: plm-data
###################################################
data("Grunfeld", package = "AER")
library("plm")
gr <- subset(Grunfeld, firm %in% c("General Electric",
  "General Motors", "IBM"))
pgr <- plm.data(gr, index = c("firm", "year"))


###################################################
### chunk number 58: plm-pool
###################################################
gr_pool <- plm(invest ~ value + capital, data = pgr, 
  model = "pooling")


###################################################
### chunk number 59: plm-FE
###################################################
gr_fe <- plm(invest ~ value + capital, data = pgr, 
  model = "within")
summary(gr_fe)


###################################################
### chunk number 60: plm-pFtest
###################################################
pFtest(gr_fe, gr_pool)


###################################################
### chunk number 61: plm-RE
###################################################
gr_re <- plm(invest ~ value + capital, data = pgr, 
  model = "random", random.method = "walhus")
summary(gr_re)


###################################################
### chunk number 62: plm-plmtest
###################################################
plmtest(gr_pool)


###################################################
### chunk number 63: plm-phtest
###################################################
phtest(gr_re, gr_fe)


###################################################
### chunk number 64: EmplUK-data
###################################################
data("EmplUK", package = "plm")
form <- log(emp) ~ log(wage) + log(capital) + log(output)


###################################################
### chunk number 65: plm-AB
###################################################
empl_ab <- pgmm(dynformula(form, list(2, 1, 0, 1)),
  data = EmplUK, index = c("firm", "year"),
  effect = "twoways", model = "twosteps",
  gmm.inst = ~ log(emp), lag.gmm = list(c(2, 99)))


###################################################
### chunk number 66: plm-AB-summary
###################################################
summary(empl_ab)     


###################################################
### chunk number 67: systemfit
###################################################
library("systemfit")
gr2 <- subset(Grunfeld, firm %in% c("Chrysler", "IBM"))
pgr2 <- plm.data(gr2, c("firm", "year"))


###################################################
### chunk number 68: SUR
###################################################
gr_sur <- systemfit(invest ~ value + capital,
  method = "SUR", data = pgr2)
summary(gr_sur, residCov = FALSE, equations = FALSE)


###################################################
### chunk number 69: nlme eval=FALSE
###################################################
## library("nlme")
## g1 <- subset(Grunfeld, firm == "Westinghouse")
## gls(invest ~ value + capital, data = g1, correlation = corAR1())



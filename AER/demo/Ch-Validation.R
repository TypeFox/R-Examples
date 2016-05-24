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
### chunk number 2: ps-summary
###################################################
data("PublicSchools")
summary(PublicSchools)


###################################################
### chunk number 3: ps-plot eval=FALSE
###################################################
## ps <- na.omit(PublicSchools)
## ps$Income <- ps$Income / 10000
## plot(Expenditure ~ Income, data = ps, ylim = c(230, 830))
## ps_lm <- lm(Expenditure ~ Income, data = ps)
## abline(ps_lm)
## id <- c(2, 24, 48)
## text(ps[id, 2:1], rownames(ps)[id], pos = 1, xpd = TRUE)


###################################################
### chunk number 4: ps-plot1
###################################################
ps <- na.omit(PublicSchools)
ps$Income <- ps$Income / 10000
plot(Expenditure ~ Income, data = ps, ylim = c(230, 830))
ps_lm <- lm(Expenditure ~ Income, data = ps)
abline(ps_lm)
id <- c(2, 24, 48)
text(ps[id, 2:1], rownames(ps)[id], pos = 1, xpd = TRUE)


###################################################
### chunk number 5: ps-lmplot eval=FALSE
###################################################
## plot(ps_lm, which = 1:6)


###################################################
### chunk number 6: ps-lmplot1
###################################################
plot(ps_lm, which = 1:6)


###################################################
### chunk number 7: ps-hatvalues eval=FALSE
###################################################
## ps_hat <- hatvalues(ps_lm)
## plot(ps_hat)
## abline(h = c(1, 3) * mean(ps_hat), col = 2)
## id <- which(ps_hat > 3 * mean(ps_hat))
## text(id, ps_hat[id], rownames(ps)[id], pos = 1, xpd = TRUE)


###################################################
### chunk number 8: ps-hatvalues1
###################################################
ps_hat <- hatvalues(ps_lm)
plot(ps_hat)
abline(h = c(1, 3) * mean(ps_hat), col = 2)
id <- which(ps_hat > 3 * mean(ps_hat))
text(id, ps_hat[id], rownames(ps)[id], pos = 1, xpd = TRUE)


###################################################
### chunk number 9: influence-measures1 eval=FALSE
###################################################
## influence.measures(ps_lm)


###################################################
### chunk number 10: which-hatvalues
###################################################
which(ps_hat > 3 * mean(ps_hat))


###################################################
### chunk number 11: influence-measures2
###################################################
summary(influence.measures(ps_lm))


###################################################
### chunk number 12: ps-noinf eval=FALSE
###################################################
## plot(Expenditure ~ Income, data = ps, ylim = c(230, 830))
## abline(ps_lm)
## id <- which(apply(influence.measures(ps_lm)$is.inf, 1, any))
## text(ps[id, 2:1], rownames(ps)[id], pos = 1, xpd = TRUE)
## ps_noinf <- lm(Expenditure ~ Income, data = ps[-id,])
## abline(ps_noinf, lty = 2)


###################################################
### chunk number 13: ps-noinf1
###################################################
plot(Expenditure ~ Income, data = ps, ylim = c(230, 830))
abline(ps_lm)
id <- which(apply(influence.measures(ps_lm)$is.inf, 1, any))
text(ps[id, 2:1], rownames(ps)[id], pos = 1, xpd = TRUE)
ps_noinf <- lm(Expenditure ~ Income, data = ps[-id,])
abline(ps_noinf, lty = 2)


###################################################
### chunk number 14: journals-age
###################################################
data("Journals")
journals <- Journals[, c("subs", "price")]
journals$citeprice <- Journals$price/Journals$citations
journals$age <- 2000 - Journals$foundingyear


###################################################
### chunk number 15: journals-lm
###################################################
jour_lm <- lm(log(subs) ~ log(citeprice), data = journals)


###################################################
### chunk number 16: bptest1
###################################################
bptest(jour_lm)


###################################################
### chunk number 17: bptest2
###################################################
bptest(jour_lm, ~ log(citeprice) + I(log(citeprice)^2),
  data = journals)


###################################################
### chunk number 18: gqtest
###################################################
gqtest(jour_lm, order.by = ~ citeprice, data = journals)


###################################################
### chunk number 19: resettest
###################################################
resettest(jour_lm)


###################################################
### chunk number 20: raintest
###################################################
raintest(jour_lm, order.by = ~ age, data = journals)


###################################################
### chunk number 21: harvtest
###################################################
harvtest(jour_lm, order.by = ~ age, data = journals)


###################################################
### chunk number 22: 
###################################################
library("dynlm")


###################################################
### chunk number 23: usmacro-dynlm
###################################################
data("USMacroG")
consump1 <- dynlm(consumption ~ dpi + L(dpi),
  data = USMacroG)


###################################################
### chunk number 24: dwtest
###################################################
dwtest(consump1)


###################################################
### chunk number 25: Box-test
###################################################
Box.test(residuals(consump1), type = "Ljung-Box")


###################################################
### chunk number 26: bgtest
###################################################
bgtest(consump1)


###################################################
### chunk number 27: vcov
###################################################
vcov(jour_lm)
vcovHC(jour_lm)


###################################################
### chunk number 28: coeftest
###################################################
coeftest(jour_lm, vcov = vcovHC)


###################################################
### chunk number 29: sandwiches
###################################################
t(sapply(c("const", "HC0", "HC1", "HC2", "HC3", "HC4"),
  function(x) sqrt(diag(vcovHC(jour_lm, type = x)))))


###################################################
### chunk number 30: ps-anova
###################################################
ps_lm <- lm(Expenditure ~ Income, data = ps)
ps_lm2 <- lm(Expenditure ~ Income + I(Income^2), data = ps)
anova(ps_lm, ps_lm2)


###################################################
### chunk number 31: ps-waldtest
###################################################
waldtest(ps_lm, ps_lm2, vcov = vcovHC(ps_lm2, type = "HC4"))


###################################################
### chunk number 32: vcovHAC
###################################################
rbind(SE = sqrt(diag(vcov(consump1))),
  QS = sqrt(diag(kernHAC(consump1))),
  NW = sqrt(diag(NeweyWest(consump1))))


###################################################
### chunk number 33: solow-lm
###################################################
data("OECDGrowth")
solow_lm <- lm(log(gdp85/gdp60) ~ log(gdp60) +
  log(invest) + log(popgrowth + .05), data = OECDGrowth)
summary(solow_lm)


###################################################
### chunk number 34: solow-plot eval=FALSE
###################################################
## plot(solow_lm)


###################################################
### chunk number 35: solow-lts
###################################################
library("MASS")
solow_lts <- lqs(log(gdp85/gdp60) ~ log(gdp60) +
  log(invest) + log(popgrowth + .05), data = OECDGrowth,
  psamp = 13, nsamp = "exact")


###################################################
### chunk number 36: solow-smallresid
###################################################
smallresid <- which(
  abs(residuals(solow_lts)/solow_lts$scale[2]) <= 2.5)


###################################################
### chunk number 37: solow-nohighlev
###################################################
X <- model.matrix(solow_lm)[,-1]
Xcv <- cov.rob(X, nsamp = "exact")
nohighlev <- which(
  sqrt(mahalanobis(X, Xcv$center, Xcv$cov)) <= 2.5)


###################################################
### chunk number 38: solow-goodobs
###################################################
goodobs <- unique(c(smallresid, nohighlev))


###################################################
### chunk number 39: solow-badobs
###################################################
rownames(OECDGrowth)[-goodobs]


###################################################
### chunk number 40: solow-rob
###################################################
solow_rob <- update(solow_lm, subset = goodobs)
summary(solow_rob)


###################################################
### chunk number 41: quantreg
###################################################
library("quantreg")


###################################################
### chunk number 42: cps-lad
###################################################
library("quantreg")
data("CPS1988")
cps_f <- log(wage) ~ experience + I(experience^2) + education
cps_lad <- rq(cps_f, data = CPS1988)
summary(cps_lad)


###################################################
### chunk number 43: cps-rq
###################################################
cps_rq <- rq(cps_f, tau = c(0.25, 0.75), data = CPS1988)
summary(cps_rq)


###################################################
### chunk number 44: cps-rqs
###################################################
cps_rq25 <- rq(cps_f, tau = 0.25, data = CPS1988)
cps_rq75 <- rq(cps_f, tau = 0.75, data = CPS1988)
anova(cps_rq25, cps_rq75)


###################################################
### chunk number 45: cps-rq-anova
###################################################
anova(cps_rq25, cps_rq75, joint = FALSE)


###################################################
### chunk number 46: rqbig
###################################################
cps_rqbig <- rq(cps_f, tau = seq(0.05, 0.95, by = 0.05),
  data = CPS1988)
cps_rqbigs <- summary(cps_rqbig)


###################################################
### chunk number 47: rqbig-plot eval=FALSE
###################################################
## plot(cps_rqbigs)


###################################################
### chunk number 48: rqbig-plot1
###################################################
plot(cps_rqbigs)



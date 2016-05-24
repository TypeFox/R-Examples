### R code from vignette source 'betareg.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(width = 70, prompt = "R> ", continue = "+  ")
library("betareg")


###################################################
### code chunk number 2: beta-distributions
###################################################
par(mfrow = c(1, 2), mar = c(4.1, 4.1, 4.1, 0.1))
dbeta2 <- function(x, mu, phi = 1) dbeta(x, mu * phi, (1 - mu) * phi)
x <- seq(from = 0.01, to = 0.99, length = 200)
xx <- cbind(x, x, x, x, x)

yy <- cbind(
  dbeta2(x, 0.10, 5),
  dbeta2(x, 0.25, 5),
  dbeta2(x, 0.50, 5),
  dbeta2(x, 0.75, 5),
  dbeta2(x, 0.90, 5)
)
matplot(xx, yy, type = "l", xlab = "y", ylab = "Density", main = expression(phi == 5),
  lty = 1, col = "black", ylim = c(0, 15))
text(0.05, 12  , "0.10")
text(0.95, 12  , "0.90")
text(0.22,  2.8, "0.25")
text(0.78,  2.8, "0.75")
text(0.50,  2.3, "0.50")

yy <- cbind(
  dbeta2(x, 0.10, 100),
  dbeta2(x, 0.25, 100),
  dbeta2(x, 0.50, 100),
  dbeta2(x, 0.75, 100),
  dbeta2(x, 0.90, 100)
)
matplot(xx, yy, type = "l", xlab = "y", ylab = "", main = expression(phi == 100),
  lty = 1, col = "black", ylim = c(0, 15))
text(0.10, 14.5, "0.10")
text(0.90, 14.5, "0.90")
text(0.25,  9.8, "0.25")
text(0.75,  9.8, "0.75")
text(0.50,  8.6, "0.50")


###################################################
### code chunk number 3: GasolineYield-betareg
###################################################
data("GasolineYield", package = "betareg")
gy_logit <- betareg(yield ~ batch + temp, data = GasolineYield)


###################################################
### code chunk number 4: GasolineYield-loglog
###################################################
gy_loglog <- betareg(yield ~ batch + temp, data = GasolineYield,
  link = "loglog")


###################################################
### code chunk number 5: GasolineYield-visualization
###################################################
redblue <- hcl(c(0, 260), 90, 40)
plot(yield ~ temp, data = GasolineYield, type = "n",
  ylab = "Proportion of crude oil converted to gasoline",
  xlab = "Temperature at which all gasoline has vaporized",
  main = "Prater's gasoline yield data")
points(yield ~ temp, data = GasolineYield, cex = 1.75, 
  pch = 19, col = rev(gray.colors(10))[as.numeric(batch)])
points(yield ~ temp, data = GasolineYield, cex = 1.75)
legend("topleft", as.character(1:10), title = "Batch",
  col = rev(gray.colors(10)), pch = 19, bty = "n")
legend("topleft", as.character(1:10), title = "Batch", pch = 1, bty = "n")
lines(150:500, predict(gy_logit, 
  newdata = data.frame(temp = 150:500, batch = "6")),
  col = redblue[2], lwd = 2, lty = 2)
lines(150:500, predict(gy_loglog, 
  newdata = data.frame(temp = 150:500, batch = "6")),
  col = redblue[1], lwd = 2)
legend("bottomright", c("log-log", "logit"),
  col = redblue, lty = 1:2, lwd = 2, bty = "n")


###################################################
### code chunk number 6: GasolineYield-betareg1
###################################################
data("GasolineYield", package = "betareg")
gy_logit <- betareg(yield ~ batch + temp, data = GasolineYield)
summary(gy_logit)


###################################################
### code chunk number 7: GasolineYield-plot (eval = FALSE)
###################################################
## set.seed(123)
## plot(gy_logit, which = 1:4, type = "pearson")
## plot(gy_logit, which = 5, type = "deviance", sub.caption = "")
## plot(gy_logit, which = 1, type = "deviance", sub.caption = "")


###################################################
### code chunk number 8: GasolineYield-plot1
###################################################
par(mfrow = c(3, 2))
set.seed(123)
plot(gy_logit, which = 1:4, type = "pearson")
plot(gy_logit, which = 5, type = "deviance", sub.caption = "")
plot(gy_logit, which = 1, type = "deviance", sub.caption = "")


###################################################
### code chunk number 9: GasolineYield-update
###################################################
gy_logit4 <- update(gy_logit, subset = -4)
coef(gy_logit, model = "precision")
coef(gy_logit4, model = "precision")


###################################################
### code chunk number 10: FoodExpenditure-lm
###################################################
data("FoodExpenditure", package = "betareg")
fe_lm <- lm(I(food/income) ~ income + persons, data = FoodExpenditure)


###################################################
### code chunk number 11: FoodExpenditure-betareg
###################################################
fe_beta <- betareg(I(food/income) ~ income + persons,
  data = FoodExpenditure)


###################################################
### code chunk number 12: FoodExpenditure-betareg2
###################################################
fe_beta2 <- betareg(I(food/income) ~ income + persons | persons,
  data = FoodExpenditure)


###################################################
### code chunk number 13: FoodExpenditure-visualization
###################################################
redblueblack <- hcl(c(0, 260, 0), c(90, 90, 0), c(40, 40, 0))
plot(I(food/income) ~ income, data = FoodExpenditure,
  xlab = "Household income", ylab = "Proportion of food expenditures",
  main = "Food expenditures data", type = "n", ylim = c(0.04, 0.57))
points(I(food/income) ~ income, data = FoodExpenditure, cex = 1.75, pch = 19,
  col = rev(gray.colors(7))[persons])
points(I(food/income) ~ income, data = FoodExpenditure, cex = 1.75)
legend("bottomleft", rev(as.character(sort(unique(FoodExpenditure$persons)))),
  title = "Persons", col = gray.colors(7), pch = 19, bty = "n")
legend("bottomleft", rev(as.character(sort(unique(FoodExpenditure$persons)))),
  title = "Persons", pch = 1, bty = "n")
lines(10:100, predict(fe_lm, 
  newdata = data.frame(income = 10:100, persons = mean(FoodExpenditure$persons))),
  col = redblueblack[3], lwd = 2, lty = 2)
lines(10:100, predict(fe_beta, 
  newdata = data.frame(income = 10:100, persons = mean(FoodExpenditure$persons))),
  col = redblueblack[2], lwd = 2, lty = 5)
lines(10:100, predict(fe_beta2, 
  newdata = data.frame(income = 10:100, persons = mean(FoodExpenditure$persons))),
  col = redblueblack[1], lwd = 2)
legend("topright", c("logit, var. disp.", "logit, fix. disp.", "lm"),
  col = redblueblack, lty = c(1, 5, 2), lwd = 2, bty = "n")


###################################################
### code chunk number 14: FoodExpenditure-lm1
###################################################
data("FoodExpenditure", package = "betareg")
fe_lm <- lm(I(food/income) ~ income + persons, data = FoodExpenditure)


###################################################
### code chunk number 15: FoodExpenditure-bptest
###################################################
library("lmtest")
bptest(fe_lm)


###################################################
### code chunk number 16: FoodExpenditure-betareg1
###################################################
fe_beta <- betareg(I(food/income) ~ income + persons,
  data = FoodExpenditure)
summary(fe_beta)


###################################################
### code chunk number 17: GasolineYield-phireg
###################################################
gy_logit2 <- betareg(yield ~ batch + temp | temp, data = GasolineYield)


###################################################
### code chunk number 18: GasolineYield-phireg-coef
###################################################
printCoefmat(summary(gy_logit2)$coefficients$precision)


###################################################
### code chunk number 19: GasolineYield-lrtest
###################################################
lrtest(gy_logit, gy_logit2)


###################################################
### code chunk number 20: FoodExpenditure-betareg2a
###################################################
fe_beta2 <- betareg(I(food/income) ~ income + persons | persons,
  data = FoodExpenditure)


###################################################
### code chunk number 21: FoodExpenditure-comparison
###################################################
lrtest(fe_beta, fe_beta2)
AIC(fe_beta, fe_beta2, k = log(nrow(FoodExpenditure)))


###################################################
### code chunk number 22: GasolineYield-loglog1
###################################################
gy_loglog <- betareg(yield ~ batch + temp, data = GasolineYield,
  link = "loglog")


###################################################
### code chunk number 23: GasolineYield-Rsquared
###################################################
summary(gy_logit)$pseudo.r.squared
summary(gy_loglog)$pseudo.r.squared


###################################################
### code chunk number 24: GasolineYield-AIC
###################################################
AIC(gy_logit, gy_logit2, gy_loglog)


###################################################
### code chunk number 25: GasolineYield-reset
###################################################
lrtest(gy_logit, . ~ . + I(predict(gy_logit, type = "link")^2))
lrtest(gy_loglog, . ~ . + I(predict(gy_loglog, type = "link")^2))


###################################################
### code chunk number 26: GasolineYield-diagnostics
###################################################
plot(abs(residuals(gy_loglog, type = "response")),
  abs(residuals(gy_logit, type = "response")))
abline(0, 1, lty = 2)


###################################################
### code chunk number 27: GasolineYield-diagnostics1 (eval = FALSE)
###################################################
## plot(abs(residuals(gy_loglog, type = "response")),
##   abs(residuals(gy_logit, type = "response")))
## abline(0, 1, lty = 2)


###################################################
### code chunk number 28: GasolineYield-loglog
###################################################
gy_loglog2 <- update(gy_loglog, link.phi = "log")
summary(gy_loglog2)$iterations


###################################################
### code chunk number 29: FoodExpenditure-links
###################################################
sapply(c("logit", "probit", "cloglog", "cauchit", "loglog"),
  function(x) logLik(update(fe_beta2, link = x)))


###################################################
### code chunk number 30: ReadingSkills-eda
###################################################
data("ReadingSkills", package = "betareg")
rs_accuracy <- format(round(with(ReadingSkills, tapply(accuracy, dyslexia, mean)), digits = 3))


###################################################
### code chunk number 31: ReadingSkills-ols
###################################################
data("ReadingSkills", package = "betareg")
rs_ols <- lm(qlogis(accuracy) ~ dyslexia * iq, data = ReadingSkills)


###################################################
### code chunk number 32: ReadingSkills-beta
###################################################
rs_beta <- betareg(accuracy ~ dyslexia * iq | dyslexia + iq,
  data = ReadingSkills, hessian = TRUE)


###################################################
### code chunk number 33: ReadingSkills-visualization
###################################################
cl1 <- hcl(c(260, 0), 90, 40)
cl2 <- hcl(c(260, 0), 10, 95)
plot(accuracy ~ iq, data = ReadingSkills, col = cl2[as.numeric(dyslexia)],
  main = "Reading skills data", xlab = "IQ score", ylab = "Reading accuracy",
  pch = c(19, 17)[as.numeric(dyslexia)], cex = 1.5)
points(accuracy ~ iq, data = ReadingSkills, cex = 1.5,
  pch = (1:2)[as.numeric(dyslexia)], col = cl1[as.numeric(dyslexia)])
nd <- data.frame(dyslexia = "no", iq = -30:30/10)
lines(nd$iq, predict(rs_beta, nd), col = cl1[1], lwd = 2)
lines(nd$iq, plogis(predict(rs_ols, nd)), col = cl1[1], lty = 2, lwd = 2)
nd <- data.frame(dyslexia = "yes", iq = -30:30/10)
lines(nd$iq, predict(rs_beta, nd), col = cl1[2], lwd = 2)
lines(nd$iq, plogis(predict(rs_ols, nd)), col = cl1[2], lty = 2, lwd = 2)
legend("topleft", c("control", "dyslexic", "betareg", "lm"),
  lty = c(NA, NA, 1:2), pch = c(19, 17, NA, NA), lwd = 2,
  col = c(cl2, 1, 1), bty = "n")
legend("topleft", c("control", "dyslexic", "betareg", "lm"),
  lty = c(NA, NA, 1:2), pch = c(1, 2, NA, NA),
  col = c(cl1, NA, NA), bty = "n")


###################################################
### code chunk number 34: ReadingSkills-ols1
###################################################
data("ReadingSkills", package = "betareg")
rs_ols <- lm(qlogis(accuracy) ~ dyslexia * iq, data = ReadingSkills)
coeftest(rs_ols)


###################################################
### code chunk number 35: ReadingSkills-beta1
###################################################
rs_beta <- betareg(accuracy ~ dyslexia * iq | dyslexia + iq,
  data = ReadingSkills, hessian = TRUE)
coeftest(rs_beta)


###################################################
### code chunk number 36: strucchange-data
###################################################
set.seed(123)
y1 <- c(rbeta(150, 0.3 * 4, 0.7 * 4), rbeta(50, 0.5 * 4, 0.5 * 4))
y2 <- c(rbeta(100, 0.3 * 4, 0.7 * 4), rbeta(100, 0.3 * 8, 0.7 * 8))


###################################################
### code chunk number 37: strucchange-gefp
###################################################
library("strucchange")
y1_gefp <- gefp(y1 ~ 1, fit = betareg)
y2_gefp <- gefp(y2 ~ 1, fit = betareg)


###################################################
### code chunk number 38: strucchange-plot1 (eval = FALSE)
###################################################
## plot(y1_gefp, aggregate = FALSE)


###################################################
### code chunk number 39: strucchange-plot2 (eval = FALSE)
###################################################
## plot(y2_gefp, aggregate = FALSE)


###################################################
### code chunk number 40: strucchange-plot (eval = FALSE)
###################################################
## plot(y1_gefp, aggregate = FALSE)
## plot(y2_gefp, aggregate = FALSE)


###################################################
### code chunk number 41: strucchange-plot1a
###################################################
plot(y1_gefp, aggregate = FALSE)


###################################################
### code chunk number 42: strucchange-plot2a
###################################################
plot(y2_gefp, aggregate = FALSE)



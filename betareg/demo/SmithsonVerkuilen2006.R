#######################################
## Example 1: 2 x 2 Factorial Design ##
#######################################

## load packages
library("betareg")
library("lmtest")

## load and attach data
## (NOTE: see ?MockJurors on conflict coding)
data("MockJurors", package = "betareg")
attach(MockJurors)

## Table 1: variable dispersion model
## (NOTE: numerical rather than numerical Hessian are always used for replication,
##  Smithson & Verkuilen erroneously compute one-sided p-values)
mj_vd <- betareg(confidence ~ verdict * conflict | verdict * conflict, hessian = TRUE)
summary(mj_vd)

## inadequate linear models, non-significant F statistic (p. 61)
summary(lm(       confidence  ~ verdict * conflict))
summary(lm(qlogis(confidence) ~ verdict * conflict))

## model selection for beta regression: null model, fixed dispersion model (p. 61)
mj_null <- betareg(confidence ~ 1 | 1)
mj_fd <-   betareg(confidence ~ verdict * conflict | 1)
lrtest(mj_null, mj_fd)
lrtest(mj_null, mj_vd)
## McFadden's pseudo-R-squared
1 - as.vector(logLik(mj_null)/logLik(mj_vd))

## Table 2 (except coefficients and totals): logit(mu), fitted variance, observed variance
tapply(predict(mj_vd, type = "link"), list(conflict, verdict), unique)
tapply(predict(mj_vd, type = "response"), list(conflict, verdict), unique)
tapply(confidence, list(conflict, verdict), mean)

## Table 3 (except coefficients and totals): log(phi), fitted variance, observed variance
tapply(log(predict(mj_vd, type = "precision")), list(conflict, verdict), unique)
tapply(predict(mj_vd, type = "variance"), list(conflict, verdict), unique)
tapply(confidence, list(conflict, verdict), var)

## Figure 3
## (NOTE: slight differences due to data errors in observations 39 and 44
##  in Smithson & Verkuilen graphics code: 0.99 instead of 0.005)
mu <- predict(mj_vd, unique(MockJurors[,1:2]), type = "response")
phi <- predict(mj_vd, unique(MockJurors[,1:2]), type = "precision")
omega <- mu * phi
tau <- phi - mu * phi

par(mfrow = c(2, 2))
hist(confidence[verdict == "two-option" & conflict == "yes"],
  ylim = c(0, 3), freq = FALSE, col = "lightgray", breaks = 0:10/10,
  ylab = "f(x)", xlab = "Confidence", main = "Two-Option Verdict + Conflict")
x <- 0:200/200
lines(x, dbeta(x, omega[1], tau[1]))

hist(confidence[verdict == "two-option" & conflict == "no"],
  ylim = c(0, 3), freq = FALSE, col = "lightgray", breaks = 0:10/10,
  ylab = "f(x)", xlab = "Confidence", main = "Two-Option Verdict + No Conflict")
lines(x, dbeta(x, omega[3], tau[3]))

hist(confidence[verdict == "three-option" & conflict == "yes"],
  ylim = c(0, 3), freq = FALSE, col = "lightgray", breaks = 0:10/10,
  ylab = "f(x)", xlab = "Confidence", main = "Three-Option Verdict + Conflict")
lines(x, dbeta(x, omega[2], tau[2]))

hist(confidence[verdict == "three-option" & conflict == "no"],
  ylim = c(0, 3), freq = FALSE, col = "lightgray", breaks = 0:10/10,
  ylab = "f(x)", xlab = "Confidence", main = "Three-Option Verdict + No Conflict")
lines(x, dbeta(x, omega[4], tau[4]))

par(mfrow = c(1, 1))

detach(MockJurors)


################################################
## Example 2: Stress, Depression, and Anxiety ##
################################################

## load, order, and attach data
data("StressAnxiety", package = "betareg")
StressAnxiety <- StressAnxiety[order(StressAnxiety$stress),]
attach(StressAnxiety)

## Table 4
sa_null   <- betareg(anxiety ~ 1 | 1,           hessian = TRUE)
sa_stress <- betareg(anxiety ~ stress | stress, hessian = TRUE)
summary(sa_null)
summary(sa_stress)
AIC(sa_null, sa_stress)
1 - as.vector(logLik(sa_null)/logLik(sa_stress))

## Figure 4
plot(0, 0, xlim = c(0, 1), ylim = c(0, 8), type = "n",
  main = "", xlab = "x", ylab = "f(x)")
abline(h = 0, col = "lightgray")
lines(density(anxiety, from = 0, to = 1))
lines(density(stress, from = 0, to = 1), lty = 2)
legend("topright", c("Anxiety", "Stress"), lty = 1:2, bty = "n")

## Figure 5
plot(jitter(anxiety) ~ jitter(stress),
  xlab = "Stress", ylab = "Anxiety",
  xlim = c(0, 1), ylim = c(0, 1))
lines(lowess(anxiety ~ stress))
lines(fitted(sa_stress) ~ stress, lty = 2)
lines(fitted(lm(anxiety ~ stress)) ~ stress, lty = 3)
legend("topleft", c("lowess", "betareg", "lm"), lty = 1:3, bty = "n")

## Figure 6
## compute predicted standard errors along stress axis
x <- 0:100/100
x_df <- data.frame(stress = x)
x_se <- sqrt(predict(sa_stress, newdata = x_df, type = "variance"))
## NOTE: Smithson & Verkuilen erroneously compute variances (instead of standard errors)
## and using 1/phi rather than phi (forgetting that they reversed the
## sign of the coefficients in the dispersion/precision model)
x_mu <- predict(sa_stress, newdata = x_df, type = "response")
x_phi <- predict(sa_stress, newdata = x_df, type = "precision")
x_false_se <- x_mu * (1 - x_mu) / (1 + 1/x_phi)

plot(jitter(residuals(sa_stress, type = "response")) ~ jitter(stress),
  xlab = "Stress", ylab = "Response Residuals",
  xlim = c(0, 1), ylim = c(-0.4, 0.4))
abline(h = c(-1, 1) * summary(lm(anxiety ~ stress))$sigma, lty = 3)
lines(x, x_se, lty = 1)
lines(x, -x_se, lty = 1)
lines(x, x_false_se, lty = 2)
lines(x, -x_false_se, lty = 2)
legend("bottomleft", c("OLS", "Beta (erroneous)", "Beta (correct)"),
  lty = 3:1, bty = "n")

detach(StressAnxiety)


#######################################################################################
## Example 3: Sequential Regression With Dyslexia and IQ Predicting Reading Accuracy ##
#######################################################################################

## load and attach data
data("ReadingSkills", package = "betareg")
attach(ReadingSkills)

## Table 5
## OLS regression
## (NOTE: typo in iq coefficient: 0.3954 instead of 0.3594)
rs_ols <- lm(qlogis(accuracy) ~ dyslexia * iq)
summary(rs_ols)
## Beta regression
## (NOTE: Smithson & Verkuilen erroneously compute one-sided p-values)
rs_beta <- betareg(accuracy ~ dyslexia * iq | dyslexia + iq, hessian = TRUE)
summary(rs_beta)

## model selection (p. 66)
rs_beta_1    <- betareg(accuracy ~ 1)
rs_beta_i    <- betareg(accuracy ~ iq | iq)
rs_beta_di   <- betareg(accuracy ~ dyslexia + iq | dyslexia + iq)
rs_beta_dxi  <- betareg(accuracy ~ dyslexia * iq | dyslexia + iq)
rs_beta_dxi2 <- betareg(accuracy ~ dyslexia * iq | dyslexia * iq)
lrtest(rs_beta_1, rs_beta_i, rs_beta_di, rs_beta_dxi, rs_beta_dxi2)

## Table 6
nd <- expand.grid(dyslexia = c("no", "yes"), iq = c(1, -1, 0))
predict(rs_beta, newdata = nd, type = "response")
predict(rs_beta, newdata = nd, type = "variance")

## Figure 1
## compute group-wise beta distribution parameters
nd <- data.frame(dyslexia = c("no", "yes"), iq = tapply(iq, dyslexia, mean))
mu <- predict(rs_beta, newdata = nd, type = "response")
phi <- predict(rs_beta, newdata = nd, type = "precision")
omega <- mu * phi
tau <- phi - mu * phi

par(mfrow = c(1, 2))
hist(accuracy[dyslexia == "yes"],
  breaks = 0:10/10, freq = FALSE, ylim = c(0, 7),
  xlab = "Reading Accuracy", ylab = "f(x)",
  main = "Dyslexics", col = "lightgray")
x <- 0:100/100
lines(x, dbeta(x, omega[2], tau[2]))

hist(accuracy[dyslexia == "no"],
  breaks = 0:10/10, freq = FALSE, ylim = c(0, 7),
  xlab = "Reading Accuracy", ylab = "f(x)",
  main = "Controls", col = "lightgray")
lines(x, dbeta(x, omega[1], tau[1]))

par(mfrow = c(1, 1))

## (NOTE: Standardized residuals seem to be scaled differently.)
plot(scale(fitted(rs_ols)), rstandard(rs_ols),
  xlab = "Standardized Predicted Value", ylab = "Standardized Residual",
  pch = 19, xlim = c(-1.5, 1.5), ylim = c(-3, 2))

detach(ReadingSkills)


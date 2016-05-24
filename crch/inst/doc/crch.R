### R code from vignette source 'crch.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(width = 70, prompt = "R> ", continue = "+  ")
library("crch")


###################################################
### code chunk number 2: crch.Rnw:342-344
###################################################
library("crch")
data("RainIbk", package = "crch")


###################################################
### code chunk number 3: crch.Rnw:360-364
###################################################
RainIbk <- sqrt(RainIbk)
RainIbk$ensmean <- apply(RainIbk[,grep('^rainfc',names(RainIbk))], 1, mean)
RainIbk$enssd <- apply(RainIbk[,grep('^rainfc',names(RainIbk))], 1, sd)
RainIbk <- subset(RainIbk, enssd > 0)


###################################################
### code chunk number 4: crch.Rnw:369-371
###################################################
plot(rain ~ ensmean, data = RainIbk, pch = 19, col = gray(0, alpha = 0.2))
abline(0,1, col = "red")


###################################################
### code chunk number 5: crch.Rnw:382-385
###################################################
CRCH <- crch(rain ~ ensmean | log(enssd), data = RainIbk, left = 0, 
  dist = "logistic")
summary(CRCH)


###################################################
### code chunk number 6: crch.Rnw:392-393
###################################################
abline(coef(CRCH)[1:2], col = "blue")


###################################################
### code chunk number 7: crch.Rnw:398-401
###################################################
plot(rain~ensmean, data = RainIbk, pch = 19, col = gray(0, alpha = 0.2))
abline(0,1, col = "red", lwd = 2)
abline(coef(CRCH)[1:2], col = "blue", lwd = 2)


###################################################
### code chunk number 8: crch.Rnw:409-411
###################################################
CR <- crch(rain ~ ensmean, data = RainIbk, left = 0, dist = "logistic")
cbind(AIC(CR, CRCH), BIC = BIC(CR, CRCH)[,2])


###################################################
### code chunk number 9: crch.Rnw:418-423
###################################################
CRCHgau <- crch(rain ~ ensmean | log(enssd), data = RainIbk, left = 0, 
  dist = "gaussian")
CRCHstud <- crch(rain ~ ensmean | log(enssd), data = RainIbk, left = 0, 
  dist = "student")
AIC(CRCH, CRCHgau, CRCHstud)


###################################################
### code chunk number 10: crch.Rnw:433-443
###################################################
a <- seq(-5, 5, 0.01)
df <- exp(tail(coef(CRCHstud), 1))
#par(mfrow = c(1,2))
plot(a, dt(a, df), type = "l", xlab = "", ylab = "density", main = "Probability density function", lwd = 2)
lines(a, dlogis(a*dt(0, df)/dlogis(0))*dt(0, df)/dlogis(0), col = 2, lwd = 2, lty = 2)
lines(a, dnorm(a*dt(0, df)/dnorm(0))*dt(0, df)/dnorm(0), col = 4, lwd = 1, lty = 1)

#plot(a, pt(a, df), type = "l", xlab = "", ylab = "density", main = "Cumulative distribution function", lwd = 2)
#lines(a, plogis(a*dt(0, df)/dlogis(0)), col = 2, lwd = 2, lty = 2)
legend("topright", lwd = c(2,2,1), lty = c(1,2,1), col = c(1,2,4), c("student-t", "scaled logistic", "scaled normal"), bty = "n")


###################################################
### code chunk number 11: crch.Rnw:457-462
###################################################
library("glmx")
BIN <- hetglm(I(rain > 0) ~ ensmean | log(enssd), data = RainIbk,
         family = binomial(link = "logit"))
TRCH <- crch(rain~ensmean | log(enssd), data = RainIbk, subset = rain > 0, 
         left = 0, dist = "logistic", truncated = TRUE)


###################################################
### code chunk number 12: crch.Rnw:469-474
###################################################
cbind("CRCH" = c(coef(CRCH, "location")/exp(coef(CRCH, "scale"))[1], 
        coef(CRCH, "scale")[2]), 
      "BIN" = coef(BIN), 
      "TRCH" = c(coef(TRCH, "location")/exp(coef(TRCH, "scale"))[1], 
        coef(TRCH, "scale")[2]))


###################################################
### code chunk number 13: crch.Rnw:481-486
###################################################
loglik <- c("Censored" = logLik(CRCH), "Two-Part" = logLik(BIN) + logLik(TRCH))
df <- c(4, 7)
aic <- -2 * loglik + 2 * df
bic <- -2 * loglik + log(nrow(RainIbk)) * df
cbind(df, AIC = aic, BIC = bic)


###################################################
### code chunk number 14: crch.Rnw:494-496
###################################################
newdata <- data.frame(ensmean = 1.8, enssd = 0.9)
predict(CRCH, newdata, type = "quantile", at = 0.5)^2


###################################################
### code chunk number 15: crch.Rnw:507-509
###################################################
p <- predict(BIN, newdata)
predict(TRCH, newdata, type = "quantile", at = (p - 0.5)/p)^2


###################################################
### code chunk number 16: crch.Rnw:516-523
###################################################
mu <- predict(CRCH, newdata, type = "location")
sigma <- predict(CRCH, newdata, type = "scale")
pclogis(sqrt(5), mu, sigma, lower.tail = FALSE, left = 0)

mu <- predict(TRCH, newdata, type = "location")
sigma <- predict(TRCH, newdata, type = "scale")
p * ptlogis(sqrt(5), mu, sigma, lower.tail = FALSE, left = 0)



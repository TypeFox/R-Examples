### R code from vignette source 'mcemGLM-vignette.Rnw'

###################################################
### code chunk number 1: mcemGLM-vignette.Rnw:186-187
###################################################
set.seed(23786)


###################################################
### code chunk number 2: mcemGLM-vignette.Rnw:189-192
###################################################
require(mcemGLM)
data("salamander")
summary(salamander)


###################################################
### code chunk number 3: mcemGLM-vignette.Rnw:204-209
###################################################
fitBernoulli <- mcemGLMM(fixed  = Mate ~ 0+Cross,
                         random = list(~ 0+Female, ~ 0+Male), 
                         data   = salamander, 
                         family = "bernoulli", 
                         vcDist = "normal")


###################################################
### code chunk number 4: mcemGLM-vignette.Rnw:223-224
###################################################
summary(fitBernoulli)


###################################################
### code chunk number 5: mcemGLM-vignette.Rnw:229-239
###################################################
ctr0 <- matrix(c(1, -1, 0, 0, 
                 1, 0, -1, 0, 
                 1, 0, 0, -1,
                 0, 1, -1, 0,
                 0, 1, 0, -1,
                 0, 0, 1, -1), 6, 4, byrow = TRUE)

rownames(ctr0) <- c("RR - RW", "RR - WR", "RR - WW", 
                    "RW - WR", "RW - WW", "WR - WW")



###################################################
### code chunk number 6: mcemGLM-vignette.Rnw:243-244
###################################################
contrasts.mcemGLMM(fitBernoulli, ctr0)


###################################################
### code chunk number 7: mcemGLM-vignette.Rnw:252-255
###################################################
par(mfrow=c(1, 2))
plot(residuals(fitBernoulli, type = "deviance"), xlab = "Observation", ylab = "di", main = "Deviance residuals")
plot(residuals(fitBernoulli, type = "deviance")~salamander$Cross, xlab = "Cross", ylab = "di", main = "Deviance residuals")


###################################################
### code chunk number 8: mcemGLM-vignette.Rnw:261-265
###################################################
predict(fitBernoulli, newdata=list(Cross = c("RR", "RW", "WR", "WW")), 
        type = "link", se.fit = TRUE)
predict(fitBernoulli, newdata=list(Cross = c("RR", "RW", "WR", "WW")), 
        type = "response", se.fit = TRUE)


###################################################
### code chunk number 9: mcemGLM-vignette.Rnw:273-276
###################################################
par(mfrow=c(1, 2))
qqnorm(ranef.mcemGLMM(fitBernoulli)[1:60], main = "Normal Q-Q Plot for \n Female Random Effects")
qqnorm(ranef.mcemGLMM(fitBernoulli)[61:120], main = "Normal Q-Q Plot for \n Male Random Effects")


###################################################
### code chunk number 10: mcemGLM-vignette.Rnw:285-287
###################################################
ts.plot(fitBernoulli$mcemEST, col=1:6, xlab = "EM Iteration", ylab = "Estimate")
legend("topright", lty=1, col = 1:6, c("RR", "RW", "WR", "WW", expression(paste(sigma,"2", "_f")), expression(paste(sigma,"2", "_m"))))


###################################################
### code chunk number 11: mcemGLM-vignette.Rnw:296-297
###################################################
ts.plot(fitBernoulli$QfunVal, xlab = "EM iteration", ylab="Q function")


###################################################
### code chunk number 12: mcemGLM-vignette.Rnw:307-312
###################################################
par(mfrow = c(2,2))
ts.plot(fitBernoulli$randeff[, 1], xlab = "MC iteration", ylab = "U_F1")
acf(fitBernoulli$randeff[, 1], lag.max = 100, main = "U_F1")
ts.plot(fitBernoulli$QfunMCMC, xlab = "MC iteration", ylab = "Q function")
acf(fitBernoulli$QfunMCMC, lag.max = 100)


###################################################
### code chunk number 13: mcemGLM-vignette.Rnw:333-339
###################################################
fitNegbin <- mcemGLMM(fixed  = count ~ base * group + age + visit, 
                 random = list( ~ 0 + id),
                 data   = epilepsy, 
                 family = "negbinom", 
                 vcDist = "t", 
                 df     = 10)


###################################################
### code chunk number 14: mcemGLM-vignette.Rnw:342-343
###################################################
summary(fitNegbin)



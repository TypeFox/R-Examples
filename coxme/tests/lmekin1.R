library(coxme)
require(nlme)
aeq <- function(x, y, ...) all.equal(as.vector(x), as.vector(y), ...)
tdata <- eortc
tdata$center2 <- factor(tdata$center)
tdata$trt2 <- factor(tdata$trt)

fit1 <- lme(y ~ trt, random= ~ 1|center2, data=tdata,
            method="ML")
fit2 <- lmekin(y ~ trt + (1|center), tdata)

#I get the same loglik but slightly different coefficients.
# (The loglik is flat on top, and we are using different starting
#  estimates.)  So we need to use a toloerance for the test.
vcomp <- function(m1, m2, ...) {
    temp1 <- as.numeric(VarCorr(m1)[,2])  #lme gives a matrix of text
#    temp2 <- m2$sigma^2 * c(unlist(VarCorr(m2)), 1) #lmekin the variance struct
    temp2 <- c(unlist(VarCorr(m2)), m2$sigma^2) #lmekin the variance struct
    aeq(temp1, sqrt(temp2), ...)
}

aeq(fit1$logLik, fit2$loglik)
aeq(fit1$sigma, fit2$sigma, tol=1e-4)
aeq(fixef(fit1), fixef(fit2), tol=1e-4)
vcomp(fit1, fit2, tol=1e-3)

# A second fit, using the matrix form
ucen <- sort(unique(tdata$center))
kmat <- diag(length(ucen))
dimnames(kmat) <- list(ucen, ucen)
fit3 <- lmekin(y ~ trt + (1|center), tdata, varlist=kmat)
aeq(fit3$coefficients, fit2$coefficients)
aeq(unlist(VarCorr(fit3)), unlist(VarCorr(fit2)))

# Force the same coefs. 
temp <- as.numeric(unclass(fit1$modelStruct$reStruct)[[1]]) 
fit3 <- lmekin(y~trt + (1|center), tdata,
               vfixed= exp(2*temp))
aeq(fit1$logLik, fit3$loglik)
aeq(fit1$sigma, fit3$sigma)
aeq(fixef(fit1), fixef(fit3))
aeq(ranef(fit1), ranef(fit3), check.attributes=FALSE)
vcomp(fit1, fit3, tol=1e-7)
aeq(residuals(fit1), residuals(fit3))

# Now a model with random slopes and intercepts
fit5 <- lme(y ~ trt, random= ~ trt|center, data=tdata,
            method="REML")
fit6 <- lmekin(y~ trt + (1+trt | center), data=tdata,
               method="REML")
aeq(fit5$logLik, fit6$loglik)
aeq(fit5$sigma, fit6$sigma, tol=1e-4)
aeq(fixef(fit5), fixef(fit6), tol=1e-4)


#
# Use an lme data set
#
mystool <- as.data.frame(ergoStool) #get rid of contrast attributes

efit1 <-  lme(effort ~ Type, data=mystool, random= ~1|Subject,
            method="ML")
efit2 <- lmekin(effort ~ Type + (1|Subject), mystool)
aeq(efit1$logLik, efit2$loglik)
aeq(fixef(efit1), fixef(efit2))
vcomp(efit1, efit2, tol=1e-5)
aeq(efit1$sigma, efit2$sigma, tol=1e-5)

efit3 <-lme(effort ~ Type, data=mystool, random= ~1|Subject,
            method="REML")
efit4 <- lmekin(effort ~ Type + (1|Subject), mystool, method="REML")

aeq(efit3$logLik, efit4$loglik)
aeq(fixef(efit3), fixef(efit4))
vcomp(efit3, efit4, tol=1e-3)
aeq(efit3$sigma, efit4$sigma, tol=1e-4)

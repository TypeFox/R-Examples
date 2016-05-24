require("DoseFinding")
data(IBScovars)
lmfit <- lm(resp~factor(dose)+gender, data=IBScovars)
cf <- coef(lmfit)[-c(6)]
vcv <- vcov(lmfit)[-c(6), -c(6)]
lmfit2 <- lm(resp~as.factor(dose)-1+gender, data=IBScovars)
cf2 <- coef(lmfit2)[-c(6)]
vcv2 <- vcov(lmfit2)[-c(6), -c(6)]
dose <- c(0:4)

## test fitting all available models
fitMod(dose[-1], cf[-1], S=vcv[-1,-1], model="linear", placAdj=TRUE,type="general")
fitMod(dose, cf2, S=vcv2, model="linear", placAdj=FALSE,type="general")
fitMod(dose[-1], cf[-1], S=vcv[-1,-1], model="quadratic", placAdj=TRUE,type="general")
fitMod(dose, cf2, S=vcv2, model="quadratic", placAdj=FALSE,type="general")
fitMod(dose, cf2, S=vcv2, model="linlog", placAdj=FALSE,
       addArgs=list(off=0.01*max(dose)),type="general")
fitMod(dose[-1], cf[-1], S=vcv[-1,-1], model="emax", placAdj=TRUE,
       bnds=defBnds(max(dose))$emax,type="general")
fitMod(dose, cf2, S=vcv2, model="emax", placAdj=FALSE,
       bnds=defBnds(max(dose))$emax,type="general")
fitMod(dose[-1], cf[-1], S=vcv[-1,-1], model="sigEmax", placAdj=TRUE,
       bnds=defBnds(max(dose))$sigEmax,type="general")
fitMod(dose, cf2, S=vcv2, model="sigEmax", placAdj=FALSE,
       bnds=defBnds(max(dose))$sigEmax,type="general")
fitMod(dose[-1], cf[-1], S=vcv[-1,-1], model="exponential",
       placAdj=TRUE, bnds=defBnds(max(dose))$exponential,type="general")
fitMod(dose, cf2, S=vcv2, model="exponential", placAdj=FALSE,
       bnds=defBnds(max(dose))$exponential,type="general")
fitMod(dose, cf2, S=vcv2, model="logistic", placAdj=FALSE,
       bnds=defBnds(max(dose))$logistic,type="general")
fitMod(dose[-1], cf[-1], S=vcv[-1,-1], model="betaMod", placAdj=TRUE,
       bnds=defBnds(max(dose))$betaMod, addArgs=list(scal=1.2*4),type="general")
fitMod(dose, cf2, S=vcv2, model="betaMod", placAdj=FALSE,
       bnds=defBnds(max(dose))$betaMod, addArgs=list(scal=1.2*4),type="general")
fitMod(dose[-1], cf[-1], S=vcv[-1,-1], model="linInt", placAdj=TRUE, type="general")
fitMod(dose, cf2, S=vcv2, model="linInt", placAdj=FALSE, type="general")
## test using starting value (instead of grid search)
fitMod(dose[-1], cf[-1], S=vcv[-1,-1], model="emax",
       placAdj=TRUE, bnds=defBnds(max(dose))$emax, start = 0.5,type="general")
fitMod(dose, cf2, S=vcv2, model="emax", placAdj=FALSE,
       bnds=defBnds(max(dose))$emax, start = 0.2,type="general")
fitMod(dose[-1], cf[-1], S=vcv[-1,-1], model="betaMod",
       placAdj=TRUE, bnds=defBnds(max(dose))$betaMod,
       addArgs=list(scal=1.2*4),type="general")
fitMod(dose, cf2, S=vcv2, model="betaMod", placAdj=FALSE,
       bnds=defBnds(max(dose))$betaMod, start = c(1, 1),
       addArgs=list(scal=1.2*4),type="general")

## test predict, vcov, coef, intervals, plot, summary
ggI <- fitMod(dose, cf2, S=vcv2, model="betaMod",
              placAdj=FALSE, bnds=defBnds(max(dose))$betaMod,
              addArgs=list(scal=1.2*4),type="general")
ggNI <- fitMod(dose[-1], cf[-1], S=vcv[-1,-1], model="betaMod",
               placAdj=TRUE, bnds=defBnds(max(dose))$betaMod,
               addArgs=list(scal=1.2*4),type="general")
predict(ggI, se.fit=TRUE, predType = "e")
predict(ggNI, se.fit=TRUE, predType = "e")
vcov(ggI)
vcov(ggNI)
plot(ggI, CI=T, plotData = "meansCI")
plot(ggNI, CI=T, plotData = "meansCI")

ggI <- fitMod(dose, cf2, S=vcv2, model="linInt",
              placAdj=FALSE,type="general")
ggNI <- fitMod(dose[-1], cf[-1], S=vcv[-1,-1], model="linInt",
               placAdj=TRUE,type="general")
predict(ggI, se.fit=TRUE, predType = "full-model")
predict(ggI, se.fit=TRUE, predType = "effect-curve")
predict(ggNI, se.fit=TRUE, predType = "full-model")
vcov(ggI)
vcov(ggNI)
plot(ggI, CI=T, plotData = "meansCI")
plot(ggNI, CI=T, plotData = "meansCI")

## even more tests for the linInt model
data(IBScovars)
## without covariates
fit <- fitMod(dose, resp, data=IBScovars, model="linInt")
plot(fit, CI=TRUE, plotData="meansCI")
fit <- fitMod(dose, resp, data=IBScovars, model="linInt",
                  addCovars=~gender)
plot(fit, CI=TRUE, plotData="meansCI")
vcov(fit)
fit <- lm(resp~as.factor(dose)-1, data=IBScovars)
cf <- coef(fit)
vc <- vcov(fit)
doseVec <- 0:4
fit <- fitMod(doseVec, cf, model="linInt", S=vc, type = "general")
plot(fit, CI=TRUE, plotData="meansCI")
vcov(fit)
fit <- lm(resp~as.factor(dose)+gender, data=IBScovars)
cf <- coef(fit)[2:5]
vc <- vcov(fit)[2:5,2:5]
doseVec <- 1:4
fit <- fitMod(doseVec, cf, model="linInt", S=vc, type = "general", placAdj=TRUE)
vcov(fit)
plot(fit, CI=TRUE, plotData="meansCI")
predict(fit, predType = "effect-curve", se.fit=TRUE)

library(partsm)

###########################
###### acf.ext1
###########################

data("gergnp")
lgergnp <- log(gergnp, base=exp(1))

out <- acf.ext1(wts=lgergnp, transf.type="orig",
                type="correlation", lag.max=12, showcat=TRUE, plot=FALSE)

out <- acf.ext1(wts=lgergnp, transf.type="perdiffsd", perdiff.coeff = c(1.004, 0.981, 1.047, 0.969),
                type="correlation", lag.max=12, showcat=TRUE, plot=FALSE)

###########################
###### fit.ar.par
###########################

## Models for the the logarithms of the Real GNP in Germany.
data("gergnp")
lgergnp <- log(gergnp, base=exp(1))

## Fit an AR(4) model with intercept and seasonal dummies.
detcomp <- list(regular=c(1,0,c(1,2,3)), seasonal=c(0,0), regvar=0)
out.ar <- fit.ar.par(wts=lgergnp, type="AR", detcomp=detcomp, p=4)

## Fit a PAR(2) model with seasonal intercepts.
detcomp <- list(regular=c(0,0,0), seasonal=c(1,0), regvar=0)
out.par <- fit.ar.par(wts=lgergnp, type="PAR", detcomp=detcomp, p=2)

###########################
###### fit.piar
###########################

## Fit a PIAR(2) model for the logarithms of the Real GNP in Germany.
data("gergnp")
lgergnp <- log(gergnp, base=exp(1))
detcomp <- list(regular=c(0,0,0), seasonal=c(1,0), regvar=0)
out <- fit.piar(wts=lgergnp, detcomp=detcomp, p=2, initvalues=NULL)

###########################
###### Fnextp.test
###########################

## Test the significance of a second order lag in a PAR model for the Real GNP in Germany.
## Including seasonal intercepts.
data("gergnp")
lgergnp <- log(gergnp, base=exp(1))
detcomp <- list(regular=c(0,0,0), seasonal=c(1,0), regvar=0)
out <- Fnextp.test(wts=lgergnp, detcomp=detcomp, p=1, type="PAR")

###########################
###### Fpar.test
###########################

## Test for periodicity in a second order PAR model for
## the logarithms of the Real GNP in Germany time series.
data("gergnp")
lgergnp <- log(gergnp, base=exp(1))
detcomp <- list(regular=c(0,0,0), seasonal=c(1,0), regvar=0)
out <- Fpar.test(wts=lgergnp, detcomp=detcomp, p=2)


###########################
###### Fpari.piar.test
###########################

## Test for the unit root 1 in a PAR(2) with seasonal intercepts for
## the logarithms of the Real GNP in Germany.
data("gergnp")
lgergnp <- log(gergnp, base=exp(1))
detcomp <- list(regular=c(0,0,0), seasonal=c(1,0), regvar=0)
out <- Fpari.piar.test(wts=lgergnp, detcomp=detcomp, p=2, type="PARI1")

###########################
###### Fsh.test
###########################

## Fsh test for the residuals of the first differences
## of the logarithms of the Real GNP in Germany
## on an AR(4) model with seasonal intercepts.
data("gergnp")
lgergnp <- log(gergnp, base=exp(1))
wts <- ts(c(NA, diff(gergnp, lag=1)), frequency=4, start=start(lgergnp))

detcomp=list(regular=c(0,0,0), seasonal=c(1,0), regvar=0)
ar4 <- fit.ar.par(wts=lgergnp, type="AR", p=4, detcomp=detcomp)
out <- Fsh.test(res=residuals(ar4@lm.ar), s=frequency(wts))

###########################
###### LRurpar.test
###########################

## Test for a single unit root in a PAR(2) model with seasonal intercepts for the
## logarithms of the Real GNP in Germany.
data("gergnp")
lgergnp <- log(gergnp, base=exp(1))
detcomp <- list(regular=c(0,0,0), seasonal=c(1,0), regvar=0)
out <- LRurpar.test(wts=lgergnp, detcomp=detcomp, p=2)

###########################
###### PAR.MVrepr
###########################

## Models for the the logarithms of the Real GNP in Germany.
data("gergnp")
lgergnp <- log(gergnp, base=exp(1))

## Fit an PAR model
detcomp <- list(regular=c(0,0,0), seasonal=c(1,0), regvar=0)
out.par <- fit.ar.par(wts=lgergnp, type="PAR", detcomp=detcomp, p=2)

## Show the matrix representation: 
out.MV <- PAR.MVrepr(out.par)
out.MV

###########################
###### PAR.MVrepr-methods
###########################

## Load data and select the deterministic components.
data("gergnp")
lgergnp <- log(gergnp, base=exp(1))
detcomp <- list(regular=c(0,0,0), seasonal=c(1,0), regvar=0)

## Multivariate representation of a PAR(2) model with sesonal intercepts.
out.par <- fit.ar.par(wts=lgergnp, type="PAR", detcomp=detcomp, p=2)
PAR.MVrepr(out.par)

## Multivariate representation of a PIAR(2) model with sesonal intercepts.
out.piar <- fit.piar(wts=lgergnp, detcomp=detcomp, p=2)
PAR.MVrepr(out.piar)

###########################
###### plotpdiff
###########################

## Load data and select the deterministic components.
data("gergnp")
lgergnp <- log(gergnp, base=exp(1))
detcomp <- list(regular=c(0,0,0), seasonal=c(1,0), regvar=0)

## Fit a PIAR(2) model with seasonal intercepts.
out.piar <- fit.piar(wts=lgergnp, detcomp=detcomp, p=2)
plotpdiff(out.piar)

###########################
###### plotpredpiar
###########################

## Load data and select the deterministic components.
data("gergnp")
lgergnp <- log(gergnp, base=exp(1))

## Fit a PIAR(2) model with seasonal intercepts.
out.pred <- predictpiar(wts=lgergnp, p=2, hpred=24)
plotpredpiar(out.pred)

###########################
###### predictpiar
###########################

## 24 step-ahead forecasts in a PIAR(2) model for the
## logarithms of the Real GNP in Germany.
data("gergnp")
lgergnp <- log(gergnp, base=exp(1))
pred.out <- predictpiar(wts=lgergnp, p=2, hpred=24)

###########################
###### show-methods
###########################

## Load data and select the deterministic components.
data("gergnp")
lgergnp <- log(gergnp, base=exp(1))
detcomp <- list(regular=c(0,0,0), seasonal=c(1,0), regvar=0)

## Fit an AR(4) model with intercept and seasonal dummies.
dcar <- list(regular=c(1,0,c(1,2,3)), seasonal=c(0,0), regvar=0)
out.ar <- fit.ar.par(wts=lgergnp, type="AR", detcomp=dcar, p=4)
show(out.ar)

## Fit a PAR(2) model with seasonal intercepts.
out.par <- fit.ar.par(wts=lgergnp, type="PAR", detcomp=detcomp, p=2)
show(out.par)

## Fnextp.test
Fnextp.out <- Fnextp.test(wts=lgergnp, detcomp=detcomp, p=1, type="PAR")
show(Fnextp.out)

## Fpar.test
Fpar.out <- Fpar.test(wts=lgergnp, detcomp=detcomp, p=2)
show(Fpar.out)

## Fsh.test
ar4 <- fit.ar.par(wts=lgergnp, type="AR", p=4, detcomp=detcomp)
Fsh.out <- Fsh.test(res=residuals(ar4@lm.ar), s=frequency(lgergnp))
show(Fsh.out)

## Fit a PIAR(2) model with seasonal intercepts.
out.piar <- fit.piar(wts=lgergnp, detcomp=detcomp, p=2)
show(out.piar)

## Fpari.piar.test
Fpari1.out <- Fpari.piar.test(wts=lgergnp, detcomp=detcomp, p=2, type="PARI1")
show(Fpari1.out)

## Fit a PIAR(2) model with seasonal intercepts.
out.piar <- fit.piar(wts=lgergnp, detcomp=detcomp, p=2)
show(out.piar)

## Test for a single unit root in a PAR(2) model with seasonal intercepts.
out.LR <- LRurpar.test(wts=lgergnp, detcomp=detcomp, p=2)
show(out.LR)

## 24 step-ahead forecasts in a PIAR(2) model.
pred.out <- predictpiar(wts=lgergnp, p=2, hpred=24)
options(digits=4)
show(pred.out)
options(digits=7)

###########################
###### summary-methods
###########################

## Load data and select the deterministic components.
data("gergnp")
lgergnp <- log(gergnp, base=exp(1))
detcomp <- list(regular=c(0,0,0), seasonal=c(1,0), regvar=0)

## Fit an AR(4) model with intercept and seasonal dummies.
dcar <- list(regular=c(1,0,c(1,2,3)), seasonal=c(0,0), regvar=0)
out.ar <- fit.ar.par(wts=lgergnp, type="AR", detcomp=dcar, p=4)
summary(out.ar)

## Fit a PAR(2) model with seasonal intercepts.
detcomp <- list(regular=c(0,0,0), seasonal=c(1,0), regvar=0)
out.par <- fit.ar.par(wts=lgergnp, type="PAR", detcomp=detcomp, p=2)
summary(out.par)

## Fnextp.test
Fnextp.out <- Fnextp.test(wts=lgergnp, detcomp=detcomp, p=1, type="PAR")
summary(Fnextp.out)

## Fpar.test
Fpar.out <- Fpar.test(wts=lgergnp, detcomp=detcomp, p=2)
summary(Fpar.out)

## Fsh.test
ar4 <- fit.ar.par(wts=lgergnp, type="AR", p=4, detcomp=detcomp)
Fsh.out <- Fsh.test(res=residuals(ar4@lm.ar), s=frequency(lgergnp))
summary(Fsh.out)

## Fit a PIAR(2) model.
out.piar <- fit.piar(wts=lgergnp, detcomp=detcomp, p=2)
summary(out.piar, digits=2)

## Fpari.piar.test
Fpari1.out <- Fpari.piar.test(wts=lgergnp, detcomp=detcomp, p=2, type="PARI1")
options(digits=3)
summary(Fpari1.out)
options(digits=3)

## Fit a PIAR(2) model with seasonal intercepts.
out.piar <- fit.piar(wts=lgergnp, detcomp=detcomp, p=2)
summary(out.piar, digits=2)

## Test for a single unit root in a PAR(2) model with seasonal intercepts.
out.LR <- LRurpar.test(wts=lgergnp, detcomp=detcomp, p=2)
options(digits=1)
summary(out.LR)
options(digits=7)

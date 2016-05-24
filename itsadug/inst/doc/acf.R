## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(itsadug)

## ------------------------------------------------------------------------
library(itsadug)
library(mgcv)
data(simdat)
# add missing values to simdat:
simdat[sample(nrow(simdat), 15),]$Y <- NA

## ------------------------------------------------------------------------
simdat <- start_event(simdat, column="Time", event=c("Subject", "Trial"), label.event="Event")
head(simdat)

## ------------------------------------------------------------------------
library(mgcv)
# example model:
m1 <- bam(Y ~ te(Time, Trial)+s(Subject, bs='re'), data=simdat)

## ---- fig.width=12, fig.height=4-----------------------------------------
par(mfrow=c(1,3), cex=1.1)

# default ACF function:
acf(resid(m1), main="acf(resid(m1))")
# resid_gam:
acf(resid_gam(m1), main="acf(resid_gam(m1))")
# acf_resid:
acf_resid(m1, main="acf_resid(m1)")

## ---- fig.width=4, fig.height=4------------------------------------------
# we also ask to plot the ACF by specifying plot (FALSE by default):
r1 <- start_value_rho(m1, plot=TRUE)

## ------------------------------------------------------------------------
acf(resid(m1), plot=FALSE)$acf[2]

## ------------------------------------------------------------------------
# example model:
m1AR1 <- bam(Y ~ te(Time, Trial)+s(Subject, bs='re'), data=simdat, rho=r1, AR.start=simdat$start.event)

## ---- fig.width=8, fig.height=4------------------------------------------
par(mfrow=c(1,2), cex=1.1)
acf(resid(m1))
acf(resid(m1AR1))

## ---- fig.width=8, fig.height=4------------------------------------------
par(mfrow=c(1,2), cex=1.1)
acf_resid(m1)
acf_resid(m1AR1)

## ---- fig.width=12, fig.height=4-----------------------------------------
par(mfrow=c(1,3), cex=1.1)
acf_resid(m1AR1, split_pred = c("Subject", "Trial"))
# alternatively, if the predictors are not found in the model we can use the data:
acf_resid(m1AR1, split_pred=list(Subject=simdat$Subject, Trial=simdat$Trial))
# ... or the AR.start information, if provided to the model:
acf_resid(m1AR1, split_pred="AR.start")

## ---- fig.width=9, fig.height=6------------------------------------------
par(cex=1.1)
acf_resid(m1AR1, split_pred = c("Subject", "Trial"), n=6)

## ---- fig.width=8, fig.height=4------------------------------------------
par(mfrow=c(1,2), cex=1.1)
# normal residuals:
normal.res <- resid(m1AR1)
acf(normal.res)
# corrected residuals:
corrected.res <- resid_gam(m1AR1)
acf(corrected.res)

## ----error=TRUE----------------------------------------------------------
# This will elicit an error:
simdat$res.m1 <- resid(m1)
# solution:
simdat$res.m1 <- NA
simdat[!is.na(simdat$Y),]$res.m1 <- resid(m1)

# This will generate an error:
simdat$res.m1AR1 <- resid_gam(m1AR1)
# ... and this too!
simdat$res.m1AR1 <- NA
simdat[!is.na(simdat$Y),]$res.m1AR1 <- resid_gam(m1AR1)
# solution:
simdat$res.m1AR1 <- NA
simdat[!is.na(simdat$Y),]$res.m1AR1 <- resid_gam(m1AR1, incl_na=TRUE)

## ---- fig.width=12, fig.height=4-----------------------------------------
par(mfrow=c(1,3), cex=1.1)
acf(resid_gam(m1AR1))
acf_plot(resid_gam(m1AR1))
acf_resid(m1AR1)

## ---- fig.width=12, fig.height=4-----------------------------------------
par(mfrow=c(1,3), cex=1.1)
# when using acf_plot one need to remove missing values manually:
acf_plot(resid_gam(m1AR1, incl_na = TRUE), 
         split_by=list(Subject=simdat[!is.na(simdat$Y),]$Subject,
                       Trial=simdat[!is.na(simdat$Y),]$Trial))
# ... acf_resid takes care of that automatically:
acf_resid(m1AR1, split_pred=c("Subject", "Trial"))
# ... also when using a list to identify time series:
acf_resid(m1AR1, split_pred=list(Subject=simdat$Subject,
                       Trial=simdat$Trial))

## ---- fig.width=4, fig.height=4------------------------------------------
tmp <- simdat[!is.na(simdat$Y),]
# default function is mean:
acf.y <- acf_plot(tmp$res.m1, 
         split_by=list(Subject=tmp$Subject, Trial=tmp$Trial),
         main="ACF with standard deviation")
points(as.numeric(names(acf.y)),acf.y, pch=16, cex=.5)
# alternatively, we could ask for SE:
acf.se <- acf_plot(tmp$res.m1AR1, 
         split_by=list(Subject=tmp$Subject, Trial=tmp$Trial),
         fun=sd, plot=FALSE)
add_bars(as.numeric(names(acf.se)), y=acf.y+acf.se, y0=acf.y-acf.se, col=NA, border=2, width=.5)
legend('topright', legend="sd", fill=NA, border=2, bty='n')

## ------------------------------------------------------------------------
simdat$Event <- NA
simdat[!is.na(simdat$Y),]$Event <- derive_timeseries(m1AR1)
str(simdat)


## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(itsadug)
infoMessages('off')

## ------------------------------------------------------------------------
library(itsadug)
library(mgcv)
data(simdat)
# select subset of data to reduce processing time:
select <- 1:18
select <- select[select %% 3 ==0]
simdat <- droplevels(simdat[simdat$Subject %in% c(sprintf("a%02d",select), sprintf("c%02d", select)),])

## ------------------------------------------------------------------------
# add start.event and Event columns:
simdat <- start_event(simdat, column="Time", event=c("Subject", "Trial"), label.event="Event")

## ------------------------------------------------------------------------
m1 <- bam(Y ~ Group + s(Time, by=Group) + s(Condition, by=Group, k=5) + ti(Time, Condition, by=Group) + s(Time, Subject, bs='fs', m=1, k=5) + s(Event, bs='re'), data=simdat, discrete=TRUE, method="fREML")

## ------------------------------------------------------------------------
r1 <- start_value_rho(m1)
m1Rho <- bam(Y ~ Group + s(Time, by=Group) + s(Condition, by=Group, k=5) + ti(Time, Condition, by=Group) + s(Time, Subject, bs='fs', m=1, k=5) + s(Event, bs='re'), data=simdat, method="fREML", AR.start=simdat$start.event, rho=r1)

## ------------------------------------------------------------------------
m2Rho <- bam(Y ~ Group + s(Time, by=Group) + s(Condition, by=Group, k=5) + ti(Time, Condition) + s(Time, Subject, bs='fs', m=1, k=5) + s(Event, bs='re'), data=simdat, method="fREML", AR.start=simdat$start.event, rho=r1)

## ------------------------------------------------------------------------
# make sure that info messages are printed to the screen:
infoMessages('on')
compareML(m1Rho, m2Rho)

## ---- include=FALSE------------------------------------------------------
infoMessages('off')

## ------------------------------------------------------------------------
AIC(m1Rho, m2Rho)

## ---- results="asis"-----------------------------------------------------
gamtabs(m1Rho, type="HTML")

## ------------------------------------------------------------------------
simdat$OFGroup <- as.ordered(simdat$Group) 
contrasts(simdat$OFGroup) <- "contr.treatment"
contrasts(simdat$OFGroup)

## ------------------------------------------------------------------------
m1Rho.OF <- bam(Y ~ OFGroup + s(Time) + s(Time, by=OFGroup) + s(Condition, k=5) + s(Condition, by=OFGroup, k=5) + ti(Time, Condition) + ti(Time, Condition, by=OFGroup) + s(Time, Subject, bs='fs', m=1, k=5) + s(Event, bs='re'), data=simdat, method="fREML", AR.start=simdat$start.event, rho=r1)

## ---- results="asis"-----------------------------------------------------
gamtabs(m1Rho.OF, type="HTML")

## ------------------------------------------------------------------------
report_stats(m1Rho.OF)

## ---- fig.width=8, fig.height=4------------------------------------------
par(mfrow=c(1,2))

# PLOT 1:
plot_diff(m1Rho, view="Time", comp=list(Group=c("Adults", "Children")), cond=list(Condition=1), rm.ranef=TRUE, ylim=c(-15,15))
plot_diff(m1Rho, view="Time", comp=list(Group=c("Adults", "Children")),  cond=list(Condition=4), add=TRUE, col='red')
# add legend:
legend('bottom', legend=c("Condition=1", "Condition=4"), col=c(1,2), lwd=1, cex=.75, bty='n')

# PLOT 2:
plot_diff(m1Rho, view="Condition", comp=list(Group=c("Adults", "Children")), cond=list(Time=1000), rm.ranef=TRUE, ylim=c(-15,15))
plot_diff(m1Rho, view="Condition", comp=list(Group=c("Adults", "Children")),  cond=list(Time=2000), add=TRUE, col='red')
# add legend:
legend('bottom', legend=c("Time=1000", "Time=2000"), col=c(1,2), lwd=1, cex=.75, bty='n')

## ---- fig.width=8, fig.height=4------------------------------------------
par(mfrow=c(1,2), cex=1.1)
plot_diff2(m1Rho, view=c("Time", "Condition"), comp=list(Group=c("Adults", "Children")), zlim=c(-15,15), rm.ranef=TRUE)

# with CI:
plot_diff2(m1Rho, view=c("Time", "Condition"), comp=list(Group=c("Adults", "Children")), zlim=c(-15,15), plotCI=TRUE, rm.ranef=TRUE,)


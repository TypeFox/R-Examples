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
simdat$Condition <- as.factor(simdat$Condition)
# Note: this is really not the best fitting model for the data:
m1 <- bam(Y ~ Group * Condition + s(Time) + s(Trial) + ti(Time, Trial) + s(Time, Subject, bs='fs', m=1, k=5), data=simdat)

## ------------------------------------------------------------------------
summary(m1)

## ---- results="asis"-----------------------------------------------------
gamtabs(m1, type="HTML")

## ---- echo=FALSE, fig.show='hide'----------------------------------------
plot_smooth(m1, view="Time", cond=list(Group="Adults"), rug=FALSE)
infoMessages('off')

## ---- fig.width=8, fig.height=4------------------------------------------
par(mfrow=c(1,2), cex=1.1)
# including random effects:
plot_smooth(m1, view="Time", cond=list(Group="Adults"))
# excluding random effects:
plot_smooth(m1, view="Time", cond=list(Group="Adults"), rm.ranef=TRUE)

## ---- fig.width=8, fig.height=4------------------------------------------
par(mfrow=c(1,2), cex=1.1)
# version 1:
plot_smooth(m1, view="Time", plot_all="Group", rm.ranef=TRUE, ylim=c(-15,15))
# version 2:
plot_smooth(m1, view="Time", cond=list(Group="Adults"), rm.ranef=TRUE, rug=FALSE, col="red", ylim=c(-15,15))
plot_smooth(m1, view="Time", cond=list(Group="Children"), rm.ranef=TRUE, rug=FALSE, col="cyan", add=TRUE)

## ---- fig.width=8, fig.height=4------------------------------------------
par(mfrow=c(1,2), cex=1.1)
# including random effects:
fvisgam(m1, view=c("Time", "Trial"), cond=list(Group="Adults"))
# excluding random effects:
fvisgam(m1, view=c("Time", "Trial"), cond=list(Group="Adults"), rm.ranef=TRUE)

## ---- fig.width=8, fig.height=4------------------------------------------
par(mfrow=c(1,2), cex=1.1)
# group children:
fvisgam(m1, view=c("Time", "Trial"), cond=list(Group="Children"), rm.ranef=TRUE, zlim=c(-15,15), main="Children")
# group Adults:
fvisgam(m1, view=c("Time", "Trial"), cond=list(Group="Adults"), rm.ranef=TRUE, zlim=c(-15,15), main="Adults")

## ---- fig.width=12, fig.height=4-----------------------------------------
par(mfrow=c(1,3), cex=1.1)
# all:
plot_parametric(m1,  pred=list(Condition=levels(simdat$Condition),Group=c('Adults', 'Children')), rm.ranef=TRUE)
# children:
plot_parametric(m1,  pred=list(Condition=levels(simdat$Condition)), cond=list(Group="Children"), rm.ranef=TRUE, main="Children")
# adults:
plot_parametric(m1,  pred=list(Condition=levels(simdat$Condition)), cond=list(Group="Adults"), rm.ranef=TRUE, main="Adults")

## ---- fig.width=12, fig.height=4-----------------------------------------
par(mfrow=c(1,3), cex=1.1)
# first smooth of time:
plot(m1, select=1, shade=TRUE, scale=0, ylim=c(-15,15))
abline(h=0)
# second smooth, of trial:
plot(m1, select=2, shade=TRUE, scale=0, ylim=c(-15,15))
abline(h=0)
# third smooth, nonlienar interaction between Time and Trial:
plot(m1, select=3, scale=0, rug=FALSE)

## ---- fig.width=12, fig.height=4-----------------------------------------
par(mfrow=c(1,3), cex=1.1)
# Partial effect of nonlinear interaction between Time and Trial:
plot(m1, select=3, scale=0, rug=FALSE)
# Partial effect of nonlinear interaction between Time and Trial:
pvisgam(m1, select=3, view=c("Time", "Trial"), main="pvisgam", dec=1)
# Summed effect of nonlinear interaction between Time and Trial:
fvisgam(m1, view=c("Time", "Trial"), rm.ranef=TRUE, main="fvisgam", dec=1)

## ---- fig.width=8, fig.height=4------------------------------------------
par(mfrow=c(1,2), cex=1.1)
plot(m1, select=4)
abline(h=0)
inspect_random(m1, select=4)

## ---- fig.width=8, fig.height=4------------------------------------------
par(mfrow=c(1,2), cex=1.1)
# PLOT 1: selection
# ... the adults in black:
inspect_random(m1, select=4, cond=list(Subject=unique(simdat[simdat$Group=="Adults", "Subject"])), col=1, ylim=c(-15, 15), xpd=TRUE)
# ... add children in red:
inspect_random(m1, select=4, cond=list(Subject=unique(simdat[simdat$Group=="Children", "Subject"])), col=2, add=TRUE, xpd=TRUE)

# PLOT 2: averages
# ... of the adults:
inspect_random(m1, select=4, fun=mean, ylim=c(-15,15), cond=list(Subject=unique(simdat[simdat$Group=="Adults", "Subject"])), lwd=2)
# ... of the children:
inspect_random(m1, select=4, fun=mean, cond=list(Subject=unique(simdat[simdat$Group=="Children", "Subject"])), lwd=2, col=2, add=TRUE)

## ---- fig.width=4, fig.height=4------------------------------------------
# All data, clustered by Group (very small dots):
plot_data(m1, view="Time", split_by="Group", cex=.5)

## ---- fig.width=4, fig.height=4------------------------------------------
simdat$Event <- interaction(simdat$Subject, simdat$Trial)
plot_modelfit(m1, view="Time", event=simdat$Event, n=8)

## ---- fig.width=4, fig.height=4------------------------------------------
data(eeg)
# simple GAMM model:
m1 <- gam(Ampl ~ te(Time, X, Y, k=c(10,5,5), 
    d=c(1,2)), data=eeg)
# get electrode postions:
electrodes <- eeg[,c('X','Y','Electrode')]
electrodes <- as.list( electrodes[!duplicated(electrodes),] )
# topo plot, by default uses fvisgam 
# and automatically selects a timestamp (270ms):
plot_topo(m1, view=c("X", "Y"), el.pos=electrodes, dec=2)


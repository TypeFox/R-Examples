
## ----setup, include=FALSE, cache=FALSE-----------------------------------
library(knitr)
options(replace.assign=TRUE,width=50)
opts_chunk$set(fig.path='figs/gph-', cache.path='cache/gph-',
               fig.align='center', dev='pdf', fig.width=5,
               fig.height=5, fig.show='hold', cache=TRUE, par=TRUE)
knit_hooks$set(par=function(before, options, envir){
if (before && options$fig.show!='none') par(mar=c(4,4,1.6,.1),
              cex.lab=.95,cex.axis=.9,mgp=c(2,.7,0),tcl=-.3)
}, crop=hook_pdfcrop)


## ----LOAD-DAAG-----------------------------------------------------------
library(DAAG, warn.conflicts=FALSE)
library(latticeExtra)


## ----nimff, eval=FALSE, echo=FALSE---------------------------------------
## plot(timef~time, data=nihills,
##      xlab="Male record times",
##      ylab="Female record times")
## mftime.lm <- lm(timef ~ time, data=nihills)
## abline(mftime.lm)
## plot(mftime.lm, which=1)


## ----nimff-do, ref.label="nimff", dev='pdf', fig.width=2.75, fig.height=2.85, eval=TRUE, echo=FALSE, out.width="0.47\\textwidth"----


## ----code-nimff, ref.label="nimff", eval=FALSE, echo=TRUE, tidy=FALSE----
## NA


## ----4sims-do-nimff, dev='pdf', fig.width=6.5, fig.height=2.875, eval=TRUE, echo=FALSE----
mftime.lm <- lm(timef ~ time, data=nihills)
gph <- plotSimScat(mftime.lm, layout=c(4,1), show="residuals")
gph <- update(gph, xlab="Record times for males (h)",
              ylab="Record times for females (h)")
print(gph)


## ----4sims-do-nimff, ref.label="4sims-do-nimff", eval=FALSE, echo=TRUE----
## mftime.lm <- lm(timef ~ time, data=nihills)
## gph <- plotSimScat(mftime.lm, layout=c(4,1), show="residuals")
## gph <- update(gph, xlab="Record times for males (h)",
##               ylab="Record times for females (h)")
## print(gph)


## ----diag-logmftime, dev='pdf', fig.width=2.85, fig.height=3, eval=TRUE, echo=FALSE, out.width="0.24\\textwidth"----
plot(mftime.lm, cex.caption=0.8, ask=FALSE)


## ----4sims-mftimesimdiag1, eval=TRUE, echo=FALSE, fig.width=6.25, fig.height=2.85----
plotSimDiags(obj=mftime.lm, which=1, layout=c(4,1))


## ----4sims-mftimesimdiag1-code, ref.label="4sims-mftimesimdiag1", eval=FALSE, echo=TRUE, tidy=FALSE----
## NA


## ----4sims-mftimesimdiag2, dev='pdf', fig.width=6.5, fig.height=2.65, eval=T, echo=FALSE----
plotSimDiags(obj=mftime.lm, which=2, layout=c(4,1))


## ----4sims-mftimesimdiag3, dev='pdf', fig.width=6.5, fig.height=2.65, eval=T, echo=FALSE----
plotSimDiags(obj=mftime.lm, which=3, layout=c(4,1))


## ----4sims-mftimesimdiag5, dev='pdf', fig.width=6.5, fig.height=2.65, eval=T, echo=FALSE----
plotSimDiags(obj=mftime.lm, which=5, layout=c(4,1))


## ----all6, eval=T, echo=T------------------------------------------------
gphs1to6 <- plotSimDiags(obj=mftime.lm, which=1:6, layout=c(4,2))


## ----plot1, eval=T, fig.width=6.5, fig.height=4.2------------------------
update(gphs1to6[[1]], layout=c(4,2))



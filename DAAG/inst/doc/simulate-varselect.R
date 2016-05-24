
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


## ----loadDAAG------------------------------------------------------------
library(DAAG, quietly=TRUE, warn.conflicts=FALSE)


## ----bestset, eval=FALSE, echo=TRUE--------------------------------------
## bestsetNoise(m=100, n=40, nvmax=3)
## bestsetNoise(m=100, n=40, method="backward", nvmax=3)


## ----exhaust, eval=TRUE, echo=FALSE, fig.width=3.25, fig.height=3.25, out.width="0.6\\textwidth"----
## Code
library(quantreg, quietly=TRUE)
library(splines, quietly=TRUE)
set.seed(37)   # Use to reproduce graph that is shown
bsnVaryNvar(m=100, nvar=3:50, nvmax=3)


## ----exhaust, eval=FALSE, echo=TRUE--------------------------------------
## ## Code
## library(quantreg, quietly=TRUE)
## library(splines, quietly=TRUE)
## set.seed(37)   # Use to reproduce graph that is shown
## bsnVaryNvar(m=100, nvar=3:50, nvmax=3)



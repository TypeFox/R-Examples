## ----echo=FALSE,message=FALSE-------------------
library(knitr)
library(poweRlaw)
options(replace.assign=FALSE,width=50)

opts_chunk$set(fig.path='knitr_figure_compare/graphics-', 
               cache.path='knitr_cache_compare/', 
               fig.align='center', 
               dev='pdf', fig.width=5, fig.height=5, 
               fig.show='hold', cache=FALSE, par=TRUE,
               out.width='0.4\\textwidth')
knit_hooks$set(crop=hook_pdfcrop)

knit_hooks$set(par=function(before, options, envir){
  if (before && options$fig.show!='none') {
    par(mar=c(3,3,2,1),cex.lab=.95,cex.axis=.9,
        mgp=c(2,.7,0),tcl=-.01, las=1)
  }}, crop=hook_pdfcrop)

set.seed(1)
palette(c(rgb(170,93,152, maxColorValue=255),
          rgb(103,143,57, maxColorValue=255),
          rgb(196,95,46, maxColorValue=255),
          rgb(79,134,165, maxColorValue=255),
          rgb(205,71,103, maxColorValue=255),
          rgb(203,77,202, maxColorValue=255),
          rgb(115,113,206, maxColorValue=255)))

## -----------------------------------------------
library("poweRlaw")
set.seed(1)
x = rpldis(10000, xmin=2, alpha=2.1)

## -----------------------------------------------
m1 = displ$new(x)
m1$setPars(estimate_pars(m1))

## -----------------------------------------------
m2 = dislnorm$new(x)
m2$setPars(estimate_pars(m2))

## ----F1,echo=1:3,fig.keep='none'----------------
plot(m2, ylab="CDF")
lines(m1)
lines(m2, col=2, lty=2)
grid()

## ----F1,echo=FALSE------------------------------
plot(m2, ylab="CDF")
lines(m1)
lines(m2, col=2, lty=2)
grid()

## -----------------------------------------------
comp = compare_distributions(m1, m2)
comp$p_two_sided

## ----echo=FALSE---------------------------------
compare_distributions(m1, m2)$p_two_sided
compare_distributions(m2, m1)$p_two_sided

## ----echo=FALSE---------------------------------
## We only care if m1 is better than m2
## m1 is clearly better
compare_distributions(m1, m2)$p_one_sided
## m2 isn't better than m1
compare_distributions(m2, m1)$p_one_sided

## -----------------------------------------------
data("moby")

## -----------------------------------------------
m1 = displ$new(moby)
m1$setXmin(estimate_xmin(m1))

## -----------------------------------------------
m2 = dislnorm$new(moby)
m2$setXmin(m1$getXmin())
m2$setPars(estimate_pars(m2))

## ----F2,echo=1:3, fig.keep='none'---------------
plot(m2, ylab="CDF")
lines(m1)
lines(m2, col=2, lty=2)
grid()

## ----F2,echo=FALSE------------------------------
plot(m2, ylab="CDF")
lines(m1)
lines(m2, col=2, lty=2)
grid()

## -----------------------------------------------
comp = compare_distributions(m1, m2)

## -----------------------------------------------
comp$p_two_sided
comp$test_statistic

## ----clean-up, include=FALSE--------------------
# R compiles all vignettes in the same session, which can be bad
rm(list = ls(all = TRUE))


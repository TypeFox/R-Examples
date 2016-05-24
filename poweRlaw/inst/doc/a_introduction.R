## ----echo=FALSE, message=FALSE-------------------------------------------
library(knitr)
library(poweRlaw)
options(replace.assign=FALSE)

opts_chunk$set(fig.path='knitr_figure_poweRlaw/graphics-', 
               cache.path='knitr_cache_poweRlaw/', 
               fig.align='center', 
               dev='pdf', fig.width=5, fig.height=5, 
               fig.show='hold', cache=FALSE, par=TRUE,
               out.width='0.4\\textwidth')
knit_hooks$set(crop=hook_pdfcrop)

knit_hooks$set(par=function(before, options, envir){
  if (before && options$fig.show!='none') {
    par(mar=c(3,4,2,1),cex.lab=.95,cex.axis=.9,
        mgp=c(3,.7,0),tcl=-.01, las=1)
  }}, crop=hook_pdfcrop)

set.seed(1)
palette(c(rgb(170,93,152, maxColorValue=255),
          rgb(103,143,57, maxColorValue=255),
          rgb(196,95,46, maxColorValue=255),
          rgb(79,134,165, maxColorValue=255),
          rgb(205,71,103, maxColorValue=255),
          rgb(203,77,202, maxColorValue=255),
          rgb(115,113,206, maxColorValue=255)))

## ----installation, eval=FALSE--------------------------------------------
#  install.packages("poweRlaw")

## ----eval=FALSE----------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("csgillespie/poweRlaw", subdir="pkg")

## ------------------------------------------------------------------------
library("poweRlaw")

## ----eval=FALSE----------------------------------------------------------
#  help(package="poweRlaw")

## ----results='hide', eval=FALSE------------------------------------------
#  vignette(package="poweRlaw")

## ----results='hide', eval=FALSE------------------------------------------
#  browseVignettes("poweRlaw")

## ----tidy=FALSE, eval=FALSE----------------------------------------------
#  ?displ

## ----eval=FALSE----------------------------------------------------------
#  example(displ)

## ----eval=FALSE----------------------------------------------------------
#  demo(package="poweRlaw")
#  data(package="poweRlaw")

## ----results='hide', eval=FALSE------------------------------------------
#  citation("poweRlaw")

## ----echo=FALSE----------------------------------------------------------
data(bootstrap_moby)
bs = bootstrap_moby
data(bootstrap_p_moby)
bs_p = bootstrap_p_moby

## ----example_word--------------------------------------------------------
data("moby")

## ----fitting-------------------------------------------------------------
m_m = displ$new(moby)

## ----tidy=FALSE----------------------------------------------------------
m_m$getXmin()
m_m$getPars()

## ------------------------------------------------------------------------
m_m$setXmin(5)
m_m$setPars(2)

## ------------------------------------------------------------------------
(est = estimate_pars(m_m))

## ----m_m, echo=1---------------------------------------------------------
(est = estimate_xmin(m_m))
m_m$setXmin(est)

## ----echo=FALSE, fig.width=8, fig.height=8, out.width="0.8\\textwidth"----
par(mfrow=c(2, 2))
plot(m_m, xlab="x", ylab="CDF", 
     pch=21, bg=1, cex=0.6, 
     panel.first=grid())
lines(m_m, col=2, lwd=2)
hist(bs$bootstraps[,2], xlab=expression(x[min]), ylim=c(0, 1600), 
     xlim=c(0, 30), main=NULL, breaks="fd")
grid()
hist(bs$bootstraps[,3], xlab=expression(alpha), 
     ylim=c(0, 500), xlim=c(1.8, 2.1), main=NULL, breaks="fd")
grid()
plot(jitter(bs$bootstraps[,2], factor=1.2), bs$bootstraps[,3], 
     xlab=expression(x[min]), ylab=expression(alpha), 
     xlim=c(0, 30), ylim=c(1.8, 2.1), cex=0.35, 
     pch=21, bg=1, panel.first=grid())

## ----m_m, echo=2, eval=FALSE, results='hide'-----------------------------
#  (est = estimate_xmin(m_m))
#  m_m$setXmin(est)

## ----fig.keep='none'-----------------------------------------------------
## Plot the data (from xmin)
plot(m_m)
## Add in the fitted distribution
lines(m_m, col=2)

## ----fig.keep='none'-----------------------------------------------------
dd = plot(m_m)
head(dd, 3)

## ----uncertainty, eval=FALSE---------------------------------------------
#  bs = bootstrap(m_m, no_of_sims=1000, threads=1)

## ------------------------------------------------------------------------
parallel::detectCores()

## ----fig.keep='none'-----------------------------------------------------
hist(bs$bootstraps[,2], breaks="fd")
hist(bs$bootstraps[,3], breaks="fd")

## ----fig.keep='none'-----------------------------------------------------
plot(jitter(bs$bootstraps[,2], factor=1.2), bs$bootstraps[,3])

## ----do_we_have_a_power,echo=FALSE, fig.width=8, fig.height=4, out.width="\\textwidth"----
par(mfrow=c(1, 3))
hist(bs_p$bootstraps[,2], xlab=expression(x[min]), ylim=c(0, 1600), 
     xlim=c(0, 45), main=NULL, breaks="fd")
grid()
hist(bs_p$bootstraps[,3], xlab=expression(alpha), 
     ylim=c(0, 500), xlim=c(1.80, 2.05), main=NULL, breaks="fd")
grid()

plot(jitter(bs_p$bootstraps[,2], factor=1.2), bs_p$bootstraps[,3], 
     xlab=expression(x[xmin]), ylab=expression(alpha), 
     xlim=c(0, 40), ylim=c(1.8, 2.05), cex=0.35, 
     pch=21, bg=1,
     panel.first= grid())

## ----eval=FALSE, tidy=FALSE----------------------------------------------
#  ## This may take a while
#  ## Use the mle to estimate the parameters
#  bs_p = bootstrap_p(m_m, no_of_sims=1000, threads=2)

## ----distribution_objects------------------------------------------------
m_m = displ$new(moby)

## ----echo=FALSE----------------------------------------------------------
m_m$setXmin(est)

## ----echo=FALSE, results='hide', message=FALSE, warning=FALSE, error=FALSE----
if(!file.exists("blackouts.txt"))
  download.file("http://goo.gl/BsqnP", destfile="blackouts.txt")
blackouts = read.table("blackouts.txt")

## ----loading_data,eval=FALSE---------------------------------------------
#  blackouts = read.table("blackouts.txt")

## ------------------------------------------------------------------------
m_bl = conpl$new(blackouts$V1)

## ----echo=FALSE----------------------------------------------------------
est = estimate_xmin(m_bl)
m_bl$setXmin(est)
plot(m_bl, panel.first=grid())
lines(m_bl, col=2)

## ----echo=FALSE, results='hide', message=FALSE, warning=FALSE, error=FALSE----
if(!file.exists("plfit.R"))
  download.file("http://tuvalu.santafe.edu/~aaronc/powerlaws/plfit.r", destfile="plfit.R")
source("plfit.R")
if(file.exists("plfit_res.rds")) {
  plfit_res = readRDS("plfit_res.rds")
} else {
  plfit_res = plfit(moby)  
  saveRDS(plfit_res, file="plfit_res.rds")
}

## ----eval=FALSE----------------------------------------------------------
#  source("http://tuvalu.santafe.edu/~aaronc/powerlaws/plfit.r")

## ----eval=FALSE----------------------------------------------------------
#  plfit_res = plfit(moby)

## ----results='hide', eval=FALSE------------------------------------------
#  estimate_xmin(m_m, pars=seq(1.5, 2.5, 0.01))

## ----cache=TRUE----------------------------------------------------------
x = rpldis(1000, 1, 2)
plfit(x)

## ----the_ctn_case--------------------------------------------------------
r = rplcon(1000, 10, 2.5)

## ------------------------------------------------------------------------
plfit(r)

## ------------------------------------------------------------------------
m_r = conpl$new(r)
(est = estimate_xmin(m_r))
m_r$setXmin(est)

## ----fig.keep='none'-----------------------------------------------------
plot(m_r)
lines(m_r, col=2)

## ----echo=FALSE----------------------------------------------------------
plot(m_r, ylab="CDF", 
     pch=21, bg=1, cex=0.5, 
     panel.first=grid(), xlim=c(10, 1000))
lines(m_r, col=2, lwd=2)

## ----clean-up, include=FALSE---------------------------------------------
# R compiles all vignettes in the same session, which can be bad
rm(list = ls(all = TRUE))

## ----echo=FALSE, eval=FALSE----------------------------------------------
#  ##Used to generate the figures for github
#  ppi = 50
#  png("../graphics/figure1.png", width=6*ppi, height=4*ppi, res=ppi)
#  setnicepar(mfrow=c(1, 2))
#  plot(m_m, xlab="x", ylab="CDF")
#  lines(m_m, col=2, lty=2)
#  grid()
#  plot(m_bl, xlab="x", ylab="CDF")
#  lines(m_bl, col=2, lty=2)
#  grid()
#  sink = dev.off()
#  
#  png("../graphics/figure2.png", width=6*ppi, height=4*ppi, res=ppi)
#  setnicepar(mfrow=c(1,2))
#  hist(bs$bootstraps[,2], xlab=expression(x[min]), ylim=c(0, 2000),
#       xlim=c(0, 30), main=NULL, breaks="fd")
#  grid()
#  hist(bs$bootstraps[,3], xlab=expression(alpha),
#       ylim=c(0, 500), xlim=c(1.8, 2.1), main=NULL,
#       breaks="fd")
#  grid()
#  sink=dev.off()


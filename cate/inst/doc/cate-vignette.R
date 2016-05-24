## ----setup, include=FALSE, cache=FALSE--------------------------------------------------
library(knitr)
# set global chunk options
opts_chunk$set(fig.path='figure/cate-', fig.align='center', fig.show='hold')
options(formatR.arrow=TRUE,width=90)
Sys.setenv(RSTUDIO_PDFLATEX = "/Library/TeX/texbin/latexmk")

## ----load-------------------------------------------------------------------------------
library(cate)
data(gender.sm)
names(gender.sm)
cbind(X = dim(gender.sm$X), Y = dim(gender.sm$Y), Z = dim(gender.sm$Z)) # matrix dimensions

## ----two-sample-t,cache=TRUE------------------------------------------------------------
t.stats <- apply(gender.sm$Y, 2, function(y, x) t.test(y~x)$statistic, gender.sm$X)

## ----two-sample-t-hist, fig.width=8, fig.height=2.5, echo=FALSE, fig.cap="Histogram of t-statistics before confounder adjustment"----
library(ggplot2)
hist.t.stats <- ggplot() + aes(x = t.stats, y = ..density..) + geom_histogram(binwidth = 0.05 / 2, colour = "black", fill = "white") + xlim(c(-1, 1)) + geom_line(aes(x = seq(-1, 1, 0.01), y = dnorm(seq(-1, 1, 0.01), median(t.stats), mad(t.stats))), col = "violetred3", size = 1) + geom_text(aes(x = 0.5, y = 5, label = paste0("N(", signif(median(t.stats), 2),",", signif(mad(t.stats), 2), "^2)")), col = "violetred3", show_guide= FALSE, size = 3) + xlab("t-statistics") + ylab("density") + theme_bw(base_size = 10)
hist.t.stats

## ----set-seed, include=FALSE------------------------------------------------------------
set.seed(1)

## ----est-confounder-num,cache=TRUE------------------------------------------------------
n <- nrow(gender.sm$Y) # number of samples
gender.data <- data.frame(gender = gender.sm$X, gender.sm$Z)
factor.num <- est.confounder.num(~ gender | . - gender + 0, 
                                 gender.data, gender.sm$Y, 
                                 method = "bcv", bcv.plot = FALSE,
                                 rmax = 30, nRepeat = 20)
factor.num$r

## ---------------------------------------------------------------------------------------
est.confounder.num(~ gender | . - gender + 0, 
                   gender.data, gender.sm$Y, method = "ed")

## ----cate-------------------------------------------------------------------------------
cate.results <- cate(~ gender | . - gender + 0, 
                     gender.data, gender.sm$Y, r = factor.num$r)
names(cate.results)

## ----confounding-test-------------------------------------------------------------------
cate.results$alpha.p.value

## ----candidates, eval=FALSE-------------------------------------------------------------
#  which(p.adjust(cate.results$beta.p.value, "bonferroni") < 0.05) # control FWER at 0.05
#  which(p.adjust(cate.results$beta.p.value, "BH") < 0.2) # control FDR at 0.2

## ----cate-arg---------------------------------------------------------------------------
args(cate)

## ----cate-nc----------------------------------------------------------------------------
cate.results.nc <- cate(~ gender | . - gender + 0, 
                        gender.data, gender.sm$Y, r = factor.num$r,
                        adj.method = "nc", nc = gender.sm$spikectl)

## ----t-hist-after, fig.height=5, fig.width=8, fig.cap="Histograms of test statistics after adjustment", echo=FALSE----
t.stats <- as.vector(cate.results$beta.t)
hist.t.stats.after <- ggplot() + aes(x = t.stats, y = ..density..) + geom_histogram(binwidth = 0.1, colour = "black", fill = "white") + xlim(c(-4, 4)) + geom_line(aes(x = seq(-4, 4, 0.01), y = dnorm(seq(-4, 4, 0.01), median(t.stats), mad(t.stats))), col = "violetred3", size = 1) + geom_text(aes(x = 2, y = 0.3, label = paste0("N(0,1)")), col = "violetred3", show_guide= FALSE, size = 3) + xlab("t-statistics") + ylab("density") + theme_bw(base_size = 10) + ggtitle("Adjusted by Robust Regression")
t.stats.nc <- as.vector(cate.results.nc$beta.t)
hist.t.stats.after.nc <- ggplot() + aes(x = t.stats.nc, y = ..density..) + geom_histogram(binwidth = 0.1, colour = "black", fill = "white") + xlim(c(-4, 4)) + geom_line(aes(x = seq(-4, 4, 0.01), y = dnorm(seq(-4, 4, 0.01), median(t.stats.nc), mad(t.stats.nc))), col = "violetred3", size = 1) + geom_text(aes(x = 2, y = 0.3, label = paste0("N(0,1)")), col = "violetred3", show_guide= FALSE, size = 3) + xlab("t-statistics") + ylab("density") + theme_bw(base_size = 10) + ggtitle("Adjusted by Negative Control")
library(gridExtra)
grid.arrange(hist.t.stats.after, hist.t.stats.after.nc)

## ----fa-demo----------------------------------------------------------------------------
mle <- factor.analysis(gender.sm$Y, r = 5)
names(mle)


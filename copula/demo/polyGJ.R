## Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


### Explore, i.e., plot of poly* (= polyG, polyJ) for all methods ##############
###                                 =====  =====

library(copula)
library(lattice)
## Animation: only if available, and the user does not want to skip it:
do.animation <- require("animation") && (!exists("dont.animate") || !dont.animate)

options(warn=1)

## expected evaluation points for estimating Gumbel and Joe copulas ###########
eep.fun <- function(family, alpha, d, n.MC=5000){
    vapply(alpha, function(alph)
       {
           th <- 1/alph
           cop <- onacopulaL(family, list(th, 1:d))
           U <- rnacopula(n.MC, cop)
           switch(family,
                  "Gumbel" =
              {
                  mean(rowSums(cop@copula@iPsi(U, th))^alph)
              },
                  "Joe" =
              {
                  U. <- (1-U)^th
                  lh <- rowSums(log1p(-U.)) # log(h(..))
                  l1_h <- log(-expm1(lh))
                  mean(exp(lh-l1_h))
              }, stop("wrong family in eep.fun()"))
       }, NA_real_)
}


### compute & plot function for a vector of alphas

plot.poly <- function(family, xlim, ylim, method, alpha, d, n.out = 128,
                      pch = 4, cex = 0.4)
{
    stopifnot(is.numeric(alpha), (len <- length(alpha)) >= 1, is.numeric(d), d == round(d),
              is.numeric(xlim), xlim > 0, is.character(method))
    cols <- colorRampPalette(c("red", "orange", "darkgreen", "turquoise", "blue"),
                             space="Lab")(len)
    switch(family,
	   "Gumbel" = {
	       FUN <- copula:::polyG
	       str <- "G"
	   },
	   "Joe" = {
	       FUN <- copula:::polyJ
	       str <- "J"
	   },
	   stop("wrong 'family'"))
    tit <- paste("poly", str, "(log(x), alpha=..., d=", d,
		 ", log=TRUE, method=\"",method,"\")", sep="")
    xx <- seq(xlim[1], xlim[2], length = n.out)
    lx <- log(xx)
    R <- sapply(alpha, function(ALP)
		FUN(lx, alpha= ALP, d=d, method=method, log=TRUE))
    matplot(xx, R, xlim=xlim, ylim=ylim, main = tit, type = "o", pch=pch, cex=cex,
            xlab="x", ylab=paste("log(poly",str,"(log(x), ...))", sep=""),
            lty = 1, lwd = 1.4, col=cols)
    label <- as.expression(lapply(1:len, function(i)
        substitute(alpha == A, list(A = alpha[i]))))
    legend("bottomright", label, bty="n", lwd=1.4, col=cols, pch=pch, pt.cex=cex)
    invisible(list(f.x = R, x = xx))
}## {plot.poly}


### animation in alpha

if(do.animation)
## animation of poly* functions
## m = number of frames
## d = dimension
## method = method for polyG
poly.ani <- function(family, m, d, method, xlim, ylim)
{
    switch(family,
           "Gumbel" = {
               fun <- copula:::polyG
               str <- "G"
           },
           "Joe" = {
               fun <- copula:::polyJ
               str <- "J"
           },
       {stop("wrong family in plot.poly")})
    alphas <- (1:m)/(m+1) # alphas
    eep <- eep.fun(family, alphas, d) # corresponding expected evaluation points for Gumbel
    x <- seq(xlim[1], xlim[2], length.out=1000)
    lx <- log(x)
    ## Return a list of xyplot objects :
    lapply(1:m, function(i) {
        if(i %% 5 == 1) print(paste(formatC(round(i/m*100), width=3),"% done",sep="")) # progress
        y <- fun(lx, alpha=alphas[i], d=d, method=method, log=TRUE)
        p <- xyplot(y~x, type="l", aspect = 1, xlab= "x",
                    ylab=paste("log(poly",str,"(log(x), ...))",sep=""),
                    xlim=xlim, ylim=ylim, key=
                    list(x=0.35, y=0.1,
                         lines=list(lty=1, col="black"),
                         text=list(paste("expected x-value for alpha=", alphas[i],sep=""))),
                    panel=function(...){
                        panel.xyplot(...)
                        panel.abline(v=eep[i]) # vertical line at expected value
                    }, main=paste("poly",str,"(log(x), alpha=",alphas[i],
                       ", d=",d,", method=",method,", log=TRUE)",sep=""))
        list(y=y, plot=p)
    })
}## {poly.ani}


### Gumbel #####################################################################

### plots for small d

family <- "Gumbel"
polyG <- copula:::polyG
polyG.meths <- eval(formals(polyG)$method, envir = asNamespace("copula"))
##
alpha <- c(0.99, 0.5, 0.01) # alphas; plot largest first, so that all values are visible
xlim <- c(1e-16, 1000)
ylim <- c(-40, 40)
(ev <- eep.fun(family, alpha=alpha, d=5)) # 4.927077 2.462882 1.025498
stopifnot(all(xlim[1] < ev, ev < xlim[2]))
##
(my.polyG.meths <- polyG.meths[!(polyG.meths %in% c("dsumSibuya", "dsSib.RmpfrM"))])
pp5 <- sapply(my.polyG.meths, function(met) {
    r <- plot.poly(family, xlim=xlim, ylim=ylim, method = met, alpha=alpha, d=5)
    Sys.sleep(2)
    r
}, simplify = FALSE)
## => all are fine, even for the much larger range than the expected values: all finite:
t(sapply(pp5, function(L) apply(L$f.x, 2, function(.) sum(!is.finite(.)))))#-> 0 0 0 {no non-finite}


### plots for large d -- quite different picture!

xlim <- c(1e-16, 200)
ylim <- c(300, 600)
(ev <- eep.fun(family, alpha, d=100)) # 96.135490 11.128448  1.060105
stopifnot(all(xlim[1] < ev, ev < xlim[2]))
pp100 <- sapply(my.polyG.meths, function(met) {
    r <- plot.poly(family, xlim=xlim, ylim=ylim, method = met, alpha=alpha, d=100)
    Sys.sleep(2)
    r
}, simplify = FALSE)
## problems -- # NA's:
t(sapply(pp100, function(L) apply(L$f.x, 2, function(.) sum(!is.finite(.)))))
##-> only "default", "direct", and "dsSib.Rmpfr" have no NA's
## but what about the *values*?


### method == "pois"
plot.poly(family, xlim=xlim, ylim=ylim, method="pois", alpha=alpha, d=100)
## => problems for small and moderate alpha

## method == "pois.direct"
plot.poly(family, xlim=xlim, ylim=ylim, method="pois.direct", alpha=alpha, d=100)
## => problems for small and moderate alpha

## method == "stirling"
plot.poly(family, xlim=xlim, ylim=ylim, method="stirling", alpha=alpha, d=100)
## => problems only for large alphas

## method == "stirling.horner"
plot.poly(family, xlim=xlim, ylim=ylim, method="stirling.horner", alpha=alpha, d=100)
## => same as "stirling"

## other methods
plot.poly(family, xlim=xlim, ylim=ylim, method="sort", alpha=alpha, d=100) # log(< 0)
plot.poly(family, xlim=xlim, ylim=ylim, method="horner", alpha=alpha, d=100) # log(< 0)
plot.poly(family, xlim=xlim, ylim=ylim, method="direct", alpha=alpha, d=100) # log(< 0)
plot.poly(family, xlim=xlim, ylim=ylim, method="dsumSibuya", alpha=alpha, d=100) # okay for large alpha *only* ??
## this is *slow* and needs many "recall"s [and of course *looks* good]!
plot.poly(family, xlim=xlim, ylim=ylim, method="dsSib.Rmpfr", alpha=alpha, d=100)


### run time comparison of the methods that worked for some parameter
### ========
set.seed(1)
x <- runif(100000, min=0.01, max=120)
lx <- log(x)

## pois: for large alpha (where it works)
system.time(y.pois <- polyG(lx, alpha=0.99, d=100, method="pois", log=TRUE))
## => 8.91s
stopifnot(all(is.finite(y.pois))) # check

## pois.direct: for large alpha (where it works)
system.time(y.pois.d <- polyG(lx, alpha=0.99, d=100, method="pois.direct", log=TRUE))
## => 6.80s
stopifnot(all(is.finite(y.pois.d))) # check

## stirling: for moderate alpha (where it works)
system.time(y.stirl <- polyG(lx, alpha=0.5, d=100, method="stirling", log=TRUE))
## => 1.92s
stopifnot(all(is.finite(y.stirl))) # check

## stirling.horner: for moderate alpha (where it works)
system.time(y.stirl.Ho <- polyG(lx, alpha=0.5, d=100, method="stirling.horner",
                                log=TRUE))[[1]]
## => 2.79s
stopifnot(all(is.finite(y.stirl.Ho))) # check

## dsumSibuya: for large alpha (where it works)
system.time(y.dsSib.log <- polyG(lx, alpha=0.99, d=100, method="dsSib.log", log=TRUE))
## => 2.28s
stopifnot(all(is.finite(y.dsSib.log))) # check

## conclusion:
## - fastest for large alpha: "dsumSibuya", "pois.direct"
## - fastest for small and moderate alpha: "stirling"
## - further methods tried: pulling out max() for "stirling" => does not increase precision

### check default method

## comparison with Maple (Digits = 100)
v1 <- polyG(log(1), alpha=0.01, d=100, log=TRUE)
v2 <- polyG(log(1), alpha=0.5 , d=100, log=TRUE)
v3 <- polyG(log(1), alpha=0.99, d=100, log=TRUE)
M.v <- c(354.52779560, 356.56733266, 350.99662083)
stopifnot(all.equal(c(v1,v2,v3), M.v))

v1 <- polyG(log(17), alpha=0.01, d=100, log=TRUE)
v2 <- polyG(log(17), alpha=0.5 , d=100, log=TRUE)
v3 <- polyG(log(17), alpha=0.99, d=100, log=TRUE)
M.v <- c(358.15179523, 374.67231305, 370.20372192)
stopifnot(all.equal(c(v1,v2,v3), M.v))

M.v <- c(362.38428102, 422.83827969, 435.36899283)
v1 <- polyG(log(77), alpha=0.01, d=100, log=TRUE)
v2 <- polyG(log(77), alpha=0.5 , d=100, log=TRUE)
v3 <- polyG(log(77), alpha=0.99, d=100, log=TRUE)
stopifnot(all.equal(c(v1,v2,v3), M.v, tolerance=1e-6))


### more detailed graphical precision comparison in d = 100

## dsumSibuya
m <- 49
ylim <- c(200, 700)
polyG.ani.dsumSibuya <- poly.ani(family, m, d=100, method="dsumSibuya",
                                 xlim=c(1e-16,200), ylim=ylim)
if(do.animation)
saveHTML(for(i in 1:m) print(polyG.ani.dsumSibuya[[i]]$plot),
         outdir=file.path(tempdir(),"G_dsumSib"))
## => works for alpha >= 0.75

## pois.direct
polyG.ani.pois.direct <- poly.ani(family, m, d=100, method="pois.direct",
                                  xlim=c(1e-16,200), ylim=ylim)
if(do.animation)
saveHTML(for(i in 1:m) print(polyG.ani.pois.direct[[i]]$plot),
         outdir=file.path(tempdir(),"G_pois.direct"))
## => works for the whole range of *expected* values, esp. for alpha >= 0.72

## stirling
polyG.ani.stirling <- poly.ani(family, m, d=100, method="stirling",
                               xlim=c(1e-16,200), ylim=ylim)
if(do.animation)
saveHTML(for(i in 1:m) print(polyG.ani.stirling[[i]]$plot),
         outdir=file.path(tempdir(),"G_stirling"))
## => works for alpha <= 0.56


## animation in alpha
polyG.ani.default <- poly.ani(family, m, d=100, method="default",
                              xlim=c(1e-16,200), ylim=ylim)
if(do.animation)
saveHTML(for(i in 1:m) print(polyG.ani.default[[i]]$plot),
         outdir=file.path(tempdir(),"G_default"))


### Joe ########################################################################

## plots for small d ###########################################################

family <- "Joe"
polyJ <- copula:::polyJ
alpha <- c(0.05, 0.5, 0.99) # alphas; plot smallest first, so that all values are visible
xlim <- c(1e-16, 1e120)
ylim <- c(0, 1200)
set.seed(1)
(ev <- eep.fun(family, alpha, d=5, n.MC=100000))
## 5.618664e+79 6.923153e+03 5.912850e-02; varies a lot for different runs!
if(!all(xlim[1] < ev, ev < xlim[2])) warning("ev outside xlim")

Jmeths <- eval(formals(polyJ)$method)
Jpp5 <-  sapply(Jmeths, function(met) {
    r <- plot.poly(family, xlim=xlim, ylim=ylim, method = met, alpha=alpha, d=5)
    Sys.sleep(2)
    r
}, simplify = FALSE)
## => "poly" does not work for any reasonable x range -- others ok
t(sapply(Jpp5, function(L) apply(L$f.x, 2, function(.) sum(!is.finite(.)))))

## plots for large d ###########################################################

set.seed(1)
xlim <- c(1e-16, 1e120)
ylim <- c(0, 30000)
system.time(ev <- eep.fun(family, alpha, d=100, n.MC=100000))# longish:
## 15.5 sec (elapsed time)
ev
## 2.119430e+78 8.011466e+02 2.040599e-04; varies a lot for different runs!

if(!all(xlim[1] < ev, ev < xlim[2])) warning("ev outside xlim")
Jpp100 <-  sapply(Jmeths, function(met) {
    r <- plot.poly(family, xlim=xlim, ylim=ylim, method = met, alpha=alpha, d=100)
    Sys.sleep(2)
    r
}, simplify = FALSE)
## => "poly" does not work for any reasonable x range -- others ok
t(sapply(Jpp100, function(L) apply(L$f.x, 2, function(.) sum(!is.finite(.)))))

## i.e. same "message" for d=5 and d=100


### run time comparison of the methods that worked for some parameter ##########
### --------

set.seed(1)
x <- runif(100000, min=0.01, max=1e100)
lx <- log(x)

## log.poly:
system.time(y.log.poly <- polyJ(lx, alpha=0.5, d=100, method="log.poly",
                                log=TRUE))[[1]]
## => 2.701s
stopifnot(all(is.finite(y.log.poly))) # check

## log1p:
system.time(y.log1p <- polyJ(lx, alpha=0.5, d=100,
                             method="log1p", log=TRUE))[[1]]
## => 4.118s
stopifnot(all(is.finite(y.log1p))) # check

## conclusion: use  log.poly  as default

## comparison with Maple (Digits = 100)
v1 <- polyJ(log(1), alpha=0.01, d=100, log=TRUE)
v2 <- polyJ(log(1), alpha=0.5, d=100, log=TRUE)
v3 <- polyJ(log(1), alpha=0.99, d=100, log=TRUE)
M.v <- c(395.73694325, 393.08027226, 386.96715831)# Maple
stopifnot(all.equal(c(v1,v2,v3), M.v))

v1 <- polyJ(log(1e20), alpha=0.01, d=100, log=TRUE)
v2 <- polyJ(log(1e20), alpha=0.5, d=100, log=TRUE)
v3 <- polyJ(log(1e20), alpha=0.99, d=100, log=TRUE)
M.v <- c(4918.2008336, 4915.3815020, 4909.1039909)# Maple
stopifnot(all.equal(c(v1,v2,v3), M.v))

v1 <- polyJ(log(1e100), alpha=0.01, d=100, log=TRUE)
v2 <- polyJ(log(1e100), alpha=0.5, d=100, log=TRUE)
v3 <- polyJ(log(1e100), alpha=0.99, d=100, log=TRUE)
M.v <- c(23154.67477009, 23151.85543852, 23145.57792740)# Maple
stopifnot(all.equal(c(v1,v2,v3), M.v))

## animation in alpha
polyJ.ani.default <- poly.ani(family, m, d=100, method="log.poly",
                              xlim=c(1e-16,1e120), ylim=ylim)
if(do.animation)
saveHTML(for(i in 1:m) print(polyJ.ani.default[[i]]$plot),
         outdir=file.path(tempdir(),"J_log.poly"))
## => rather extreme but seems to be fine

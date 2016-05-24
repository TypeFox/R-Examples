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


require(copula)

#### Testing /  exploring  coeffG(), the coefficients  a_k  for
#### the Gumbel copula's generator derivatives and copula density

### Part 1: Investigate dsumSibuya() ###########################################

## Use dsumSibuya() the way it's used from coeffG() :
dsSib <- function(d, alpha, method, log=TRUE) {
    stopifnot(length(d) == 1, d >= 1, is.numeric(alpha), is.character(method))
    dsumSibuya(x=d, n = 1:d, alpha=alpha, method=method, log=log)
}
(ds.MethsA <- eval(formals(dsumSibuya)$method))
ds.Meths  <- ds.MethsA[ds.MethsA != "RmpfrM"]
(ds.Meths1 <- ds.Meths [ds.Meths  != "Rmpfr"])
dsSibA <- function(d, alpha, methods = ds.Meths, log=TRUE) {
    stopifnot(length(d) == 1, d >= 1, is.numeric(alpha))
    structure(
              vapply(methods, dsumSibuya, FUN.VALUE = rep.int(NA_real_, d),
                     x=d, n = 1:d, alpha=alpha, log=log),
              alpha=alpha, log=log)
}

dsSibMpfr <- function(d, alpha, minPrec = 21, fac.prec = 33/32, log=TRUE, verbose=TRUE) {
    ## fac.prec = 33/32 = 1.0..  ==> try to not waste -- only get "minimally needed" precision
    stopifnot(length(d) == 1, d >= 1, is.numeric(alpha))
    r <- dsumSibuya(x=d, n = 1:d, alpha=alpha, method = "RmpfrM", log=log,
                    mpfr.ctrl = list(minPrec= minPrec, fac = fac.prec, verbose=verbose))
    list(dsumSib = as.numeric(r), prec = getPrec(r))
}

p.dsSib <- function(dsSmat, type="l", ...) {
    stopifnot(is.matrix(dsSmat), (p <- ncol(dsSmat)) >= 1,
              is.numeric(alp <- attr(dsSmat, "alpha")))
    d <- nrow(dsSmat)
    matplot (dsSmat, type=type, xlab = quote(k), ylab = "",
             main = paste("dsumSibuya(x= (d= )", d,
             ", n= (k= ) 1:d, alpha=", formatC(alp),", log = ",attr(dsSmat,"log"),")",
             sep=""), col=1:p, lty=1, ...)
    legend("topright", colnames(dsSmat),
           pch = if(type %in%c("o","b","p")) paste(1:p),
           col=1:p, lty=1, bty="n")
}

str(m50..1 <- dsSibA(50, alpha = 0.1))
stopifnot(is.matrix(m50..1))
p.dsSib(m50..1)## -- now Rmpfr has no NaN anymore!

str(r50.1 <- dsSibMpfr(50, 0.1))# accurate "prec" estimates
plot(r50.1$prec)

str(r50.01 <- dsSibMpfr(50, 0.01))# accurate "prec" estimates
plot(r50.01$prec)

str(r60.01 <- dsSibMpfr(60, 0.01))# accurate "prec" estimates
plot(r60.01$prec)

str(r70.2 <- dsSibMpfr(70, 0.2))
plot(r70.2$prec)

str(r70.4 <- dsSibMpfr(70, 0.4))
plot(r70.4$prec)

str(r80.4 <- dsSibMpfr(80, 0.4))
plot(r80.4$prec)

str(r80.1 <- dsSibMpfr(80, 0.1, verbose=FALSE))
plot(r80.1$prec)

str(r120.1 <- dsSibMpfr(120, 0.1, verbose=FALSE))
plot(r120.1$prec)

system.time(# takes a few minutes!
str(r150.1 <- dsSibMpfr(150, 0.1, verbose=FALSE))
)
plot(r150.1$prec) # up to 600 something

system.time(# takes a few minutes!
str(r150.01 <- dsSibMpfr(150, 0.01, verbose=FALSE))
)
##    user  system elapsed
## 172.563   0.004 173.086
plot(r150.01$prec) # up to ~ 1200

system.time(
str(r100.001 <- dsSibMpfr(100, 0.001, verbose=FALSE))
)
##    user  system elapsed
## 103.226   0.356 103.900
plot(r100.001$prec) # up to ~ 1100

system.time(
str(r90.02 <- dsSibMpfr(90, 0.02, verbose=FALSE))
)
##   user  system elapsed
## 46.479   0.148  46.764
plot(r90.02$prec) # up to ~ 600

system.time(
str(r40.1em4 <- dsSibMpfr(40, 1e-4))
)
 ##   user  system elapsed
 ## 17.745   0.036  17.847
plot(r40.1em4$prec) # up to ~ 600

## quite fast
str(r150.99 <- dsSibMpfr(150, 0.99))
plot(r150.99$prec)# ?? all have prec = 150
plot(r150.99$dsumSib)
## all methods are ok for large alpha :
p.dsSib(dsSibA(150, alpha = 0.99))

if(FALSE) # at D-Math, ETH Zurich:
setwd("/u/maechler/R/Pkgs/copula/demo")
save(list = ls(patt="^r[0-9]"), file = "dsSibMpfr_set.rda")
attach("dsSibMpfr_set.rda")

## TODO:  Data frame / Matrix  -- col.  (d, k, alpha, dsumSib, prec)
## ----- and find "simple" model    prec ~ f(d, k, alpha)


p.dsSib(dsSibA(80, alpha = 0.4), type="b")

p.dsSib(ds70 <- dsSibA(70,  alpha = 0.1))## many more recalls
## more extreme:
p.dsSib(ds1c <- dsSibA(100, alpha = 0.01))## many many more recalls

##--> For small alpha and "large" d,  even Rmpfr0
## is not good enough...  and we may need quite high precision

## Look at the values *before* log(.) :
dsumSibuya(50, 1:50, alpha=0.1, method="Rmpfr")
dd <- dsumSibuya(50, 1:50, alpha=mpfr(0.1, 200), method="Rmpfr")# now -- higher prec. -- "all" fine!
plot(dd, log="y")


### Part 2 #####################################################################

coeffG <- copula:::coeffG


### step (1): look at the a_k's, check if they can be evaluated ################

## Now, explore things seriously :

## use all methods for a set of alpha and d.vec
cGmeths <- function() {
    mm <- eval(formals(coeffG)$method)
    ## for now - do not use deprecated,..
    mm[!(mm %in% c("dsumSibuya", "dsSib.RmpfrM"))]
}

cGmeths()

asN <- function(x, name=deparse(substitute(x))[1]) {
    names(x) <- paste(name, vapply(x, format, ""), sep="=")
    x
}

cG.all <- function(alpha, d.s, meths = NULL) {
    if(is.null(meths)) meths <- cGmeths()
    stopifnot(is.numeric(alpha), is.numeric(d.s), is.character(meths))
    ##
    ## the asN(.) below ensure nice dimnames
    sapply(asN(d.s, "d"), function(d) {
        cat("\nd = ", d,"\n--------\n\n")
        sapply(asN(alpha), function(al) {
            cat("alpha = ", format(al), "\n")
            sapply(meths, coeffG, d=d, alpha=al, log=TRUE)
        }, simplify = "array")
    }, simplify=FALSE)
}

alph.vec <- c(.1, .3, .5, .7, .8, .9, .99, .995)## = 1 - tau
## desirable, but too long for demo:
## d.vec <- c(5,10*(1:10), 20*(6:10))
##--> smaller and fewer for now:
d.vec <- c(5,10*c(1:3, 5, 7, 10, 15))
d.vec <- c(5,10*c(1:3, 5))

options(warn = 1)# show them immediately

ak.all <- cG.all(alpha = alph.vec, d.s = d.vec)
## --> > many warnings, only for  d >= 30
## d =  30, alpha in {0.1, 0.3}
## d >= 40: for all alpha

stopifnot("array" == sapply(ak.all, class))
str(head(ak.all, 3))
ak.all$`d=20`[,,"alpha=0.99"]
## TODO improve the checks, now we have the dsSib.Rmpfr version
chk1 <- function(ak.mat, tol = 1e-7) {
    stopifnot(is.matrix(ak.mat), (d <- nrow(ak.mat)) >= 2)
    n.meth <- ncol(ak.mat)
    med <- apply(ak.mat, 1, median, na.rm=TRUE)
    apply(ak.mat, 2, all.equal, target=med, tolarance=tol)
}

chk1(ak.all$`d=20`[,,"alpha=0.3"])
## chk1(ak.all$`d=90`[,,"alpha=0.3"])

chk.all <-
    lapply(ak.all, function(ak.arr) apply(ak.arr, 3, chk1))
chk.all # quite interesting -- but is the median the truth ??

ak.all$`d=50`[,,"alpha=0.1"] ##--> we see that some "Rmpfr" results got NaN !!!

source(system.file("test-tools-1.R", package="Matrix"))#-> relErr()
relE.ak <- function(ak.mat, tol = 1e-7, trNam = "dsSib.Rmpfr") {
    stopifnot(is.matrix(ak.mat), (d <- nrow(ak.mat)) >= 2,
              is.numeric(true <- ak.mat[, trNam]))
    n.meth <- ncol(ak.mat)
    apply(ak.mat[,trNam != colnames(ak.mat)], 2, relErr, target= true)
}

relE.ak(ak.all$`d=20`[,,"alpha=0.3"])
relE.all <- lapply(ak.all, function(ak.arr) apply(ak.arr, 3, relE.ak))
print(relE.all, digits = 4) # --- wow!


## For d = 5,..85  this is fine (unless for large alpha (!) :

(a.k <- coeffG(100, 0.55, method = "horner"))
## => just works [but in the "extreme area", the numbers are not quite correct,
##    e.g., a.k[53] = 4.325e+83 and Maple says 4.627673570e83]

## conclusion: large alpha's [small theta's] cause the problems!!!
## ==========

## An example showing that for  "dsumSibuya" the problem is exactly *small* alphas:
plot (a.k.H <- coeffG(100, 0.01, method = "horner"), type = "l", lwd=3, log="y")
lines(a.k.J <- coeffG(100, 0.01, method = "dsSib.log"), col=2, type ="o")
lines(a.k.s <- coeffG(100, 0.01, method = "sort"), col=3, type ="l")
lines(a.k.d <- coeffG(100, 0.01, method = "direct"),
      col=adjustcolor("blue"), type ="l", lwd=4)


set.seed(1)

n <- 50
d <- 100
tau <- 0.2
theta <- copGumbel@iTau(tau)
alpha <- 1/theta

## animate this
library(animation)
library(lattice)

m <- 50 # frames
plot.list <- vector("list", m)
alpha.list <-  (1:m)/(m+1)
d <- 100
for(i in 1:m){
    coeffs <- coeffG(d, alpha.list[i], log=TRUE, method = "dsumSibuya")
    plot.list[[i]] <-
        xyplot(coeffs~1:d, type="l", xlim = c(-3,104), ylim = c(-303,374),
               xlab = "k", ylab = expression(log(a[k])), aspect = 1,
               main = substitute(expression(alpha == alpha.),
               list(alpha. = alpha.list[i])))
}
saveHTML(for(i in 1:m) print(plot.list[[i]]))

## conclusion: seems to be better for large alpha [seems to work for alpha >= 0.78,
##             including our test case by a whisker... not totally satisfactory so far]


plot(coeffs~1:100, type="l", xlim = c(-3,104), ylim = c(-303,374),
                             xlab = "k", ylab = expression(log(a[k])), aspect = 1,
                             main = substitute(expression(alpha == alpha.),
                             list(alpha. = alpha.list[i])))

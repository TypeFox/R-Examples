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
source(system.file("Rsource", "tstFit-fn.R", package="copula", mustWork=TRUE))
source(system.file("Rsource", "utils.R",     package="copula", mustWork=TRUE))
##-> assertError(), ... showProc.time()

(doExtras <- copula:::doExtras())


uu <- array(c(9, 7, 8, 3, 2,   4, 1, 5, 6, 10,
              6, 9, 1, 7, 3,   2, 5, 8, 4, 10), dim = c(10L, 2L)) / 11
set.seed(7)
u3 <- cbind(uu, round(runif(10),2))

### t-copula instead of normal -- minimal set for testing here:
## d = 2
(f1 <- fit1(tCopula(df.fixed=TRUE), x = uu))
stopifnot(identical(f1, fit1(tCopula(df.fixed=TRUE),
			     x = data.frame(uu))))# *WITH* a warning
## did not work with data.frame before 2012-08-12

## for df.fixed=FALSE, have 2 parameters ==> cannot use "fit1":
	 (f2.t <- fitCopula(tCopula(), uu, method="itau"))#
	 (f2.r <- fitCopula(tCopula(), uu, method="irho"))#
if(doExtras) {
    print(f2.m <- fitCopula(tCopula(), uu, method=  "ml"))# gives SE for 'df' {from optim()}
    print(f2.M <- fitCopula(tCopula(), uu, method= "mpl"))# no SE for 'df' (for now ..)
}
showProc.time()

## d = 3 : -------------
## ok with df.fixed
tC3f <- tCopula(c(.2,.7, .8), dim=3, dispstr="un", df.fixed=TRUE)
print(f3 <- fitCopula(tC3f, u3, method="itau"))

tC3 <- tCopula(c(.2,.7, .8), dim=3, dispstr="un")
(f3.t <- fitCopula(tC3, u3, method="itau"))
(f3.r <- fitCopula(tC3, u3, method="irho"))
if(doExtras) {
    print(f3.m <- fitCopula(tC3, u3, method=  "ml"))
    print(f3.M <- fitCopula(tC3, u3, method= "mpl"))
}

showProc.time()

set.seed(17)
d <- 5 # dimension
nu <- 4 # degrees of freedom
## define and sample the copula, build pseudo-observations
ec4 <- ellipCopula("t", dim=d, df=nu, df.fixed=TRUE) # <- copula with param NA
(r <- iTau(ec4, tau <- c(0.2, 0.4, 0.6)))
P <- c(r[2], r[1], r[1], r[1], # upper triangle (w/o diagonal) of corr.matrix
             r[1], r[1], r[1],
                   r[3], r[3],
                         r[3])
assertError( setTheta(ec4, value = P) )
## rather need "un" dispersion: Now with smarter tCopula():
(uc4 <- tCopula(dim=d, df=nu, disp = "un", df.fixed=TRUE))
validObject(copt4 <- setTheta(uc4, value = P))
U. <- pobs(rCopula(n=1000, copula=copt4))
pairs(U., gap=0) # => now correct dependency
(cU <- cor(U., method="kendall")) # => correct:
stopifnot(cor(P, cU[lower.tri(cU)]) > 0.99)



### Fitting  multivariate incl margins --- mvdc --------------------------------------
### ===========================================

set.seed(121)
gumbelC <- gumbelCopula(3, dim=2)
gMvGam <- mvdc(gumbelC, c("gamma","gamma"), param = list(list(2,3), list(4,1)))
gMvGam # now nicely show()s -- the *AUTO*-constructed parameter names
stopifnot(identical(gMvGam@paramMargins,
                    list(list(shape = 2, rate = 3),
                         list(shape = 4, rate = 1))))
X <- rMvdc(16000, gMvGam)
plot(X, cex = 1/4)

persp  (gMvGam, dMvdc, xlim = c(0,4), ylim=c(0,8)) ## almost discrete ????
contour(gMvGam, dMvdc, xlim = c(0,2), ylim=c(0,8))
points(X, cex = 1/16, col=adjustcolor("blue", 0.5))

if(FALSE)# unfinished --- TODO maybe move below ('doExtras')!
fMv <- fitMvdc(X, gMvGam)

pFoo <- function(x, lower.tail=TRUE, log.p=FALSE)
     pnorm((x - 5)/20, lower.tail=lower.tail, log.p=log.p)
dFoo <- function(x, lower.tail=TRUE, log.p=FALSE)
     1/20* dnorm((x - 5)/20, lower.tail=lower.tail, log.p=log.p)
qFoo <- qunif

## 'Foo' distribution has *no* parameters:
mv1 <- mvdc(gumbelC, c("gamma","Foo"), param= list(list(3,1), list()))
validObject(mv1)
stopifnot(nrow(R <- rMvdc(3, mv1)) == 3, ncol(R) == 2)
## a wrong way:
assertError(
  mvW <- mvdc(gumbelC, c("gamma","Foo"), param= list(list(3,1), list(NULL)))
)
## must not valid: stopifnot(!isTRUE(validObject(mvW, test=TRUE)))

showProc.time()

## An example which fails (and should)? --
## From: Suzanne Li  @...queensu.ca, Date: Fri, 9 Aug 2013
gumbel <- archmCopula(family = "gumbel",dim = 2)

set.seed(47)# MM {not important, but still want sure reproducibility}
u <- cbind(runif(10),runif(10))
## this now (newly) gives an error:
assertError(fgu <- fitCopula(gumbel, u, method = "ml"))
copGumbel@paraInterval # -> [1, Inf) = exp([0, Inf))
par <- 2^c((0:32)/16, 2+(1:10)/8)
llg <- sapply(par, function(p) loglikCopula(param=p, u, gumbel))
if(dev.interactive()) plot(par, llg, type="b", col=2)
stopifnot(diff(llg) < 0) # so the maximum is for par = 2^0 = 1
## is it just this problem?
## well, the likelihood was not ok for large param=theta; now is:
lrgP <- 100*seq(8, 100, by = 3)
llrg <- vapply(lrgP, function(P) loglikCopula(param = P, u, gumbel), NA_real_)
stopifnot(is.finite(llrg), diff(llrg) < 0, llrg < -11990)## no longer NaN

## Now is it because we *really* should use  elme()  and the "nacopula" families?
## No, this fails too: "outside interval" error {instead of just -Inf !}:
(copG <- onacopulaL("Gumbel", list(NaN,1:2)))
## Estimation -> error for now:
try(efm <- emle(u, copG))


## A simple example with *negative* correlation;
## Want *more helpful* error message here:
set.seed(7)
u1 <- seq(0,1, by=1/128)[2:127]; u2 <- -u1 + round(rnorm(u1)/4,2); u <- pobs(cbind(u1,u2))
plot(u)
msg <- tryCatch(fitCopula(gumbelCopula(), data = u), error=function(e)e$message)
## check error message __FIXME__ want "negative correlation not possible"
## or NO ERROR and a best fit to tau=0 [and the same for other Archimedean families!]
msg



if(!doExtras && !interactive()) q(save="no") ## so the following auto prints
##--------------------------------------------------------------------------

## d = 2 :
## ----- catching fitCopula() error: 'Lapack routine dgesv: system is exactly singular: U[2,2] = 0'
rtx <- tstFit1cop(tCopula(df.fixed=TRUE), tau.set=c(.4, .8),
                  n.set= c(10, 25), N=64)

## for df.fixed=FALSE, have 2 parameters ==> cannot use "fit1":
## ....
## .... TODO

showProc.time()

## The other example  with 'df' (fixed / free):
(tevc <- tevCopula(iTau(tevCopula(), 0.75)))
set.seed(1); str(x <- rCopula(1000, tevc))
plot(x, main = "1000 samples of tevCopula(iTau(tevCopula(), 0.75))")
fitCopula(tevCopula(),		    x, method="irho")# warning
fitCopula(tevCopula(df.fixed=TRUE), x, method="itau")# fine
fitCopula(tevCopula(df.fixed=TRUE), x)# two warnings ==> do not estimate.var:
fitCopula(tevCopula(df.fixed=TRUE), x, estimate.variance=FALSE)
fitCopula(tevCopula(df.fixed=TRUE), x, method="ml")
fitCopula(tevCopula(),		    x)
fitCopula(tevCopula(), 		    x, estimate.variance=FALSE)
try( ## 'df' is not estimated, but it should
fitCopula(tevCopula(), 		    x, method="ml")
)

set.seed(7)
try(
rtevx <- tstFit1cop(tevCopula(, df.fixed=TRUE),
                    tau.set= c(.5, .75), n.set=c(10, 25), N=32)
)##--> singular linear system (Lapack ...)
## now "non-finite finite-difference {in optim()}

## for df.fixed=FALSE, have 2 parameters ==> cannot use "fit1":
## ....
## .... TODO
showProc.time()

data(rdj)
rdj <- rdj[,2:4]
dim(u <- pobs(rdj))# 1262 3
fc <- frankCopula(dim=3)
ffc <- fitCopula(fc, u) ## (failed in 0.999-4 {param constraints})
ffc
summary(ffc)
stopifnot(all.equal(unname(coef(ffc)),
                    2.866564929, tolerance = 1e-5))

showProc.time()


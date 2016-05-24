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


### (Nested) Archimedean Copulas -----------------------------------------------

require(copula)

if(!dev.interactive(orNone=TRUE)) pdf("copula-play.pdf")


### testing psi

myCop <- setTheta(copAMH, value = 0.5) # is maybe more natural

## Care: copula *does* define psi() already!
setGeneric("psi.", function(cop) standardGeneric("psi."))
setMethod(psi., "acopula",
          function(cop) { function(t) cop@psi(t, theta = cop@theta) })
psi.(myCop) # is a function
psi.(myCop)(0:4)
curve(psi.(myCop)(x), 0, 4)
##' but this can also be done directly [ => same curve "on top" :]
curve(myCop@psi(x, theta = myCop@theta),  0, 4, col = 2, add = TRUE)


### testing Kendall's tau

p.Tau <- function(cop, n = 201, xlim = pmin(paraI, 50), ...) {
    stopifnot(is(cop, "acopula"))
    paraI <- cop@paraInterval
    theta <- seq(xlim[1], xlim[2], length.out = n)
    tit <- substitute(tau[NAME](theta), list(NAME = cop@name))
    plot(theta, cop@tau(theta), type = "l", main = tit, ...)
    abline(h = c(0,1), lty = 3, col = "gray20")
}

p.Tau(copAMH)
p.Tau(copClayton)
p.Tau(copFrank, xlim = c(0, 80), ylim= 0:1) # fast via debye_1()
p.Tau(copGumbel)
p.Tau(copJoe, ylim = 0:1, yaxs="i")


### test function ##############################################################

##' @title stopifnot() plus output
##' @param expr
##' @param prefix
##' @param true
##' @return
##' @author Martin Maechler
checkifnot <- function(expr, prefix = "check if", true = "[Ok]")
{
    c0  <- function(...) cat(..., sep = "")
    ## match.call(): not "calling" expr too early:
    c0(prefix, deparse(match.call()[[2]])[1],": ")
    stopifnot(expr)
    c0(true,"\n")
}


##' @title Perform a set of checks on a Archimedean copula object (with theta set)
##' @param cop acopula
##' @param theta1 parameter theta1
##' @param thetavec vector of parameters
##' @param i10 values where psi is evaluated
##' @param nRnd number of generated V0's and V01's
##' @param u01 values where psiinv is evaluated
##' @param lambdaLvec vector of lower tail-dependence coefficients
##' @param lambdaUvec vector of upper tail-dependence coefficients
##' @return list of measurements
##' @author Marius Hofert, Martin Maechler
tstCop <- function(cop, theta1 = cop@theta, thetavec = cop@theta, i10 = 1:10,
                   nRnd = 50, u01 = (1:63)/64, # exact binary fractions
                   lambdaLvec = NA_real_, lambdaUvec = NA_real_)
{
    stopifnot(is(cop, "acopula"))
    cat0 <- function(...) cat(..., "\n", sep = "")
    theta0 <- cop@theta
    CT <- list()

    ### (1) cop name

    cat0(sprintf("(1) copula family: %10s, theta0 = %g",
                 cop@name, theta0))

    ### (2) generator

    ### (2.1) psi and iPsi

    cat("\n(2) values of psi at i10:\n")
    CT <- c(CT, list(psi = system.time(p.i <- cop@psi(i10,theta = theta0))))
    print(p.i)
    checkifnot(identical(numeric(0), cop@iPsi(numeric(0), theta = theta0)))
    checkifnot(cop@iPsi(0, theta = theta0) == Inf)
    cat0("\nvalues of iPsi at u01:")
    CT <- c(CT, list(psiI = system.time(pi.t <-
                     cop@iPsi(u01, theta = theta0))))
    print(pi.t)
    CT[["psiI"]] <- CT[["psiI"]] +
        system.time(pi.pi <- cop@iPsi(p.i,theta = theta0))
    CT[["psi" ]] <- CT[["psi" ]] +
        system.time(p.pit <- cop@psi(pi.t, theta = theta0))
    cat0("check if iPsi(psi(i10))==i10: ", all.equal(pi.pi, i10))
    cat0("check if psi(iPsi(u01))==u01: ", all.equal(p.pit, u01))

    ### (2.2) absdPsi

    ## absdPsi with degree = 10
    cat0("\nvalues of absdPsi with degree=10 at i10:")
    CT <- c(CT, list(absdPsi = system.time(p.D <- cop@absdPsi(i10,theta = theta0,
                     degree = 10))))
    print(p.D)
    cat0("check if all values are nonnegative")
    stopifnot(is.vector(p.D), all(p.D >= 0))
    cat("check absdPsi(Inf,theta,degree=10) = 0 and the class of absdPsi(0,theta,degree=10): ")
    at.0 <- cop@absdPsi(0, theta = theta0, degree = 10)
    stopifnot(cop@absdPsi(Inf, theta = theta0, degree = 10) == 0,
              is.numeric(at.0), !is.nan(at.0))
    cat0("[Ok]")
    ## absdPsi with degree = 10 and MC
    cat("\nvalues of absdPsi with degree=10 and MC at i10:\n")
    CT <- c(CT, list(absdPsi = system.time(p.D <- cop@absdPsi(i10,theta = theta0,
                     degree = 10, n.MC = 1000))))
    print(p.D)
    cat0("check if all values are nonnegative")
    stopifnot(all(p.D >= 0))
    cat("check absdPsi(Inf,theta,degree=10,n.MC=1000) = 0 and the class of absdPsi(0,theta,degree=10,n.MC=1000): ")
    at.0 <- cop@absdPsi(0, theta = theta0, degree = 10, n.MC = 1000)
    stopifnot(cop@absdPsi(Inf, theta = theta0, degree = 10, n.MC = 1000)==0,
              is.numeric(at.0), !is.nan(at.0))
    cat0("[Ok]")

    ### (2.3) absdiPsi

    cat0("\nvalues of absdiPsi at u01:")
    CT <- c(CT, list(absdiPsi. = system.time(absdiPsi. <-
                     cop@absdiPsi(u01, theta = theta0))))
    print(absdiPsi.)
    stopifnot(all(absdiPsi. >= 0, is.numeric(absdiPsi.), !is.nan(absdiPsi.)))
    cat("check the class of absdiPsi(0,theta): ")
    at.0 <- cop@absdiPsi(0, theta = theta0)
    stopifnot(is.numeric(at.0),!is.nan(at.0))
    cat0("[Ok]")

    ### (3) parameter interval

    cat("\n(3) parameter interval:\n")
    print(cop@paraInterval)
    cat0("theta1=",theta1)
    cat0("nesting condition for theta0 and theta1 fulfilled: ",
         cop@nestConstr(theta0,theta1))

    ### (4) V0, dV0, V01, dV01

    ## V0
    CT <- c(CT, list(V0 = system.time(V0 <- cop@V0(nRnd,theta0))))
    cat0("\n(4) ",nRnd," generated V0's:")
    print(summary(V0))
    ## dV0
    cat("\nvalues of dV0 at i10:\n")
    CT <- c(CT, list(dV0 = system.time(dV0.i <- cop@dV0(i10,theta0))))
    print(dV0.i)
    ## V01
    CT <- c(CT, list(V01 = system.time(V01 <- cop@V01(V0,theta0,theta1))))
    cat0("\n",nRnd," generated V01's:")
    print(summary(V01))
    nt <- length(thetavec)
    ## dV01
    cat("\nvalues of dV01 at i10:\n")
    CT <- c(CT, list(dV01 = system.time(dV01.i <- cop@dV01(i10,V0=1,theta0=theta0,
                     theta1=theta1))))
    print(dV01.i)

    ### (5) cCopula {was "cacopula"}

    cat("\n(5) values of cCopula(cbind(v,rev(v)), cop) for v=u01:\n")
    cop. <- onacopulaL(cop@name, list(theta0, 1:2))
    CT <- c(CT, list(cCopula. = system.time(cac <- cCopula(cbind(u01,rev(u01)),
                     cop=cop.))))
    stopifnot(is.vector(cac), length(cac) == length(u01), 0 <= cac, cac <= 1)
    print(cac)

    ### (6) dCopula (log = TRUE) {was dnacopula()}

    u <- matrix(runif(400),ncol=20)
    ocop.2d  <- onacopulaL(cop@name,list(theta0,1:2))
    ocop.20d <- onacopulaL(cop@name,list(theta0,1:20))

    ## d = 2
    cat("\n(6) check dCopula(*, log = TRUE) for u being a random (20x2)-matrix:\n")
    CT <- c(CT, list(dCopula. =
                     system.time(lD <- dCopula(u[,1:2], ocop.2d, log = TRUE))))
    print(lD); stopifnot(is.numeric(lD), is.finite(lD)); cat0("[Ok]")
    cat("check at (0,0.5) and (1,0.5):\n")
    stopifnot(dCopula(cbind(0:1,0.5), ocop.2d, log = FALSE) == 0,
	      dCopula(cbind(0:1,0.5), ocop.2d, log = TRUE ) == -Inf)
    cat0("[Ok]")

    ## d = 20, n.MC = 0
    cat("\n check dCopula(*, log = TRUE) for u being a random (20x20)-matrix:\n")
    CT <- c(CT, list(dCopula. =
                     system.time(lD. <- dCopula(u, ocop.20d, log = TRUE))))
    print(lD.); stopifnot(is.numeric(lD.), is.finite(lD.)); cat0("[Ok]")

    ## d = 20, n.MC > 0
    cat("\n check dCopula(*, log = TRUE) and MC for u being a random (20x20)-matrix:\n")
    CT <- c(CT, list(dCopula. =
                     system.time(lD.. <- dCopula(u, ocop.20d, n.MC = 1000, log = TRUE))))
    print(lD..); stopifnot(is.numeric(lD..), is.finite(lD..)); cat0("[Ok]")

    ## d = 20, check if n.MC > 0 is close to n.MC = 0
    stopifnot(all.equal(lD., lD.., tolerance=0.5))

    ### (7) K

    check.K.u01 <- function(K){
	d.K <- diff(K)
	if(any(neg <- d.K < 0)){ # happens for AMH, Clayton, and Frank (near 1)
            if(any(Neg <- abs(d.K[neg]) > 1e-15* abs(K[-1][neg]))) {
                warning("K(.) is 'substantially' non-monotone for K() / diff(K) =",
                        immediate.=TRUE)
                print(cbind(K = K[-1][Neg], diff.K = d.K[Neg]))
            }
	}
	stopifnot(is.numeric(K), length(K) == length(u01), 0 <= K, K <= 1)
    }

    ## K for d = 2
    cat("\n(7) values of K for d = 2 at u01:\n")
    CT <- c(CT, list(K = system.time(K. <- pK(u01, cop, d = 2))))
    check.K.u01( print(K.) )
    cat("check if K(0) = 0 and K(1) = 1: ")
    stopifnot(pK(0, cop, d = 2)==0,
              pK(1, cop, d = 2)==1)
    cat0("[Ok]")
    ## K for d = 10
    cat("\nvalues of K for d = 10 at u01:\n")
    CT <- c(CT, list(K = system.time(K. <- pK(u01, cop, d = 10))))
    check.K.u01( print(K.) )
    cat("check if K(0) = 0 and K(1) = 1: ")
    stopifnot(pK(0, cop, d = 10)==0,
              pK(1, cop, d = 10)==1)
    cat0("[Ok]")
    ## K for d = 10 and MC
    cat("\nvalues of K for d = 10 and MC at u01:\n")
    CT <- c(CT, list(K = system.time(K. <- pK(u01, cop, d = 10, n.MC = 1000))))
    check.K.u01( print(K.) )
    cat("check if K(0)=0 and K(1)=1: ")
    stopifnot(pK(0, cop, d = 10, n.MC = 1000)==0,
              pK(1, cop, d = 10, n.MC = 1000)==1)
    cat0("[Ok]")

    ### (8) tau, iTau

    cat("\n(8) tau at thetavec:\n")
    CT <- c(CT, list(tau = system.time(ta <- cop@tau(thetavec))))
    print(ta)
    CT <- c(CT, list(tauI = system.time(ta.I <- cop@iTau(ta))))
    cat0("check if iTau(tau(thetavec))==thetavec: ",
         all.equal(ta.I, thetavec))
    lambdaLvec <- rep(as.double(lambdaLvec), length.out= nt)
    lambdaUvec <- rep(as.double(lambdaUvec), length.out= nt)

    ### (9) lambdaL, lambdaLInv

    cat("\n(9) lambdaL at thetavec:\n")
    CT <- c(CT, list(lambdaL = system.time(lT <- cop@lambdaL(thetavec))))
    CT <- c(CT, list(lT.I = system.time(lT.I <- cop@lambdaLInv(lT))))
    print(lT)
    cat0("check if lambdaLInv(lambdaL(thetavec))==lambdaLvec: ",
         all.equal(lT.I, lambdaLvec))

    ### (10) lambdaU, lambdaUInv

    cat("\n(10) lambdaU at thetavec:\n")
    CT <- c(CT, list(lambdaU = system.time(uT <- cop@lambdaU(thetavec))))
    CT <- c(CT, list(uT.I = system.time(uT.I <- cop@lambdaUInv(uT))))
    print(uT)
    cat0("check if lambdaUInv(lambdaU(thetavec))==lambdaUvec: ",
         all.equal(uT.I, lambdaUvec))

    ### (11) dDiag

    cat("\n(11) dDiag at u01 for d=10:\n")
    CT <- c(CT, list(dDiag = system.time(dDiag. <- cop@dDiag(u01, theta=theta0,
                     d=10))))
    print(dDiag.)
    stopifnot(is.numeric(dDiag.), all(dDiag. > 0))
    cat0("[Ok]")

    class(CT) <- "proc_time_list"
    CT
}

##' print() method for the tstCop() results
print.proc_time_list <- function (x, ...) {
    stopifnot(is.list(x), !is.null(nx <- names(x)))
    cat("proc.time()s:                 user system elapsed\n")
    ##    2 4 6 8 0 2 4 6 8 0 2 4 6 89|1 3 |1 3 56|1 3 5 7
    ##            1         2        2
    for(nm in nx)
        if(!all(x[[nm]] == 0, na.rm=TRUE)) {
            ## use 'Time ..' as that works with 'R CMD Rdiff'
            m <- 1000*x[[nm]]
            cat(sprintf("Time [ms] for %13s :%5.0f %6.0f %7.0f\n",
                        ##            2 4 6 8 0 2 4 6 8 0| (20 + (13-4)) = 29
                        nm, m[1], m[2], m[3]))
            ## cat(nm,":\n"); print(x[[nm]], ...)
        }
    invisible(x)
}


### copAMH #####################################################################

myAMH <- setTheta(copAMH, 0.7135001)
thetavec <- c(0.1,0.3,0.5,0.7,0.9)
set.seed(1)
tstCop(myAMH, 0.9429679, thetavec = thetavec)

### copClayton #################################################################

myClayton <- setTheta(copClayton, 0.5)
thetavec <- c(0.5,1,2,5,10)
tstCop(myClayton, 2, thetavec, lambdaL = thetavec, lambdaU = NA)

### copFrank ###################################################################

myFrank <- setTheta(copFrank, 1.860884)
thetavec <- c(0.5,1,2,5,10)
set.seed(11)
tstCop(myFrank, 5.736283, thetavec)

## with a slightly more extensive test:
tau.th <- c(0.055417, 0.11002, 0.21389, 0.4567, 0.66578)
tau.F <- myFrank@tau(thetavec)
stopifnot(all.equal(tau.th, tau.F, tolerance = 0.0001),
          all.equal(.9999, copFrank@tau(copFrank@iTau(0.9999))),
	  all.equal(myFrank@iTau(tau.F, tol = 1e-14), thetavec, tolerance=1e-11))


### copGumbel ##################################################################

myGumbel <- setTheta(copGumbel, 1.25)
thetavec <- c(1,2,4,6,10)
(tG <- tstCop(myGumbel,2, thetavec, lambdaL = NA, lambdaU = thetavec))
u <- seq(0,1, length=32 + 1)[-c(1,32+1)]
u <- as.matrix(expand.grid(u,u))
myGumbel@dacopula(u, theta=1.25)


### copJoe #####################################################################

myJoe <- setTheta(copJoe, 1.25)
thetavec <- c(1.1,2,4,6,10)
set.seed(111)
tstCop(myJoe, 2, thetavec, lambdaL = NA, lambdaU = thetavec)


### Regression tests ------------------------------------

chkPsi <- function(copula, t = c(0, 2^c(-1000,-500, -200,-10*(10:0)), 2:3, 2^(2:40),Inf)) {
    stopifnot(is(copula, "Copula"))
    if(is.unsorted(t)) t <- sort(t)
    psf  <- psi(copula, t)
    ## and also an equidistant t --> to check convexity
    ps.eq <- psi(copula, t. <- seq(0, 20, length=1+2^7))
    stopifnot(is.finite(psf), 0 <= psf, psf <= 1,
              psf[1] == 1, diff(psf) <= 0,
              is.na (pN <- psi(copula, c(NA, NaN))),
              is.nan(pN[2]),
              0 <= ps.eq, ps.eq <= 1, diff(ps.eq) <= 0,
              ## convexity (in light of finite accuracy arithmetic):
              diff(ps.eq, diff=2) >= - 4*.Machine$double.eps *ps.eq[-(1:2)]
              )

    ## for plotting:
    it <- sort.list(tt <- c(t,t.))
    invisible(list(x=tt[it], y= c(psf, ps.eq)[it]))
}

### Negative tau (and dim = 2):

taus <- c(-1,0,1); names(taus) <- paste0("tau=",taus)
taus

## Frank: --------------------------------------------------------
vapply(taus, function(tau) iTau(frankCopula(), tau), 1.)
##    tau=-1     tau=0     tau=1
## -1.81e+16  0.00e+00  7.21e+16
## ~= - Inf        0      + Inf

r <- chkPsi(frankCopula(-2))
plot(r, type="o")
plot(r, type="o", log="xy")
chkPsi(frankCopula( -800))# failed before 2014-06
chkPsi(frankCopula(-2000))# (ditto)
chkPsi(frankCopula(-1e10))# (ditto)


## Clayton: ------------------------------------------------------

vapply(taus, function(tau) iTau(claytonCopula(), tau), 1.)
## tau=-1  tau=0  tau=1
##     -1      0    Inf

stopifnot(all.equal(-2/3, iTau(claytonCopula(), -1/2)))

tools::assertError(chkPsi(claytonCopula(-1.1))) # par. out of bound
chkPsi(claytonCopula(-1)) ## all failed before 2014-05
chkPsi(claytonCopula(-.5))
chkPsi(claytonCopula(-1/8))
chkPsi(claytonCopula(-2^-10))

## AMH:

## --> ../vignettes/rhoAMH-dilog.Rnw

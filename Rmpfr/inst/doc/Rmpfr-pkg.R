### R code from vignette source 'Rmpfr-pkg.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(SweaveHooks= list(fig=function() par(mar=c(5.1, 4.1, 1.1, 2.1))),
        width = 75,
        digits = 7, # <-- here, keep R's default!
        prompt = "R> ",
        continue="   ")
Sys.setenv(LANGUAGE = "en")
if(.Platform$OS.type != "windows")
  Sys.setlocale("LC_MESSAGES","C")


###################################################
### code chunk number 2: diagnose-lib
###################################################
if(nzchar(Sys.getenv("R_MM_PKG_CHECKING"))) print( .libPaths() )
stopifnot(require("sfsmisc"))


###################################################
### code chunk number 3: exp-1
###################################################
exp(1)


###################################################
### code chunk number 4: exp-1-dig-17
###################################################
print(exp(1), digits = 17)


###################################################
### code chunk number 5: exp-1-mp
###################################################
require("Rmpfr") # after having installed the package ...
(one <- mpfr(1, 120))
exp(one)


###################################################
### code chunk number 6: factorial-1
###################################################
ns <- 1:24 ; factorial(ns)


###################################################
### code chunk number 7: factorial-full
###################################################
noquote(sprintf("%-30.0f", factorial(24)))


###################################################
### code chunk number 8: factorial-mpfr
###################################################
ns <- mpfr(1:24, 120) ; factorial(ns)


###################################################
### code chunk number 9: chooseM-ex-fake (eval = FALSE)
###################################################
## chooseMpfr.all(n = 80)


###################################################
### code chunk number 10: chooseM-run
###################################################
capture.and.write(# <- in package 'sfsmisc': ~/R/Pkgs/sfsmisc/R/misc-goodies.R
chooseMpfr.all(n = 80)
                  , 5, 2, middle = 4, i.middle = 13)


###################################################
### code chunk number 11: ex1
###################################################
(0:7) / 7  #  k/7,  for  k= 0..7  printed with R's default precision
options(digits= 16)
(0:7) / 7  #  in full  double precision accuracy
options(digits=  7) # back to default
str(.Machine[c("double.digits","double.eps", "double.neg.eps")], digits=10)
2^-(52:53)


###################################################
### code chunk number 12: n-digs
###################################################
53 * log10(2)


###################################################
### code chunk number 13: ex1
###################################################
x <- mpfr(0:7, 80)/7 # using 80 bits precision
x
7*x
7*x  - 0:7


###################################################
### code chunk number 14: Const-names
###################################################
formals(Const)$name


###################################################
### code chunk number 15: Const-ex
###################################################
Const("pi")
Const("log2")


###################################################
### code chunk number 16: pi-1000
###################################################
system.time(Pi <- Const("pi", 1000 *log2(10)))
Pi


###################################################
### code chunk number 17: pi-fn-Gauss-HB
###################################################
piMpfr <- function(prec=256, itermax = 100, verbose=TRUE) {
    m2 <- mpfr(2, prec) # '2' as mpfr number
    ## -> all derived numbers are mpfr (with precision 'prec')
    p <- m2 + sqrt(m2) # 2 + sqrt(2) = 3.414..
    y <- sqrt(sqrt(m2)) # 2^ {1/4}
    x <- (y+1/y) / m2
    it <- 0L
    repeat {
	p.old <- p
	it <- it+1L
	p <- p * (1+x) / (1+y)
	if(verbose) cat(sprintf("it=%2d, pi^ = %s, |.-.|/|.|=%e\n",
				it, formatMpfr(p, min(50, prec/log2(10))), 1-p.old/p))
	if (abs(p-p.old) <= m2^(-prec))
	    break
	if(it > itermax) {
	    warning("not converged in", it, "iterations") ; break
	}
	## else
	s <- sqrt(x)
	y <- (y*s + 1/s) / (1+y)
	x <- (s+1/s)/2
    }
    p
}

piMpfr()# indeed converges  *quadratically* fast
## with relative error
relErr <- 1 - piMpfr(256, verbose=FALSE) / Const("pi",260)
## in bits :
asNumeric(-log2(abs(relErr)))


###################################################
### code chunk number 18: Math2-def
###################################################
getGroupMembers("Math2")
showMethods("Math2", classes=c("mpfr", "mpfrArray"))


###################################################
### code chunk number 19: round-ex
###################################################
i7 <- 1/mpfr(700, 100)
c(i7, round(i7, digits = 6), signif(i7, digits = 6))


###################################################
### code chunk number 20: roundMpfr-ex
###################################################
roundMpfr(i7, precBits = 30)
roundMpfr(i7, precBits = 15)


###################################################
### code chunk number 21: asNumeric-meth
###################################################
showMethods(asNumeric)


###################################################
### code chunk number 22: format-ex
###################################################
cbind( sapply(1:7, function(d) format(i7, digits=d)) )


###################################################
### code chunk number 23: Math-group
###################################################
getGroupMembers("Math")


###################################################
### code chunk number 24: Matrix-ex
###################################################
head(x <- mpfr(0:7, 64)/7) ; mx <-  x
dim(mx) <- c(4,2)


###################################################
### code chunk number 25: mpfrArr-ex
###################################################
dim(aa <- mpfrArray(1:24, precBits = 80, dim = 2:4))


###################################################
### code chunk number 26: pr-mpfrArr-fake (eval = FALSE)
###################################################
## aa


###################################################
### code chunk number 27: pr-mpfrArr-do
###################################################
capture.and.write(aa, 11, 4)


###################################################
### code chunk number 28: crossprod
###################################################
mx[ 1:3, ] + c(1,10,100)
crossprod(mx)


###################################################
### code chunk number 29: apply-mat
###################################################
apply(7 * mx, 2, sum)


###################################################
### code chunk number 30: Ei-curve
###################################################
getOption("SweaveHooks")[["fig"]]()
curve(Ei,  0, 5, n=2001);  abline(h=0,v=0, lty=3)


###################################################
### code chunk number 31: Li2-1
###################################################
if(mpfrVersion() >= "2.4.0")  ## Li2() is not available in older MPFR versions
  all.equal(Li2(1), Const("pi", 128)^2/6, tol = 1e-30)


###################################################
### code chunk number 32: Li2-curve
###################################################
getOption("SweaveHooks")[["fig"]]()
if(mpfrVersion() >= "2.4.0")
   curve(Li2, -2, 13,   n=2000); abline(h=0,v=0, lty=3)


###################################################
### code chunk number 33: erf-curves
###################################################
getOption("SweaveHooks")[["fig"]]()
curve(erf, -3,3, col = "red", ylim = c(-1,2))
curve(erfc, add = TRUE, col = "blue")
abline(h=0, v=0, lty=3); abline(v=c(-1,1), lty=3, lwd=.8, col="gray")
legend(-3,1, c("erf(x)", "erfc(x)"), col = c("red","blue"), lty=1)


###################################################
### code chunk number 34: integrateR-dnorm
###################################################
integrateR(dnorm,0,2000)
integrateR(dnorm,0,2000, rel.tol=1e-15)
integrateR(dnorm,0,2000, rel.tol=1e-15, verbose=TRUE)


###################################################
### code chunk number 35: integ-exp-double
###################################################
(Ie.d <- integrateR(exp,            0     , 1, rel.tol=1e-15, verbose=TRUE))


###################################################
### code chunk number 36: integ-exp-mpfr
###################################################
(Ie.m <- integrateR(exp, mpfr(0,200), 1, rel.tol=1e-25, verbose=TRUE))
(I.true <- exp(mpfr(1, 200)) - 1)
## with absolute errors
as.numeric(c(I.true - Ie.d$value,
             I.true - Ie.m$value))


###################################################
### code chunk number 37: integ-poly-double
###################################################
if(require("polynom")) {
    x <- polynomial(0:1)
    p <- (x-2)^4 - 3*(x-3)^2
    Fp <- as.function(p)
    print(pI <- integral(p)) # formally
    print(Itrue <- predict(pI, 5) - predict(pI, 0)) ## == 20
} else {
    Fp <- function(x) (x-2)^4 - 3*(x-3)^2
    Itrue <- 20
}
(Id <- integrateR(Fp, 0,      5))
(Im <- integrateR(Fp, 0, mpfr(5, 256),
                  rel.tol = 1e-70, verbose=TRUE))
## and the numerical errors, are indeed of the expected size:
256 * log10(2) # - expect ~ 77 digit accuracy for mpfr(*., 256)
as.numeric(Itrue - c(Im$value, Id$value))



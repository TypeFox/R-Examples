library(robustX)

covNN.1 <- robustX:::covNNC1  ## the original definition (2003)

data(iris)
system.time(cN1 <- covNN.1(iris[-5]))
system.time(cN  <- covNNC (iris[-5]))# faster indeed

## report.and.stop.if.not.all.equal
report.stopifnot.all.eq <- function(a,b, tol, ...) {
    call <- sys.call()
    ae <- all.equal(a,b, tol=tol, ...)
    call[[1]] <- quote(all.equal)
    if(!isTRUE(ae))
	stop(sprintf("Not %s:\n%s\n\n", deparse(call),
		     paste(ae, collapse="\n")),
	     call.=FALSE)
    ## else
    invisible(TRUE)
}

UN <- function(L) lapply(L, unname)

chk.NN.new.old <- function(cNew, cNold, tol = 2e-15, tol.1 = 20*tol) {
    stopifnot(is.list(cNold$innc), length(n.i <- names(cNold$innc)) == 4)
    report.stopifnot.all.eq(UN(cNew [1:4]),
                            UN(cNold[1:4]), tol=tol.1)
    report.stopifnot.all.eq(cNew $innc[n.i],
                            cNold$innc[n.i], tol=tol)
}

try( # testing
    chk.NN.new.old(cN, cN1, tol=0)
)
if(R.version$arch == "x86_64")# very badly fails on 32-bit Windows (why ???)
    chk.NN.new.old(cN, cN1)

## for n = 500, you *do* see it
n <- 500
set.seed(12)
X <- rbwheel(n, 7, spherize=TRUE)

lattice::splom(X, cex=.1)
system.time(cNX1 <- covNN.1(X))# 0.82
system.time(cNX  <- covNNC (X))# 0.66
system.time(cM  <- covMcd (X))# 0.151 - !

try( # testing
    chk.NN.new.old(cNX, cNX1, tol=0)
)
if(R.version$arch == "x86_64")# very badly fails on 32-bit Windows (why ???)
    chk.NN.new.old(cNX, cNX1)

kappa(cM $cov)# 1990.8..
kappa(cNX$cov)#    4.4858
kappa(cov(X)) #    1.0478

## ---- d = 1 :
X1 <- cbind(c(1:6, 1000))

var(X1)
##          [,1]
## [1,] 141861.8
## if 1000 was not an outlier:
var(1:7) ## 4.666667

covNNC(X1)$cov ## -- really not at all robust:
##          [,1]
## [1,] 121595.8

covMcd(X1)$cov
##          [,1]
## [1,] 7.790004

MASS::cov.rob(X1)$cov
##      [,1]
## [1,]  3.5
BACON(X1)$cov
##      [,1]
## [1,]  3.5


if(FALSE) ## FIXME:
    covOGK(X1)$cov

library(Matrix)
library(expm)

source(system.file("test-tools.R", package = "expm"), keep.source=FALSE)
doExtras

tst.sqrtm <- function(m, tol = 1e-12, zap.Im.tol = 1e-10) {
    r.m <- sqrtm(m)## should now work
    s <- r.m %*% r.m
    if(is.complex(s) && all(abs(Im(s)) < mean(abs(s)) * zap.Im.tol))
        s <- Re(s)
    all.equal(m, s, tolerance=tol)
}


options(verbose = TRUE) # -> get some messages from logm.Higham

### ---- Small exact : ----------
 L2 <- cbind(1, 0:1)
lL2 <- cbind(0:1, 0)
(L3 <- rbind(cbind(1,cbind(0:1,0)),1))
(lL3 <- cbind(rbind(0, cbind((2:1)/2,0:1)), 0))

assertError(logm(L2, method="Eigen"))
assertError(logm(L3, method="Eigen"))

logm.Higham08 <- expm:::logm.Higham08
l.L2 <- logm.Higham08(L2)
l.L3 <- logm.Higham08(L3)

all.equal(l.L2, lL2, tolerance=0)# 5.64 e-14 (32-bit *and* 64-bit)
all.equal(l.L3, lL3, tolerance=0)# 2.40 e-15   (ditto)
stopifnot(all.equal(l.L2, lL2, tolerance= 1000e-16),
          all.equal(l.L3, lL3, tolerance=   80e-16))
showProc.time()

### --------- More & larger randomly generated examples : -----------------
set.seed(101)

EA <- expm.Higham08(A <- matrix(round(rnorm(25),1), 5))
all.equal(EA, expm.Higham08(logm.Higham08(EA)), tolerance=0)
## "Mean relative difference: 1.020137e-13"
stopifnot(all.equal(EA, expm.Higham08(logm.Higham08(EA)), tolerance=1e-12))

S <- crossprod(A)
all.equal(S, sqrtm(S) %*% sqrtm(S), tolerance=0)
## "Mean relative difference: 2.26885e-15"
stopifnot(all.equal(S, sqrtm(S) %*% sqrtm(S), tolerance=1e-14))
showProc.time()

set.seed(3)

## n = 50 is already "too" slow (well: logm.Higham08(.) needs 2.2 sec
## --> CPU measurements below
for(n in c(2:5, 10:11, if(doExtras) 30)) {
    cat("n = ",n,": ")
    for(kk in seq_len(if(doExtras) 30 else 10)) {
        ## Testing  logm()
        EA <- expm.Higham08(A <- matrix(round(rnorm(n^2),2), n,n))
        stopifnot(all.equal(EA, expm.Higham08(logm.Higham08(EA)), tolerance=1e-12))
        cat(" ")
        ## Testing  sqrtm() --- for positive definite *and* arbitrary
        stopifnot(tst.sqrtm(A))# A is completely random
        S <- crossprod(A) + rnorm(n^2) / 1000
        stopifnot(tst.sqrtm(S))
        cat(".")
    }
    cat("\n")
}

showProc.time()


### CPU-measurements for  logm()
options(verbose = FALSE)# printing costs ..
set.seed(5)
if(doExtras) {
    n <- 50
    sim <- 32
} else {
    n <- 21
    sim <- 8
}

cpuT <- numeric(sim)
for(k in seq_len(sim)) {
    EA <- expm.Higham08(A <- matrix(rnorm(n^2), n,n))
    cat(".")
    cpuT[k] <- system.time(LEA <- logm.Higham08(EA))[1]
    stopifnot(all.equal(EA, expm.Higham08(LEA), tolerance=1e-12))
}; cat("\n")
summary(cpuT)
## cmath-5 {Feb.2009}; Michi's original code:
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
##   1.794   2.249   2.389   2.388   2.515   2.831

showProc.time()


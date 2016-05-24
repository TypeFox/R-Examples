library(robustbase)
## library(MASS)## MASS::lqs

source(system.file("xtraR/test_LTS.R", package = "robustbase"))
##          ../inst/test_LTS.R

y20 <- c(2:4, 8, 12, 22, 28, 29, 33, 34, 38, 40, 41, 47:48, 50:51, 54, 56, 59)

test_location <- function() {
    ## Improve: print less, and test equality explicitly
    Y <- y20
    print(ltsReg(y=Y))
    print(ltsReg(y=Y, intercept=TRUE))
    print(ltsReg(y=Y, intercept=FALSE))
    print(ltsReg(y=Y, alpha=1))
    print(ltsReg(Y ~ 1))
    print(ltsReg(Y ~ 0))# = Y ~ 1 - 1 :  empty model (no coefficients)
    print(ltsReg(Y ~ 1, alpha=1))
}

test_rsquared <- function() {
    x1 <- y20
    y1 <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 3.5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5)
    ll1 <- ltsReg(x1,y1, alpha = 0.8)
    ## print() ing is platform-dependent, since only ~= 0
    stopifnot(all.equal(unname(coef(ll1)), c(1,0), tolerance=1e-12),
              ll1$scale < 1e-14)
    print(ltsReg(y1,x1, alpha = 0.8))
    print(ltsReg(y1,x1, alpha = 0.8, intercept = FALSE))
}

options(digits = 5)
set.seed(101) # <<-- sub-sampling algorithm now based on R's RNG and seed

doLTSdata()
if(FALSE) { ## FIXME: These *FAIL* !
doLTSdata(nrep = 12, time = FALSE)
doLTSdata(nrep = 12, time = FALSE, method = "MASS")
}

test_rsquared()
test_location()


cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''

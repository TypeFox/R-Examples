context("cvTuning - repeated CV")


## load packages
library("cvTools", quietly=TRUE)

## set seed for reproducibility
set.seed(1234)

## generate data for tests
n <- 20
x <- rnorm(n)
y <- x + rnorm(n)
x <- as.matrix(x)
xy <- data.frame(x, y)

## set up function call to lmrob() and lts()
lmrobCall <- call("lmrob", formula = y~x)
ltsCall <- call("ltsReg", alpha=0.75)

## tuning parameters for lmrob() and lts()
tuning.psi <- c(3.443689, 4.685061)
lmrobTuning <- list(tuning.psi = tuning.psi)
alpha <- c(0.5, 0.75)
ltsTuning <- list(alpha=alpha)

## set up cross-validation folds
K <- 5
R <- 2
folds <- cvFolds(n, K, R)


## run tests

test_that("returned object has class \"cvTuning\" and correct dimensions", {
        ## MM-regression
        lmrobCV <- cvTuning(lmrobCall, data=xy, y=xy$y, tuning=lmrobTuning, 
            cost=rtmspe, folds=folds)
        
        expect_is(lmrobCV, "cvTuning")
        lmrobRTMSPE <- lmrobCV$cv
        expect_is(lmrobRTMSPE, "data.frame")
        expect_equal(dim(lmrobRTMSPE), c(length(tuning.psi), 2))
        lmrobSE <- lmrobCV$se
        expect_is(lmrobSE, "data.frame")
        expect_equal(dim(lmrobSE), c(length(tuning.psi), 2))
        expect_false(any(is.na(lmrobSE)))
        
        ## reweighted and raw LTS
        ltsCV <- cvTuning(ltsCall, x=x, y=y, tuning=ltsTuning, 
            cost=rtmspe, folds=folds, predictArgs=list(fit="both"))
        
        expect_is(ltsCV, "cvTuning")
        ltsRTMSPE <- ltsCV$cv
        expect_is(ltsRTMSPE, "data.frame")
        expect_equal(dim(ltsRTMSPE), c(length(alpha), 3))
        ltsSE <- ltsCV$se
        expect_is(ltsSE, "data.frame")
        expect_equal(dim(ltsSE), c(length(alpha), 3))
        expect_false(any(is.na(ltsSE)))
    })

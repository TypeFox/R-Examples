context("cvFit - repeated CV")


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

## fit models via lmrob() and lts()
lmrobFit <- lmrob(y~x, data=xy)
ltsFit <- ltsReg(x, y, alpha=0.75)
ltsFit$call[[1]] <- as.name("ltsReg")

## set up cross-validation folds
K <- 5
R <- 2
folds <- cvFolds(n, K, R)


## run tests

test_that("returned object has class \"cv\" and correct dimensions", {
        ## MM-regression
        lmrobCV <- cvFit(lmrobFit, data=xy, y=xy$y, cost=rtmspe, folds=folds)
        
        expect_is(lmrobCV, "cv")
        lmrobRTMSPE <- lmrobCV$cv
        expect_is(lmrobRTMSPE, "numeric")
        expect_equal(length(lmrobRTMSPE), 1)
        lmrobReps <- lmrobCV$reps
        expect_is(lmrobReps, "matrix")
        expect_equal(dim(lmrobReps), c(R, 1))
        lmrobSE <- lmrobCV$se
        expect_is(lmrobSE, "numeric")
        expect_equal(length(lmrobSE), 1)
        expect_false(any(is.na(lmrobSE)))
        
        ## reweighted and raw LTS
        ltsCV <- cvFit(ltsFit, x=x, y=y, cost=rtmspe, folds=folds, 
            predictArgs=list(fit="both"))
        
        expect_is(ltsCV, "cv")
        ltsRTMSPE <- ltsCV$cv
        expect_is(ltsRTMSPE, "numeric")
        expect_equal(length(ltsRTMSPE), 2)
        expect_equal(ltsFit$reps, NULL)
        ltsReps <- ltsCV$reps
        expect_is(ltsReps, "matrix")
        expect_equal(dim(ltsReps), c(R, 2))
        ltsSE <- ltsCV$se
        expect_is(ltsSE, "numeric")
        expect_equal(length(ltsSE), 2)
        expect_false(any(is.na(ltsSE)))
    })

context("perryFit - repeated CV")


## load packages
library("perry", quietly=TRUE)

## set seed for reproducibility
set.seed(1234)

## generate data for tests
n <- 20
x <- rnorm(n)
y <- x + rnorm(n)
yy <- cbind(y1=y, y2=y+rnorm(n))
x <- as.matrix(x)
xy <- data.frame(x, y)

## set up cross-validation folds
K <- 5
R <- 2
folds <- cvFolds(n, K, R)

## fit univariate models via lmrob() and lts()
lmrobFit <- lmrob(y~x, data=xy)
ltsFit <- ltsReg(x, y, alpha=0.75)
ltsFit$call[[1]] <- as.name("ltsReg")

## fit multivariate model via lm()
lmFit <- lm(yy~x)


## run tests

test_that("univariate response yields correct \"perry\" object", {
        ## MM-regression
        lmrobCV <- perryFit(lmrobFit, data=xy, y=xy$y, 
            splits=folds, cost=rtmspe)
        
        expect_is(lmrobCV, "perry")
        # check prediction error
        lmrobPE <- lmrobCV$pe
        expect_is(lmrobPE, "numeric")
        expect_equal(length(lmrobPE), 1)
        # check standard error
        lmrobSE <- lmrobCV$se
        expect_is(lmrobSE, "numeric")
        expect_equal(length(lmrobSE), 1)
        expect_false(any(is.na(lmrobSE)))
        # check replications
        lmrobReps <- lmrobCV$reps
        expect_is(lmrobReps, "matrix")
        expect_equal(dim(lmrobReps), c(R, 1))
        # check predictions
        lmrobYHat <- lmrobCV$yHat
        expect_is(lmrobYHat, "list")
        expect_equal(length(lmrobYHat), R)
        
        ## reweighted and raw LTS
        ltsCV <- perryFit(ltsFit, x=x, y=y, splits=folds, 
            predictArgs=list(fit="both"), cost=rtmspe)
        
        expect_is(ltsCV, "perry")
        # check prediction error
        ltsPE <- ltsCV$pe
        expect_is(ltsPE, "numeric")
        expect_equal(length(ltsPE), 2)
        expect_equal(ltsFit$reps, NULL)
        # check standard error
        ltsSE <- ltsCV$se
        expect_is(ltsSE, "numeric")
        expect_equal(length(ltsSE), 2)
        expect_false(any(is.na(ltsSE)))
        # check replications
        ltsReps <- ltsCV$reps
        expect_is(ltsReps, "matrix")
        expect_equal(dim(ltsReps), c(R, 2))
        # check predictions
        ltsYHat <- ltsCV$yHat
        expect_is(ltsYHat, "list")
        expect_equal(length(ltsYHat), R)
    })

test_that("multivariate response yields correct \"perry\" object", {
        ## multivariate LS regression
        lmCV <- perryFit(lmFit, data=lmFit$model, y=yy, splits=folds)
        
        expect_is(lmCV, "perry")
        # check prediction error
        lmPE <- lmCV$pe
        expect_is(lmPE, "numeric")
        expect_equal(length(lmPE), 1)
        # check standard error
        lmSE <- lmCV$se
        expect_is(lmSE, "numeric")
        expect_equal(length(lmSE), 1)
        expect_false(any(is.na(lmSE)))
        # check replications
        lmReps <- lmCV$reps
        expect_is(lmReps, "matrix")
        expect_equal(dim(lmReps), c(R, 1))
        # check predictions
        lmYHat <- lmCV$yHat
        expect_is(lmYHat, "list")
        expect_equal(length(lmYHat), R)
    })

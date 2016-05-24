context("perryTuning - repeated CV")


## load packages
library("perry", quietly=TRUE)

## set seed for reproducibility
set.seed(1234)

## generate data for tests
n <- 20
x <- rnorm(n)
y <- x + rnorm(n)
x <- as.matrix(x)
xy <- data.frame(x, y)

## set up cross-validation folds
K <- 5
R <- 2
folds <- cvFolds(n, K, R)

## set up function call to lmrob() and lts()
lmrobCall <- call("lmrob", formula = y~x)
ltsCall <- call("ltsReg", alpha=0.75)

## tuning parameters for lmrob() and lts()
tuning.psi <- c(3.443689, 4.685061)
lmrobTuning <- list(tuning.psi = tuning.psi)
alpha <- c(0.5, 0.75)
ltsTuning <- list(alpha=alpha)


## run tests

test_that("univariate response yields correct \"perryTuning\" object", {
        ## MM-regression
        lmrobCV <- perryTuning(lmrobCall, data=xy, y=xy$y, tuning=lmrobTuning, 
            splits=folds, cost=rtmspe)
        
        expect_is(lmrobCV, "perryTuning")
        # check prediction error
        lmrobPE <- lmrobCV$pe
        expect_is(lmrobPE, "data.frame")
        expect_equal(dim(lmrobPE), c(length(tuning.psi), 2))
        # check standard error
        lmrobSE <- lmrobCV$se
        expect_is(lmrobSE, "data.frame")
        expect_equal(dim(lmrobSE), c(length(tuning.psi), 2))
        expect_false(any(is.na(lmrobSE)))
        # check replications
        lmrobReps <- lmrobCV$reps
        expect_is(lmrobReps, "data.frame")
        expect_equal(dim(lmrobReps), c(length(tuning.psi)*R, 2))
        # check predictions
        lmrobYHat <- lmrobCV$yHat
        expect_is(lmrobYHat, "list")
        expect_equal(length(lmrobYHat), length(tuning.psi))
        for(yHat in lmrobYHat) {
            expect_is(yHat, "list")
            expect_equal(length(yHat), R)
        }
        
        ## reweighted and raw LTS
        ltsCV <- perryTuning(ltsCall, x=x, y=y, tuning=ltsTuning, splits=folds, 
            predictArgs=list(fit="both"), cost=rtmspe)
        
        expect_is(ltsCV, "perryTuning")
        # check prediction error
        ltsPE <- ltsCV$pe
        expect_is(ltsPE, "data.frame")
        expect_equal(dim(ltsPE), c(length(alpha), 3))
        # check standard error
        ltsSE <- ltsCV$se
        expect_is(ltsSE, "data.frame")
        expect_equal(dim(ltsSE), c(length(alpha), 3))
        expect_false(any(is.na(ltsSE)))
        # check replications
        ltsReps <- ltsCV$reps
        expect_is(ltsReps, "data.frame")
        expect_equal(dim(ltsReps), c(length(tuning.psi)*R, 3))
        # check predictions
        ltsYHat <- ltsCV$yHat
        expect_is(ltsYHat, "list")
        expect_equal(length(ltsYHat), length(tuning.psi))
        for(yHat in ltsYHat) {
            expect_is(yHat, "list")
            expect_equal(length(yHat), R)
        }
    })

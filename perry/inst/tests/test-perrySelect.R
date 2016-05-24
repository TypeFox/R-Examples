context("perrySelect - model selection")


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

## fit models
lmFit <- lm(y~x, data=xy)
lmrobFit <- lmrob(y~x, data=xy)
ltsFit <- ltsReg(x, y, alpha=0.75)
ltsFit$call[[1]] <- as.name("ltsReg")

## cross-validation
lmCV <- perry(lmFit, splits=folds, cost=rtmspe)
lmrobCV <- perry(lmrobFit, splits=folds, cost=rtmspe)
ltsCV <- perry(ltsFit, splits=folds, fit="both", cost=rtmspe)


## run tests

test_that("two \"perry\" objects yield correct \"perrySelect\" object", {
        cv <- perrySelect(LS=lmCV, MM=lmrobCV)
        fits <- c("LS", "MM")
        
        expect_is(cv, "perrySelect")
        expect_equal(fits(cv), factor(fits, levels=fits))
        # check prediction error
        cvPE <- cv$pe
        expect_is(cvPE, "data.frame")
        expect_equal(dim(cvPE), c(length(fits), 2))
        # check standard error
        cvSE <- cv$se
        expect_is(cvSE, "data.frame")
        expect_equal(dim(cvSE), c(length(fits), 2))
        expect_false(any(is.na(cvSE)))
        # check replications
        cvReps <- cv$reps
        expect_is(cvReps, "data.frame")
        expect_equal(dim(cvReps), c(length(fits)*R, 2))
        # check predictions
        cvYHat <- cv$yHat
        expect_is(cvYHat, "list")
        expect_equal(length(cvYHat), length(fits))
        for(yHat in cvYHat) {
            expect_is(yHat, "list")
            expect_equal(length(yHat), R)
        }
    })

test_that("reshaping yields correct \"perrySelect\" object", {
        cv <- perrySelect(LS=lmCV, LTS=ltsCV, .reshape=TRUE)
        fits <- c("LS", "LTS.reweighted", "LTS.raw")
        
        expect_is(cv, "perrySelect")
        expect_equal(fits(cv), factor(fits, levels=fits))
        # check prediction error
        cvPE <- cv$pe
        expect_is(cvPE, "data.frame")
        expect_equal(dim(cvPE), c(length(fits), 2))
        # check standard error
        cvSE <- cv$se
        expect_is(cvSE, "data.frame")
        expect_equal(dim(cvSE), c(length(fits), 2))
        expect_false(any(is.na(cvSE)))
        # check replications
        cvReps <- cv$reps
        expect_is(cvReps, "data.frame")
        expect_equal(dim(cvReps), c(length(fits)*R, 2))
        # check predictions
        cvYHat <- cv$yHat
        expect_is(cvYHat, "list")
        expect_equal(length(cvYHat), length(fits))
        for(yHat in cvYHat) {
            expect_is(yHat, "list")
            expect_equal(length(yHat), R)
        }
    })

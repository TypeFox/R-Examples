context("perryReshape - reshape prediction error results")


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

## fit model
ltsFit <- ltsReg(x, y, alpha=0.75)
ltsFit$call[[1]] <- as.name("ltsReg")

## cross-validation
ltsCV <- perry(ltsFit, splits=folds, fit="both", cost=rtmspe)


## run tests

test_that("reshaping yields correct \"perrySelect\" object", {
        cv <- perryReshape(ltsCV)
        fits <- peNames(ltsCV)
        
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

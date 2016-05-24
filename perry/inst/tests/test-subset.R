context("subset - extracting prediction error results")


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

# 50% subsets
fit50 <- ltsReg(Y ~ ., data = coleman, alpha = 0.5)
cv50 <- perry(fit50, splits = folds, fit = "both", 
    cost = rtmspe, trim = 0.1)

# 75% subsets
fit75 <- ltsReg(Y ~ ., data = coleman, alpha = 0.75)
cv75 <- perry(fit75, splits = folds, fit = "both", 
    cost = rtmspe, trim = 0.1)

# combine results into one object
cv <- perrySelect("0.5" = cv50, "0.75" = cv75)


## run tests

test_that("subsetting \"perry\" object yields correct subset", {
        select <- "reweighted"
        sub <- subset(cv50, select=select)
        
        expect_is(sub, "perry")
        expect_equal(peNames(sub), select)
        # check prediction error
        subPE <- sub$pe
        expect_is(subPE, "numeric")
        expect_equal(length(subPE), length(select))
        # check standard error
        subSE <- sub$se
        expect_is(subSE, "numeric")
        expect_equal(length(subSE), length(select))
        expect_false(any(is.na(subSE)))
        # check replications
        subReps <- sub$reps
        expect_is(subReps, "matrix")
        expect_equal(dim(subReps), c(length(select)*R, 1))
        # check predictions
        subYHat <- sub$yHat
        expect_is(subYHat, "list")
        expect_equal(length(subYHat), R)
        for(yHat in subYHat) expect_is(yHat, "numeric")
    })

test_that("subsetting \"perrySelect\" object yields correct subset", {
        subset <- 1
        select <- "reweighted"
        sub <- subset(cv, subset=subset, select=select)
        
        expect_is(sub, "perrySelect")
        expect_equal(fits(sub), fits(cv)[subset, drop=TRUE])
        expect_equal(peNames(sub), select)
        # check prediction error
        subPE <- sub$pe
        expect_is(subPE, "data.frame")
        expect_equal(dim(subPE), c(length(subset), length(select)+1))
        # check standard error
        subSE <- sub$se
        expect_is(subSE, "data.frame")
        expect_equal(dim(subSE), c(length(subset), length(select)+1))
        expect_false(any(is.na(subSE)))
        # check replications
        subReps <- sub$reps
        expect_is(subReps, "data.frame")
        expect_equal(dim(subReps), c(length(subset)*R, length(select)+1))
        # check predictions
        subYHat <- sub$yHat
        expect_is(subYHat, "list")
        expect_equal(length(subYHat), length(subset))
        for(compYHat in subYHat) {
            expect_is(compYHat, "list")
            expect_equal(length(compYHat), R)
            for(yHat in compYHat) expect_is(yHat, "numeric")
        }
    })

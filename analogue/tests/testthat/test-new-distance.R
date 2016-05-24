## Tests for the new distance compiled code

## load packages
library("testthat")
library("analogue")

context("Testing distance compiled code")

## simple example using dummy data
train <- data.frame(matrix(abs(runif(200)), ncol = 10))
rownames(train) <- LETTERS[1:20]
colnames(train) <- as.character(1:10)
fossil <- data.frame(matrix(abs(runif(100)), ncol = 10))
colnames(fossil) <- as.character(1:10)
rownames(fossil) <- letters[1:10]

## Distance methods to check
METHODS <- c("euclidean", "SQeuclidean","chord", "SQchord",
             "bray", "chi.square", "SQchi.square", "information",
             "chi.distance", "manhattan", "kendall", "gower",
             "alt.gower", "mixed")

## test methods for x and y
test_that("distance matches compiled versions for x and y", {

    ## default settings
    expect_equal(distance(train, fossil),
                 oldDistance(train, fossil))

    ## check all the methods
    for (m in METHODS) {
        ##writeLines(paste("Method:", m))
        expect_equal(distance(train, fossil, method = m),
                     oldDistance(train, fossil, method = m))
    }

})

## test methods for x only
test_that("distance matches compiled versions for x only", {

    ## default settings
    expect_equal(distance(train), oldDistance(train))

    ## check all the methods
    for (m in METHODS) {
        ##writeLines(paste("Method:", m))
        expect_equal(distance(train, method = m),
                     oldDistance(train, method = m))
    }

})

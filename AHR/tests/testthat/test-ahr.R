library(testthat)
library(AHR)

context("ahr")

## TODO: test cases for step.integrate
## TODO: test extreme cases
## TODO: test ahrAJ

test_that("AHR estimates of two-sample problem sum to 1", {
    get.data <- function(n, a, b) {       
        T <- rweibull(n, shape=b, scale=a)
        C <- runif(n, 0, 5)
        
        Y <- pmin(T, C)
        D <- T <= C

        data.frame(Y=Y, D=D, X=rep.int(0, n), V=rep.int(0, n))
    }

    data0 <- get.data(100, 1, 1)
    data1 <- get.data(100, 1, 1.5)

    data <- rbind(data0, data1)
    Z <- rep(1:2, 100)
    x <- ahrWKM(1, Surv(Y, D) ~ Z, data)

    expect_true(sum(x$theta) == 1)
})

test_that("AHR estimates of three-sample problem sum to 1", {
    get.data <- function(n, a, b) {       
        T <- rweibull(n, shape=b, scale=a)
        C <- runif(n, 0, 5)
        
        Y <- pmin(T, C)
        D <- T <= C

        data.frame(Y=Y, D=D, X=rep.int(0, n), V=rep.int(0, n))
    }

    data0 <- get.data(100, 1, 1)
    data1 <- get.data(100, 1, 1.5)
    data2 <- get.data(100, 1, 2)

    data <- rbind(data0, data1, data2)
    Z <- rep(1:3, 100)
    x <- ahrWKM(1, Surv(Y, D) ~ Z, data)

    expect_true(sum(x$theta) == 1)
})

test_that("covariance matrix is symmetric", {
    get.data <- function(n, a, b) {       
        T <- rweibull(n, shape=b, scale=a)
        C <- runif(n, 0, 5)
        
        Y <- pmin(T, C)
        D <- T <= C

        data.frame(Y=Y, D=D, X=rep.int(0, n), V=rep.int(0, n))
    }

    data0 <- get.data(100, 1, 1)
    data1 <- get.data(100, 1, 1.5)
    data2 <- get.data(100, 1, 2)

    data <- rbind(data0, data1, data2)
    Z <- rep(1:3, 100)
    x <- ahrWKM(1, Surv(Y, D) ~ Z, data)

    expect_true(isSymmetric(x$cov.theta))
})

## special case of "rows and columns sum to 0"
test_that("covariance matrix of two-sample problem has special form", {
    get.data <- function(n, a, b) {       
        T <- rweibull(n, shape=b, scale=a)
        C <- runif(n, 0, 5)
        
        Y <- pmin(T, C)
        D <- T <= C

        data.frame(Y=Y, D=D, X=rep.int(0, n), V=rep.int(0, n))
    }

    data0 <- get.data(100, 1, 1)
    data1 <- get.data(100, 1, 1.5)

    data <- rbind(data0, data1)
    Z <- rep(1:2, 100)
    x <- ahrWKM(1, Surv(Y, D) ~ Z, data)

    expect_true(all.equal(x$cov.theta[1,1], x$cov.theta[2,2]))
    expect_true(all.equal(x$cov.theta[1,1], -x$cov.theta[1,2]))
})

test_that("covariance matrix of three-sample problem has rank 2", {
    get.data <- function(n, a, b) {       
        T <- rweibull(n, shape=b, scale=a)
        C <- runif(n, 0, 5)
        
        Y <- pmin(T, C)
        D <- T <= C

        data.frame(Y=Y, D=D, X=rep.int(0, n), V=rep.int(0, n))
    }

    data0 <- get.data(100, 1, 1)
    data1 <- get.data(100, 1, 1.5)
    data2 <- get.data(100, 1, 2)

    data <- rbind(data0, data1, data2)
    Z <- rep(1:3, 100)
    x <- ahrWKM(1, Surv(Y, D) ~ Z, data)
    
    expect_true(all.equal(svd(x$cov.theta)$d[3], 0))
})

test_that("rows and columns of covariance matrix sum to 0", {
    get.data <- function(n, a, b) {       
        T <- rweibull(n, shape=b, scale=a)
        C <- runif(n, 0, 5)
        
        Y <- pmin(T, C)
        D <- T <= C

        data.frame(Y=Y, D=D, X=rep.int(0, n), V=rep.int(0, n))
    }

    data0 <- get.data(100, 1, 1)
    data1 <- get.data(100, 1, 1.5)
    data2 <- get.data(100, 1, 2)

    data <- rbind(data0, data1, data2)
    Z <- rep(1:3, 100)
    x <- ahrWKM(1, Surv(Y, D) ~ Z, data)

    expect_true(all.equal(rep(0, 3), colSums(x$cov.theta)))
    expect_true(all.equal(rep(0, 3), rowSums(x$cov.theta)))
})

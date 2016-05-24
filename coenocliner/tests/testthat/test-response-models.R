## Tests for response model functions

## Load packages
library("testthat")
library("coenocliner")

context("Testing Response Model Functions")

### 1 gradient

## Set up parameters
x <- seq(from = 4, to = 6, length = 100)
px <- list(opt = 4.5, tol = 0.25, h = 20)
G <- Gaussian(x, px = px)

test_that("Gaussian works in single gradient mode", {
    expect_that(G, is_a("numeric"))
    expect_that(length(G), equals(length(x)))
})

test_that("Gaussian throws errors when px have wrong names", {
    names(px) <- c("foo", "tol", "h")
    expect_error(Gaussian(x, px = px))
    names(px) <- c("opt", "foo", "h")
    expect_error(Gaussian(x, px = px))
    names(px) <- c("opt", "tol", "foo")
    expect_error(Gaussian(x, px = px))
    names(px) <- c("foo", "bar", "foobar")
    expect_error(Gaussian(x, px = px))
    ## reset
    names(px) <- c("opt", "tol", "h")
})

test_that("Gaussian throws error if supplied wrong number of px", {
    expect_error(Gaussian(x, px = px[1:2]))
})

test_that("Gaussian throws error on different length params", {
    px <- list(opt = 4.5, tol = 0.25, h = rep(20, 2))
    expect_error(Gaussian(x, px = px))
})

### 2 gradients

## Set up parameters
x <- seq(from = 4, to = 6, length = 100)
y <- seq(from = 35, to = 120, length = 100)
px <- list(opt = 4.5, tol = 0.25, h = 20)
py <- list(opt = 3, tol = 0.5)
G <- Gaussian(x, y, px = px, py = py, corr = 0.5)

test_that("Gaussian works in two gradients mode", {
    expect_that(G, is_a("numeric"))
    expect_that(length(G), equals(length(x)))
})

test_that("Gaussian 2D throws errors when params have wrong names", {
    names(px) <- c("foo", "tol", "h")
    expect_error(Gaussian(x, y, px = px, py = py))
    names(px) <- c("opt", "foo", "h")
    expect_error(Gaussian(x, y, px = px, py = py))
    names(px) <- c("opt", "tol", "foo")
    expect_error(Gaussian(x, y, px = px, py = py))
    names(px) <- c("foo", "bar", "foobar")
    expect_error(Gaussian(x, y, px = px, py = py))

    names(px) <- c("opt", "tol")
    names(py) <- c("foo", "tol")
    expect_error(Gaussian(x, y, px = px, py = py))
    names(py) <- c("opt", "foo")
    expect_error(Gaussian(x, y, px = px, py = py))
    names(py) <- c("foo", "bar")
    expect_error(Gaussian(x, y, px = px, py = py))
    # reset
    names(py) <- c("opt", "tol")
})

test_that("Gaussian 2D throws error if supplied wrong number of params", {
    expect_error(Gaussian(x, y, px = px[1:2], py = py))
    expect_error(Gaussian(x, y, px = px, py = py[1]))
    expect_error(Gaussian(x, y, px = px[1:2], py = py[1]))
    expect_error(Gaussian(x, y, px = px, py = py[c(1,1,2)]))
})

test_that("Gaussian throws error if x and y of different length", {
    expect_error(Gaussian(x, head(y), px = px, py = py))
    expect_error(Gaussian(head(x), y, px = px, py = py))
})

test_that("Gaussian 2D throws error on different length params", {
    expect_error(Gaussian(x, y,
                          px = list(opt = 4.5, tol = 0.25, h = rep(20, 2)),
                          py = py))
    expect_error(Gaussian(x, y,
                          px = px,
                          py = list(opt = 4.5, tol = rep(0.25, 2))))
})

#### Beta

### 1 gradient

## Set up parameters
x <- seq(from = 4, to = 6, length = 100)
px <- list(A0 = 70, m = 10, r = 40, alpha = 2, gamma = 2)
B <- Beta(x, px = px)

test_that("Beta works in single gradient mode", {
    expect_that(B, is_a("numeric"))
    expect_that(length(B), equals(length(x)))
})

test_that("Beta throws errors when px have wrong names", {
    names(px) <- c("foo", "m", "r", "alpha", "gamma")
    expect_error(Beta(x, px = px))
    names(px) <- c("A0", "foo", "r", "alpha", "gamma")
    expect_error(Beta(x, px = px))
    names(px) <- c("A0", "m", "foo", "alpha", "gamma")
    expect_error(Beta(x, px = px))
    names(px) <- c("A0", "m", "r", "foo", "gamma")
    expect_error(Beta(x, px = px))
    names(px) <- c("A0", "m", "r", "alpha", "foo")
    expect_error(Beta(x, px = px))
    names(px) <- c("foo", "bar", "foobar", "Alpha", "Gamma")
    expect_error(Beta(x, px = px))
    ## reset
    names(px) <- c("A0", "m", "r", "alpha", "gamma")
})

test_that("Beta throws error if supplied wrong number of px", {
    expect_error(Beta(x, px = px[1:2]))
})

test_that("Beta throws error on different length params", {
    px <- list(A0 = 70, m = 10, r = 40, alpha = 2, gamma = rep(2, 10))
    expect_error(Beta(x, px = px))
})

### 2 gradients

## Set up parameters
x <- seq(from = 4, to = 6, length = 100)
y <- seq(from = 1, to = 100, length = 100)
px <- list(A0 = 70, m = 10, r = 40, alpha = 2, gamma = 2)
py <- list(m = 10, r = 40, alpha = 2, gamma = 2)
B <- Beta(x, y, px = px, py = py)

test_that("Beta works in two gradients mode", {
    expect_that(B, is_a("numeric"))
    expect_that(length(B), equals(length(x)))
})

test_that("Beta throws errors when px or py have wrong names", {
    names(px) <- c("foo", "m", "r", "alpha", "gamma")
    expect_error(Beta(x, y, px = px, py = py))
    names(px) <- c("A0", "foo", "r", "alpha", "gamma")
    expect_error(Beta(x, y, px = px, py = py))
    names(px) <- c("A0", "m", "foo", "alpha", "gamma")
    expect_error(Beta(x, y, px = px, py = py))
    names(px) <- c("A0", "m", "r", "foo", "gamma")
    expect_error(Beta(x, y, px = px, py = py))
    names(px) <- c("A0", "m", "r", "alpha", "foo")
    expect_error(Beta(x, y, px = px, py = py))
    names(px) <- c("foo", "bar", "foobar", "Alpha", "Gamma")
    expect_error(Beta(x, y, px = px, py = py))
    ## reset
    names(px) <- c("A0", "m", "r", "alpha", "gamma")

    ## test names on py
    names(py) <- c("foo", "r", "alpha", "gamma")
    expect_error(Beta(x, y, py = py, py = py))
    names(py) <- c("m", "foo", "alpha", "gamma")
    expect_error(Beta(x, y, px = px, py = py))
    names(py) <- c("m", "r", "foo", "gamma")
    expect_error(Beta(x, y, px = px, py = py))
    names(py) <- c("m", "r", "alpha", "foo")
    expect_error(Beta(x, y, px = px, py = py))
    names(py) <- c("bar", "foobar", "Alpha", "Gamma")
    expect_error(Beta(x, y, px = px, py = py))
    ## reset
    names(py) <- c("m", "r", "alpha", "gamma")
})

test_that("Beta throws error if supplied wrong number of py", {
    expect_error(Beta(x, y, px = px, py = py[1:2]))
})

test_that("Beta throws error on different length params", {
    px <- list(A0 = 70, m = 10, r = 40, alpha = 2, gamma = rep(2, 10))
    py <- list(m = 10, r = 40, alpha = 2, gamma = rep(2, 10))
    expect_error(Beta(x, y, px = px, py = py))
})

## FIXME: factor this into a set of tests
## x <- seq(3.5, 7, length = 30)
## y <- seq(1, 10, length = 30)
## xy <- expand.grid(x = x, y = y)

## parx <- list(opt = c(5,6), tol = c(0.5,0.3), h = c(50, 75))
## pary <- list(opt = c(5,7), tol = c(1.5, 1.5))

## exx <- expand(xy[,1], opt = c(5), tol = c(0.5), h = c(50))
## exy <- expand(xy[,2], opt = c(5), tol = c(1.5))
## px <- as.list(data.frame(exx[,-1]))
## py <- as.list(data.frame(exy[,-1]))

## args <- list(x = xy[,1], y = xy[,2], px = px, py = py, corr = 0.5)
## tmp <- Gaussian(xy[,1], xy[,2], px = px, py = py, corr = 0.5)
## tmp2 <- do.call("Gaussian", args)

## persp(x, y, matrix(tmp, ncol = length(x)), theta = 45, phi = 30)

## exx <- expand(xy[,1], opt = c(5,6), tol = c(0.5,0.3), h = c(50, 75))
## exy <- expand(xy[,2], opt = c(5,7), tol = c(1.5,1.5))
## px <- as.list(data.frame(exx[,-1]))
## py <- as.list(data.frame(exy[,-1]))

## args <- list(x = xy[,1], y = xy[,2], px = px, py = py, corr = 0.5)
## tmp <- Gaussian(xy[,1], xy[,2], px = px, py = py, corr = 0.5)
## tmp2 <- do.call("Gaussian", args)
## mat <- matrix(tmp, ncol = 2)

## persp(x, y, matrix(mat[,1], ncol = length(x)), theta = 45, phi = 30)
## persp(x, y, matrix(mat[,2], ncol = length(x)), theta = 45, phi = 30)

## sim <- coenocline(xy, params = list(px = parx, py = pary),
##                   responseModel = "gaussian", expectation = TRUE,
##                   extraParams = list(corr = 0.5))

## persp(x, y, matrix(sim[,1], ncol = length(x)), theta = 45, phi = 30)
## persp(x, y, matrix(sim[,2], ncol = length(x)), theta = 45, phi = 30)

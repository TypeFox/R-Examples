## Tests for main coenocline() function

## Load packages
library("testthat")
library("coenocliner")

context("Testing coenocline() functionality")

## set up for tests
x <- seq(from = 4, to = 6, length = 100)
opt <- c(3.75, 4, 4.55, 5, 5.5) + 0.5
tol <- rep(0.25, 5)
h <- rep(20, 5)

## simulate
set.seed(1)
sim <- coenocline(x,
                  responseModel = "gaussian",
                  params = cbind(opt = opt, tol = tol, h = h),
                  countModel = "poisson")

test_that("coenocline() returns an integer matrix", {
    expect_that(sim, is_a("coenocline"))
    expect_that(sim, is_a("matrix"))
    expect_that(typeof(sim) == "integer", is_true())
})

test_that("coenocline() returns matrix with correct number of species", {
    expect_that(NCOL(sim), equals(length(opt)))
})

test_that("coenocline() returns matrix with correct number of samples (rows)", {
    expect_that(NROW(sim), equals(length(x)))
})

## simulate
x <- seq(from = 4, to = 6, length = 100)
y <- seq(from = 1, to = 100, length = 100)
optx <- c(3.75, 4, 4.55, 5, 5.5) + 0.5
opty <- c(5, 50, 75, 10, 60)
tolx <- rep(0.25, 5)
toly <- rep(2, 5)
h <- rep(30, 5)

set.seed(1)
sim <- coenocline(cbind(x, y),
                  responseModel = "gaussian",
                  params = list(px = cbind(opt = optx, tol = tolx, h = h),
                                py = cbind(opt = opty, tol = toly)),
                  countModel = "poisson")

test_that("coenocline() returns an integer matrix with x and y gradients", {
    expect_that(sim, is_a("coenocline"))
    expect_that(sim, is_a("matrix"))
    expect_that(typeof(sim) == "integer", is_true())
})

test_that("coenocline() returns matrix with correct number of species with x and y gradients", {
    expect_that(NCOL(sim), equals(length(opt)))
})

test_that("coenocline() returns matrix with correct number of samples with x and y gradients", {
    expect_that(NROW(sim), equals(length(x)))
})

test_that("coenocline() works with parameters as lists", {

    lp <- list(px = list(opt = optx, tol = tolx, h = h),
               py = list(opt = opty, tol = toly))
    set.seed(1)
    sim2 <- coenocline(cbind(x, y),
                       responseModel = "gaussian",
                       params = lp,
                       countModel = "poisson")

    expect_that(sim2, is_a("coenocline"))
    expect_that(sim2, is_a("matrix"))
    expect_that(NCOL(sim2), equals(length(opt)))
    expect_that(NROW(sim2), equals(length(x)))
    expect_that(sim2, is_identical_to(sim))
})

test_that("coenocline() works with gradient values supplied as a list", {

    lp <- list(px = list(opt = optx, tol = tolx, h = h),
               py = list(opt = opty, tol = toly))
    set.seed(1)
    sim3 <- coenocline(list(x = x, y = y),
                       responseModel = "gaussian",
                       params = lp,
                       countModel = "poisson")

    expect_that(sim3, is_a("coenocline"))
    expect_that(sim3, is_a("matrix"))
    expect_that(NCOL(sim3), equals(length(opt)))
    expect_that(NROW(sim3), equals(length(x)))
    expect_that(sim3, is_identical_to(sim))
})

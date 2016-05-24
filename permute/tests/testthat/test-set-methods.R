library("testthat")
library("permute")

context("Testing set<- methods")

test_that("default methods for get functions", {
    v <- 1:10
    expect_error(setNperm(v) <- 10, regexp = "No default method")
    expect_error(setMaxperm(v) <- 10, regexp = "No default method")
    expect_error(setMinperm(v) <- 10, regexp = "No default method")
    expect_error(setComplete(v) <- TRUE, regexp = "No default method")
    expect_error(setBlocks(v) <- TRUE, regexp = "No default method")
    expect_error(setAllperms(v) <- 1:10, regexp = "No default method")
    expect_error(setObserved(v) <- TRUE, regexp = "No default method")
    expect_error(setPlots(v) <- Plots(), regexp = "No default method")
    expect_error(setWithin(v) <- Within(), regexp = "No default method")
    expect_error(setStrata(v) <- gl(2, 5), regexp = "No default method")
    expect_error(setRow(v) <- 4, regexp = "No default method")
    expect_error(setCol(v) <- 5, regexp = "No default method")
    expect_error(setDim(v) <- c(2,3), regexp = "No default method")
    expect_error(setType(v) <- "series", regexp = "No default method")
    expect_error(setMirror(v) <- TRUE, regexp = "No default method")
    expect_error(setMake(v) <- TRUE, regexp = "No default method")
    expect_error(setConstant(v) <- "series", regexp = "No default method")
})

test_that("how() set methods throw errors where not appropriate for use", {
    h <- how()
    expect_error(setDim(h) <- c(2,3), regexp = "can not be used directly on '\"how\"' objects")
    expect_error(setMirror(h) <- TRUE, regexp = "can not be used directly on '\"how\"' objects")
    expect_error(setConstant(h) <- TRUE, regexp = "can not be used directly on '\"how\"' objects")
    expect_error(setType(h) <- "series", regexp = "can not be used directly on '\"how\"' objects")
    expect_error(setRow(h) <- 2, regexp = "can not be used directly on '\"how\"' objects")
    expect_error(setCol(h) <- 3, regexp = "can not be used directly on '\"how\"' objects")
})

test_that("set within for class how works", {
    h <- how()
    setWithin(h) <- Within(type = "series")
    expect_identical(getType(h, which = "within"), "series")
})

test_that("set plots for class how works", {
    h <- how()
    setPlots(h) <- Plots(type = "series")
    expect_identical(getType(h, which = "plots"), "series")
})

test_that("test setMinperm work", {
    h <- how()
    setMinperm(h) <- 999
    expect_is(h, "how")
    expect_equal(getMinperm(h), 999)
})

test_that("test setMaxperm work", {
    h <- how()
    nperm <- 99999
    setMaxperm(h) <- nperm
    expect_is(h, "how")
    expect_equal(getMaxperm(h), nperm)
})

test_that("test setStrata works", {
    h <- how()
    f <- gl(5,5)
    setStrata(h) <- f
    expect_is(h, "how")
    expect_identical(getStrata(h), f)

    plots <- getPlots(h)
    f <- gl(4,5)
    setStrata(plots) <- f
    expect_identical(getStrata(plots), f)
    setPlots(h) <- plots
    expect_identical(getStrata(h), f)
})

test_that("test setRow<- works", {
    f <- gl(9, 25)
    h <- how(within = Within(type = "grid", nrow = 5, ncol = 5),
             plots  =  Plots(type = "grid", nrow = 3, ncol = 3,
                             strata = f))
    w <- getWithin(h)
    setRow(w) <- 2
    expect_identical(getRow(w), 2L)
    expect_is(w, "Within")

    p <- getPlots(h)
    setRow(p) <- 4
    expect_identical(getRow(p), 4L)
    expect_is(p, "Plots")
})

test_that("test setCol<- works", {
    f <- gl(9, 25)
    h <- how(within = Within(type = "grid", nrow = 5, ncol = 5),
             plots  =  Plots(type = "grid", nrow = 3, ncol = 3,
                             strata = f))
    w <- getWithin(h)
    setCol(w) <- 2
    expect_identical(getCol(w), 2L)
    expect_is(w, "Within")

    p <- getPlots(h)
    setCol(p) <- 4
    expect_identical(getCol(p), 4L)
    expect_is(p, "Plots")
})

test_that("test setDim<- works", {
    f <- gl(9, 25)
    h <- how(within = Within(type = "grid", nrow = 5, ncol = 5),
             plots  =  Plots(type = "grid", nrow = 3, ncol = 3,
                             strata = f))
    w <- getWithin(h)
    setDim(w) <- c(2, 4)
    expect_identical(getDim(w), c(2L, 4L))
    expect_is(w, "Within")

    p <- getPlots(h)
    setDim(p) <- c(4, 3)
    expect_identical(getDim(p), c(4L, 3L))
    expect_is(p, "Plots")
})

test_that("test setType<- works", {
    f <- gl(9, 25)
    h <- how(within = Within(type = "grid", nrow = 5, ncol = 5),
             plots  =  Plots(type = "series", strata = f))
    w <- getWithin(h)
    setType(w) <- "free"
    expect_is(w, "Within")
    expect_identical(getType(w), "free")

    p <- getPlots(h)
    expect_error(setType(p) <- "strata",
                 regexp = "Invalid permutation type")
    setType(p) <- "none"
    expect_is(p, "Plots")
    expect_identical(getType(p), "none")
})

test_that("test setMirror<- works", {
    f <- gl(9, 25)
    h <- how(within = Within(type = "grid", nrow = 5, ncol = 5),
             plots  =  Plots(type = "series", strata = f))
    w <- getWithin(h)
    setMirror(w) <- TRUE
    expect_is(w, "Within")
    expect_true(getMirror(w))

    p <- getPlots(h)
    setMirror(p) <- TRUE
    expect_is(p, "Plots")
    expect_true(getMirror(p))
})

test_that("test setConstant<- works", {
    f <- gl(9, 25)
    h <- how(within = Within(type = "grid", nrow = 5, ncol = 5),
             plots  =  Plots(type = "series", strata = f))
    w <- getWithin(h)
    setConstant(w) <- TRUE
    expect_is(w, "Within")
    expect_true(getConstant(w))

    p <- getPlots(h)
    expect_error(setConstant(p) <- TRUE,
                 regexp = "setConstant` does not apply to '\"Plots\"' objects.")
})

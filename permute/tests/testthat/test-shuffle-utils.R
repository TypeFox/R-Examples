library("testthat")
library("permute")

context("Testing shuffle-util functions")

test_that("shuffelStrata works for various 'type's", {
    f <- gl(4,5)
    perm <- shuffleStrata(f, type = "free")
    expect_equal(length(perm), length(f))
    expect_true(all(perm <= 20))
    expect_true(all(perm >= 1))
    perm <- shuffleStrata(f, type = "series")
    expect_equal(length(perm), length(f))
    expect_true(all(perm <= 20))
    expect_true(all(perm >= 1))
    perm <- shuffleStrata(f, type = "grid", ncol = 2, nrow = 2)
    expect_equal(length(perm), length(f))
    expect_true(all(perm <= 20))
    expect_true(all(perm >= 1))
    expect_error(shuffleStrata(f, type = "unknown"),
                 regexp = "Invalid permutation type")
})

test_that("shuffleGrid works with only ncol and nrow args", {
    nr <- 5
    nc <- 5
    perm <- shuffleGrid(ncol = nc, nrow = nr)
    expect_equal(length(perm), prod(nr, nc))
})

test_that("shuffleGrid works with mirror & flip combinations", {
    nr <- 5
    nc <- 5
    perm <- shuffleGrid(ncol = nc, nrow = nr, mirror = TRUE)
    expect_equal(length(perm), prod(nr, nc))
    perm <- shuffleGrid(ncol = nc, nrow = nr, mirror = TRUE)
    expect_equal(length(perm), prod(nr, nc))
    perm <- shuffleGrid(ncol = nc, nrow = nr, mirror = TRUE,
                        flip = c(TRUE, FALSE))
    expect_equal(length(perm), prod(nr, nc))
    perm <- shuffleGrid(ncol = nc, nrow = nr, mirror = TRUE,
                        flip = c(FALSE, TRUE))
    expect_equal(length(perm), prod(nr, nc))
})

test_that("shuffleGrid works with start.row & start.col", {
    nr <- 5
    nc <- 5
    perm <- shuffleGrid(ncol = nc, nrow = nr,
                        start.row = 2)
    expect_equal(length(perm), prod(nr, nc))
    expect_true(all(perm <= 25))
    expect_true(all(perm >=1))
    perm <- shuffleGrid(ncol = nc, nrow = nr,
                        start.col = 4)
    expect_equal(length(perm), prod(nr, nc))
    expect_true(all(perm <= 25))
    expect_true(all(perm >=1))
    perm <- shuffleGrid(ncol = nc, nrow = nr,
                        start.col = 4, start.row = 2)
    expect_equal(length(perm), prod(nr, nc))
    expect_true(all(perm <= 25))
    expect_true(all(perm >=1))
})

library("testthat")
library("permute")

context("Testing permute() function")

test_that("permute() returns ith row of all perms", {
    h <- how()
    v <- 1:4
    ch <- check(v, h, quietly = TRUE)
    h <- ch$control
    p1 <- permute(1, n = length(v), control = h)
    expect_equal(length(p1), 4L)
    expect_is(p1, "integer")
    expect_equal(p1, c(1L, 2L, 4L, 3L))
    expect_equal(p1, getAllperms(h)[1, ])

    p21 <- permute(21, n = length(v), control = h)
    expect_equal(length(p21), 4L)
    expect_is(p21, "integer")
    expect_equal(p21, c(4L, 2L, 3L, 1L))
    expect_equal(p21, getAllperms(h)[21, ])
})

test_that("permute() returns a random permutation if no $allperms", {
    h <- how()
    v <- 1:10                           # want something big so no allperms
    p <- permute(10, n = length(v), control = h)
    expect_equal(length(p), length(v))
    expect_true(all(p >= 1))
    expect_true(all(p <= 10))

    setComplete(h) <- TRUE
    expect_warning(permute(10, n = length(v), control = h),
                   regexp = "Returning a random permutation")
})

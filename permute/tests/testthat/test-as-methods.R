library("testthat")
library("permute")

context("Testing as.foo() methods")

test_that("as.matrix allPerms method", {
    ap <- allPerms(1:3)
    m <- as.matrix(ap)
    expect_is(m, "matrix")
    expect_false(inherits(m, "allPerms"))
})

test_that("as.matrix permutationMatrix method", {
    perms <- shuffleSet(10, nset = 10)
    m <- as.matrix(perms)
    expect_is(m, "matrix")
    expect_false(inherits(m, "permutationMatrix"))
})

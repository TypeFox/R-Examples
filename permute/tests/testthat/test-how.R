library("testthat")
library("permute")

context("Testing how()")

test_that("how() works with explicit NULL blocks arg", {
    ## Example of failure from Jari github #8
    h <- how(blocks = NULL)
    expect_that(h, is_a("how"))
})

test_that("print method for how", {
    expect_output(print(how()), regexp = "Permutation Design:")
})

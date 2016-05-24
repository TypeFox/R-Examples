context("Testing rand names works")

test_that("Random names work", {
    expect_is(rand_names(1), 'data.frame')
})

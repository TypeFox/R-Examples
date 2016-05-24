context("UNFv6: Factors")

test_that("Factors treated as character", {
    expect_equal(unf6(c('1','2','3')), 
                 unf6(factor(c('1','2','3'))))
})

test_that("Factors can be treated as integer", {
    expect_equal(unf6(1:3), 
                 unf6(factor(c('a','b','c')), factor_as_character=FALSE))
    expect_equal(unf6(1:3), 
                 unf6(factor(c('1','2','3')), factor_as_character=FALSE))
})

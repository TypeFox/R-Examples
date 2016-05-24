context("Test of general functionality")

test_that("citation reports year", {
    txt <- citation("sensR")
    expect_true(as.logical(grep("year", txt)))
})


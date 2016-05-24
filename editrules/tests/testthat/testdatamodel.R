require(testthat)
context("datamodel")

test_that("datamodel works",{
    expect_identical(
        datamodel(editarray(c(
            "x %in% c('a','b')",
            "y %in% c('x','y')"
        ))),
        data.frame(variable = c('x','x','y','y'),value=c('a','b','x','y'))
   )
})







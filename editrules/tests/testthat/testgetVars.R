
require(testthat)

context("getVars")

test_that("getVars.editmatrix works",{
        expect_identical(
                getVars(editmatrix(c( "x+3*y == 2*z", "x > 2"))),
                c("x","y","z")
        )
})

test_that("getVars.editarray conforms to type argument",{
   expect_true(is.null(
      getVars(
         editarray(expression(if ( A == 'a') B %in% letters[1:3])),
         type='num'
      )
   ))
})




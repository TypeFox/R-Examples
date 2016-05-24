context("is.error")

test_that("is.error", {
 expect_true(is.error(try(stop("foo"), silent=TRUE)))
 expect_false(is.error(try("foo")))
})  

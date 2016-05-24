context("Accessor")

test_that("Class Accessor", {
  suppressWarnings({
    test <- defineClass("TestAcessor", contains = "Accessor", {
      
      x <- 1
      .y <- 2
      
      doSomething <- function() {
        .self$.y <- .self$.y + 1
        invisible(.self)
      }
      
      get <- function() {
        .self$.y
      }
      
    })
  })
  
  instance <- test()
  expect_is(instance$doSomething(), "TestAcessor")
  expect_equal(instance$get(), 3)
  expect_equal(instance$x, 1)
  expect_equal(instance$x <- 2, 2)
  expect_equal(instance$x, 2)
  
  removeClass("TestAcessor")
  
})


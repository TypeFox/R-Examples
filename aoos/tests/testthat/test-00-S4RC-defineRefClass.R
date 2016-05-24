context("defineRefClass")

test_that("Class wrapper", {
  
  suppressWarnings(
    Test <- defineRefClass({
      
      Class <- "Test"
      
      x <- "character"
      y <- "numeric"
      
      doSomething <- function() {
        .self$y <- y + 1
        invisible(.self)
      }
      
    }))
  
  instance <- new("Test", x = "Working", y = 0)
  
  expect_equal(instance$y, 0)
  expect_equal(instance$x, "Working")
  expect_is(instance$doSomething(), "Test")
  expect_equal(instance$y, 1)
  
  # Inheritance
  
  Character <- setClass("Character", contains = "character")
  
  suppressWarnings({
    SubTest <- defineRefClass({
      
      Class <- "SubTest"
      contains <- "Test"
      
      x <- "Character"
      xx <- "Test"
      
      doSomething <- function() {
        .self$y <- y + as.numeric(x)
        invisible(.self)
      }
      
    })
  })
    
  instance <- SubTest(x = Character("2"), y = 5)
  
  expect_equal(instance$y, 5)
  expect_equal(instance$x, Character("2"))
  expect_is(instance$doSomething(), "SubTest")
  expect_equal(instance$y, 7)
  expect_is(instance$xx, "Test")
  
  instance$xx$y <- 0
  instance$xx$x <- "Working"
  expect_equal(instance$xx$y, 0)
  expect_equal(instance$xx$x, "Working")
  expect_is(instance$xx$doSomething(), "Test")
  expect_equal(instance$xx$y, 1)
  
  removeClass("Test")
  removeClass("SubTest")
  
})

test_that("Private members for refClasses", {
  
  suppressWarnings({
    Test <- defineRefClass({
      Class <- "Test"
      contains <- "Private"
      
      .p <- "numeric"
      
      getP <- function() .self$.p
      setP <- function(v) .self$.p <- v
      
    })
  })
    
  test <- Test()
  # in testthat the correct methods can not be found. So the lines do not produce an error.
  # I will make an example in the doc which produces error.
#   print(showMethods("$"))
#   print(showMethods("[["))
#   expect_error(test$.p)
#   expect_error(test$.self)
  
  expect_equal(test$setP(2), 2)
  expect_equal(test$getP(), 2)
  expect_error(test[[".p"]])
  expect_error(test[[".p"]] <- 2)
  
  removeClass("Test")
  
})

test_that("refClass with empty fields", {
  
  suppressWarnings({
    Test <- defineRefClass({
      Class <- "Test3"  
    })
  })
  
  test <- Test()
  
  expect_true(inherits(test, "Test3"))
  
  removeClass("Test3")
  
})
context("Basics:")

test_that("class-markup", {
  
  suppressWarnings(
    test <- defineClass("test", {
      .privateObject <- NULL
      publicObject <- public("BAM!")
      get <- public(function() .privateObject)
      set <- public(function(value) .privateObject <<- value)
      doSomethingWithPublicObject <- public(function() publicObject())
    })
    )
  
  tmp <- test()
  expect_is(tmp$get(), "NULL")
  expect_equal(tmp$publicObject(), "BAM!")
  expect_error(tmp$.privateObject)
  expect_equal(tmp$set(2), 2)
  expect_equal(tmp$get(), 2)
  
  expect_equal(tmp$publicObject("jiggle"), "jiggle")
  expect_equal(tmp$doSomethingWithPublicObject(), "jiggle")
  
  # New instance with new environment:
  suppressWarnings(tmp2 <- test())
  expect_is(tmp2$get(), "NULL")
  expect_equal(tmp2$publicObject(), "BAM!")
  
})

test_that("default-call", {
  expect_is(suppressWarnings(defineClass("tmp")), "function")
  expect_true(isClass("tmp"))
})

test_that("privacy works", {
  class <- suppressWarnings({
    defineClass("class", {
      private <- private(2)
      get <- public(function() 1)
    })
  })
  
  tmp <- class()
  
  expect_equal(tmp$get(), 1)
  expect_error(tmp$private)
  expect_error(tmp$something <- 1)
  
})


test_that("'init' function is executed", {
  suppressWarnings({
    test <- defineClass("test", {
      
      name <- ""
      
      init <- function(name) {
        self$name(name)
      } 
      
    })
    
    defineClass("TestInit", {
      field <- ""
      init <- function() {
        self$field("b")
      }
    })
    
  })
  
  expect_equal(test("jaj")$name(), "jaj")
  expect_equal(new("test")$name(), "")
  expect_equal(new("TestInit")$field(), "b")
  
})

test_that("Naming of constructor functions", {
  constTest <- suppressWarnings({
    defineClass("test", {
      x <- public()
      })
  })
  
  tmp <- constTest()
  tmp1 <- new("test")
  
  tmp$x(2)
  
  expect_is(tmp, "test")
  expect_is(tmp1, "test")
  expect_is(tmp1$x(), "NULL")
  expect_equal(tmp$x(), 2)
  
})

test_that("class-markup v2.0", {
  
  suppressWarnings(
    test <- defineClass("test2", {
      .privateObject <- NULL
      publicObject <- "BAM!"
      get <- function() .privateObject      
      set <- function(value) .privateObject <<- value
      doSomethingWithPublicObject <- function() publicObject()
    })
  )
  
  tmp <- test()
  expect_is(tmp$get(), "NULL")
  expect_equal(tmp$publicObject(), "BAM!")
  expect_error(tmp$.privateObject)
  expect_equal(tmp$set(2), 2)
  expect_equal(tmp$get(), 2)
  
  expect_equal(tmp$publicObject("jiggle"), "jiggle")
  expect_equal(tmp$doSomethingWithPublicObject(), "jiggle")
  
  # New instance with new environment:
  suppressWarnings(tmp2 <- test())
  expect_is(tmp2$get(), "NULL")
  expect_equal(tmp2$publicObject(), "BAM!")
  
})

test_that("init function is private by default", {
  PrivateInitTest <- defineClass("PrivateInitTest", contains = "Accessor", {
    x <- 1
    init <- function() {
      .self$x <- 2
    }
  })
  
  expect_equal(PrivateInitTest()$x, 2)
  expect_error(PrivateInitTest()$init())
  
  # make init public explicitly
  PublicInitTest <- defineClass("PublicInitTest", contains = "Accessor", {
    x <- 1
    init <- public(function() {
      .self$x <- 2
    })
  })
  
  instance <- PublicInitTest()
  instance$x <- 3
  expect_equal(instance$x, 3)
  instance$init()
  expect_equal(instance$x, 2)
})

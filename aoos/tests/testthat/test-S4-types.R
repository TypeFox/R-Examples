context("S4Types")
test_that("Type with defaults", {
  
  Test(x = 1, y = list(z = 1, a = list(1, 3))) %type% {
    stopifnot(.Object@x > 0)
    .Object
  }
  
  expect_error(Test(x = 0))
  expect_true(Test()@x == 1)
  expect_true(Test(2)@x == 2)
  expect_true(identical(Test()@y$z, 1))
  expect_true(typeof(Test()) == "S4")
  
  removeClass("Test")
  
})

test_that("Type with ANY", {
  
  Test(x = numeric(), y = NULL) %type% .Object
  
  expect_true(is.null(Test()@y))
  expect_true(is.numeric(Test()@x))
  x <- Test()
  x@y <- Test()
  expect_true(is(x@y, "Test"))
  
  removeClass("Test")
  
})

test_that("Class without slot", {
  
  setClass("Empty", prototype = prototype(), where = environment())
  
  Empty : Test() %type% .Object
  expect_true(is(Test(), "Test"))
  
  removeClass("Test")
  
})

test_that("Type inheritance", {
  
  Test(x = 1, y = list()) %type% {
    stopifnot(.Object@x > 0)
    .Object
  }
  
  Test:Child(z = " ") %type% {
    stopifnot(nchar(.Object@z) > 0)
    .Object
  }
  
  expect_error(Child(x = 0))
  expect_true(Child()@x == 1)
  expect_true(identical(Child()@y, list()))
  expect_true(typeof(Test()) == "S4")
  expect_true(is(Child(), "Child"))
  expect_true(inherits(Child(), "Test"))
  expect_true(identical(Child(z = "char")@z, "char"))
  expect_error(Child(z = ""))
  expect_equal(Child(x = 5)@x, 5)
  
  removeClass("Child")
  removeClass("Test")
  
  # inheriting from more than one thing
  Test2(z = "") %type% .Object
  Test : Test2 : Child() %type% .Object
  
  expect_true(Child()@x == 1)
  expect_true(identical(Child()@y, list()))
  expect_true(inherits(Child(), "Test"))
  expect_true(inherits(Child(), "Test2"))
  expect_true(identical(Child(z = "char")@z, "char"))
  
})

test_that("Types can inherit S3 classes", {
  
  numeric : Test(x = 1, .Data = 2) %type% {
    .Object
  }
  
  expect_true(Test() == 2)
  expect_true(Test(4, 3) == 3)
  
  numeric : Test(x = 1, .Data = 2) %type% {
    .Object
  }
  
  expect_true(Test() == 2)
  expect_true(Test(4, 3) == 3)
  
  numeric : Test(x = 1) %type% .Object
  expect_equal(Test()@.Data, numeric())
  expect_true(Test(4, 3) == 3)
  
})

test_that("Type with VIRTUAL", {
  
  VIRTUAL:Type() %type% .Object
  
  names.Type <- function(x) {
    slotNames(x)
  }
  
  Type:Test(x = 1) %type% .Object  
  
  expect_true(names(Test(x = 2)) == "x")
  expect_true(inherits(Test(), "Type"))
  expect_error(new("Type"))
  
  removeClass("Test")
  removeClass("Type")
  
})

test_that("Type with quoted class names", {
  
  'numeric' : "Test"(names = character()) %type% .Object
  expect_true(inherits(Test(), "numeric"))
  expect_is(Test(), "Test")
  
})

test_that("Type with explicit class names", {
  
  Test(x ~ numeric, y = list(), z) %type% .Object
  
  expect_error(Test(x = 0)) # z is missing
  expect_true(Test(1, list(), NULL)@x == 1)
  expect_true(identical(Test(1, list(), NULL)@y, list()))
  expect_true(typeof(Test(1, list(), NULL)) == "S4")
  
  removeClass("Test")
  
})

test_that("Types can deal with class unions", {
  
  'numeric | character' : Test(
    x ~ 'numeric | character | list'
  ) %type% .Object
  
  expect_is(Test(1, 2)@x, "numeric")
  expect_is(Test("", "")@x, "character")
  expect_is(Test(list())@x, "list")
  
  removeClass("Test")
  removeClass("numericORcharacter")
  removeClass("numericORcharacterORlist")
  
})


context("Inheritance:")

test_that("Inheritance", {
  
  suppressWarnings({
    parent <- defineClass("parent", {
      publicMember <- publicValue("?!:")
      privateMember <- private(NULL)
      get <- publicFunction(function() paste(publicMember(), privateMember))
    })
    
    child <- defineClass("child", contains = "parent", {
      set <- publicFunction(function(value) {
        privateMember <<- value
        invisible(privateMember)
      })
    })
  })
  
  tmp <- child()
  expect_equal(tmp$publicMember(), "?!:")
  expect_error(tmp$privateMember) # don't know how to implement that
  expect_equal(tmp$get(), "?!: ")
  expect_equal(tmp$set("s"), "s")
  expect_equal(tmp$get(), "?!: s")
  
})

test_that("Replacing fields I", {
  
  suppressWarnings({
    parent <- defineClass("parent", {
      privateMember <- private(NULL)
      get <- publicFunction(function() privateMember)
    })
    
    child <- defineClass("child", contains = "parent", {
      privateMember <- private("value")
    })  
  })
    
  tmp <- child()
  expect_equal(tmp$get()@.Data, "value")
  
})

test_that("Replacing fields II", {
  
  suppressWarnings({
    parent <- defineClass("parent", {
      get <- function() foo()
      foo <- function() 1
    })
    
    child <- defineClass("child", contains = "parent", {
      foo <- function() 2
    })
  })
    
  tmp <- child()
  expect_equal(tmp$foo(), 2)
  expect_equal(tmp$get(), 2)
  
})

test_that("Inheritance of standard S4 classes", {
  
  setClass("Parent", contains = "VIRTUAL")
  
  setGeneric("testMethod", function(x) "default")
  setMethod("testMethod", signature = "Parent", function(x) 1)
  
  suppressWarnings({
    
    child1 <- defineClass("Child1", contains = "Parent", {
      get <- function() foo()
      foo <- function() 1
    })
    
    child2 <- defineClass("Child2", contains = "Child1", { })
    
  })
  
  tmp1 <- child1()
  tmp2 <- child2()
  
  expect_equal(testMethod(tmp1), 1)
  expect_equal(testMethod(tmp2), 1)
  
})

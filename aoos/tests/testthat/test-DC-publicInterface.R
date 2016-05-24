context("Public interface")
test_that("Validity function for set", {
  
  test <- suppressWarnings({
    test <- defineClass("test", {
      x <- public(5, function(x) inherits(x, what = "numeric"))
    })
  })
  
  tmp <- test()
  
  expect_equal(tmp$x(), 5)
  expect_is(tmp$x(), "numeric")
  expect_error(tmp$x("s"))
  
})

test_that("publicValue", {
  x <- publicValue()
  expect_is(x(), "NULL")
  expect_equal(x(1), 1)
  expect_equal(x(), 1)
  expect_is(public(function() 1), "publicFunction")
  expect_is(public(1), "publicValue")
})

test_that("Handling of reference classes as public member", {
  
  refObj <- suppressWarnings({
    defineClass("refObj", {
      method <- public(function() cat("...\n"))
      value <- public(1)
    })
  })
  
  refObj2 <- suppressWarnings({
    defineClass("refObj2", {
      refObj <- public(refObj())
    })
  })
  
  ro <- refObj2()
  
  # ro$refObj$method()
  ro$refObj$value()
  ro$refObj$value(2)
  ro$refObj$value()
  is(ro$refObj, "refObj")
  
})

test_that("Access methods reference fields", {
  
  suppressWarnings(
    test <- defineClass("test", {
      method <- public(function() {
        "..."
      })
    })
  )
  
  suppressWarnings(
    test1 <- defineClass("test1", {
      refObj <- public(test())
      method <- public(function() {
        refObj$method()
      })
    })
  )
  
  suppressWarnings(
    test2 <- defineClass("test2", {
      refObj <- public(test1())
      method <- public(function() {
        refObj$method()
      })
    })
  )
  
  expect_equal(test()$method(), test1()$method())
  
  instance <- test2()  
  expect_equal(instance$refObj$refObj$method(), test()$method())
  expect_equal(instance$refObj$method(), test()$method())
  expect_equal(instance$method(), test()$method())
  
})


test_that("Enforce privacy", {
  
  test <- suppressWarnings({
    test <- defineClass("test", {
      x <- private(5)
    })
  })
  
  tmp <- test()
  
  expect_error(tmp$x)
  
})

test_that("Enforce public", {
  
  test <- suppressWarnings({
    test <- defineClass("testEnforcePublic", {
      .x <- public(5)
    })
  })
  
  tmp <- test()
  
  expect_equal(tmp$.x(), 5)
  
})
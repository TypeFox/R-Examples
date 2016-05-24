context("argument check")

test_that("checkArg", {
  f = function(x) checkArg(x, cl="integer", len=1)
  f(1L)
  expect_error(f(1))
  expect_error(f(1:2))
  f = function(x) checkArg(x, cl="integer", min.len=2)
  f(1:2)
  expect_error(f(1L))
  f = function(x) checkArg(x, cl="integer", max.len=2)
  f(1:2)
  expect_error(f(1:3))
  f = function(x) checkArg(x, cl="integer", len=2, na.ok=FALSE)
  f(1:2)
  expect_error(f(c(3L, NA)))
  f = function(x) checkArg(x, cl="numeric", lower=2)
  f(4:6)
  expect_error(f(1:3))
  f = function(x) checkArg(x, cl="numeric", upper=2)
  f(1)
  expect_error(f(3), "less than or equal 2")
  f = function(x) checkArg(x, cl="numeric", lower=1, upper=2)
  f(1); f(1.5); f(2)
  expect_error(f(0), "greater than or equal 1")
  expect_error(f(3), "less than or equal 2")

  f = function(x) checkArg(x, cl="numeric", lower=1, upper=5)
  f(1:5)
  expect_error(f(0:5), "greater than or equal 1")
  expect_error(f(1:6), "less than or equal 5")

  f = function(x) checkArg(x, formals=c("foo", "bar"))
  f(function(foo, bar) 1)
  f(function(foo, bar, bla) 1)
  expect_error(f(1), "must be of class function not: numeric")
  expect_error(f(function(blubb) 1), "must have first formal args")
  expect_error(f(function(foo) 1), "must have first formal args")

  checkArg(1, "vector")
  checkArg(1L, "vector")
  checkArg(TRUE, "vector")
  checkArg("a", "vector")
  checkArg(list(), "vector")
  checkArg(list(1), "vector")
})


test_that("checkArg with choices", {
  f = function(x) checkArg(x, choices=c("a", "b"))
  f("a")
  f("b")
  expect_error(f(c("a", "b")), "must be")
  expect_error(f(1), "must be")
  expect_error(f(NULL), "must be")
  expect_error(f(NA))

  f = function(x) checkArg(x, choices=list(NULL, 1L, data.frame()))
  f(1L)
  f(NULL)
  f(data.frame())
  expect_error(f(1), "must be")
  expect_error(f(list(1)), "must be")
})


test_that("checkArg with subset", {
  f = function(x) checkArg(x, subset=c("a", "b"))
  f("a")
  f("b")
  f(c("a", "b"))
  f(character(0))

  expect_error(f(1), "must be")
  expect_error(f(NA), "must be")

  f = function(x) checkArg(x, subset=list(NULL, 1L, data.frame()))
  f(1L)
  f(NULL)
  f(data.frame())
  f(list(NULL, data.frame()))
  expect_error(f(1), "must be")
  expect_error(f(list(1)), "must be")
})

test_that("checkArg with missing arg", {
  f = function(x) checkArg(x, "numeric")
  expect_error(f(), "Argument x must not be missing!")
})

# FIXME no idea why this does not run in "CMD check"
if (interactive()) {
test_that("checkArg with classes / s3 and s4", {
  x = 1
  class(x) = c("foo2", "foo1")
  checkArg(x, "foo1")
  checkArg(x, "foo1", s4=FALSE)
  checkArg(x, "foo1", s4=TRUE)
  checkArg(x, "foo2")
  checkArg(x, "foo2", s4=FALSE)
  checkArg(x, "foo2", s4=TRUE)
  mys41 = setClass("mys41", representation(x="numeric"))
  mys42 = setClass("mys42", contains="mys41", representation(y="numeric"))
  obj1 = mys41(x=3)
  obj2 = mys42(x=3, y=4)
  checkArg(obj1, "mys41", s4=TRUE)
  checkArg(obj2, "mys41", s4=TRUE)
  checkArg(obj2, "mys42", s4=TRUE)
})
}

test_that("checkArg with multiple classes", {
  checkArg(1, c("numeric", "list"))
  checkArg(1, c("numeric", "foo"))
  checkArg(1L, c("integer", "list"))
  checkArg(1L, c("integer", "foo"))
  checkArg(1L, c("numeric", "list"))
  checkArg(1L, c("numeric", "foo"))
})


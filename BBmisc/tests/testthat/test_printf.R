context("print*f variants")

test_that("messagef", {
  expect_message(messagef("xxx%ixxx", 123), "xxx123xxx")
})

test_that("catf", {
  expect_output(catf("xxx%ixxx", 123), "xxx123xxx")
})

test_that("catf into file", {
  fn = tempfile()
  catf("xxx%ixxx", 123, file=fn)
  s = readLines(fn)
  expect_equal(s, "xxx123xxx")
  unlink(fn)
})



test_that("warningf", {
  expect_warning(warningf("xxx%ixxx", 123), "xxx123xxx")
  f = function() warningf("123")
  # "Warning: " not caught by gives_warning
  expect_warning(f(), "123")
})

test_that("stopf", {
  expect_error(stopf("xxx%ixxx", 123), "xxx123xxx")
  f = function() stopf("123")
  # because try is called in throws_error 
  # (and prints a bit differently of course!!!!)
  # we get an extra space before the :
  expect_error(f(), "123")
})



context("time")

test_that("time works", {
  skip_on_cran()

  aa <- GET("http://httpbin.org/get/") %>% time(verbose = FALSE)
  bb <- GET("https://api.github.com/") %>% time(verbose = FALSE)
  cc <- GET("http://google.com/") %>% time(verbose = FALSE)

  expect_is(aa, "http_time")
  expect_is(aa$times, "list")
  expect_is(aa$averages, "list")
  expect_is(aa$request, "response")
  expect_equal(aa$request$url, "http://httpbin.org/get/")

  expect_is(bb, "http_time")
  expect_is(bb$times, "list")
  expect_is(bb$averages, "list")
  expect_is(bb$request, "response")
  expect_equal(bb$request$url, "https://api.github.com/")

  expect_is(cc, "http_time")
  expect_is(cc$times, "list")
  expect_is(cc$averages, "list")
  expect_is(cc$request, "response")
  expect_equal(cc$request$url, "http://www.google.com/")
})

test_that("time - count parameter works", {
  skip_on_cran()

  n2 <- time(GET("http://httpbin.org/get/"), count = 2, verbose = FALSE)
  n7 <- time(GET("http://httpbin.org/get/"), count = 7, verbose = FALSE)

  expect_is(n2, "http_time")
  expect_equal(length(n2$times), 2)

  expect_is(n7, "http_time")
  expect_equal(length(n7$times), 7)
})

test_that("time - delay parameter works", {
  skip_on_cran()

  expect_gt(
    system.time(time(GET("http://httpbin.org/get/"), count = 3, verbose = FALSE, delay = 1))[[3]],
    system.time(time(GET("http://httpbin.org/get/"), count = 3, verbose = FALSE, delay = 0.5))[[3]]
  )
})

test_that("time - flood parameter works", {
  skip_on_cran()

  expect_gt(
    system.time(time(GET("http://httpbin.org/get/"), count = 3, verbose = FALSE))[[3]],
    system.time(time(GET("http://httpbin.org/get/"), count = 3, verbose = FALSE, flood = TRUE))[[3]]
  )
})

test_that("time - fails well", {
  skip_on_cran()

  expect_error(time(), "argument \".request\" is missing")
  expect_error(time(GET("http://httpbin.org/get/"), count = 1, verbose = FALSE),
               "count parameter must be greater than 1")
})


context("TCP")

test_that("We can ping localhost", {

  ## Chances are, there is nothing here
  pr <- ping_port("127.0.0.1", port = 4695, count = 1)
  expect_equal(pr, NA_real_)

  ## Start web server
  try(tools::startDynamicHelp(start = FALSE), silent = TRUE)
  suppressMessages(tools::startDynamicHelp())
  pr <- ping_port("127.0.0.1", port = tools:::httpdPort, count = 1)
  expect_true(is.double(pr))
  expect_true(length(pr) == 1)
  expect_true(pr < 1000)

  ## Shut down web server
  tools::startDynamicHelp(start = FALSE)
})

test_that("We can ping a remote host", {

  ## There is surely nothing here
  pr <- ping_port("igraph.org", port = 4695, count = 1)
  expect_equal(pr, NA_real_)

  ## There is surely something here
  pr <- ping_port("httpbin.org", count = 1)
  expect_true(is.double(pr))
  expect_true(length(pr) == 1)
  expect_true(pr < 1000)
})

test_that("We don't wait too long", {

  ## TODO

})

test_that("We don't wait for the resolver", {

  ## TODO

})

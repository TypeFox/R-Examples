context("ipad iterator")

test_that("ipad functions properly when no fill is necessary", {
  it <- iterators::iter(1:9)
  it_ipad <- ipad(it)
  expect_equal(as.list(islice(it_ipad, end=9)), as.list(1:9))
})

test_that("ipad functions properly with the default fill argument", {
  it <- iterators::iter(1:9)
  it_ipad <- ipad(it)
  expect_equal(as.list(islice(it_ipad, end=10)), c(as.list(1:9), NA))
})

test_that("ipad functions properly with a specified fill argument", {
  it <- iterators::iter(1:9)
  it_ipad <- ipad(it, fill=TRUE)
  expect_equal(as.list(islice(it_ipad, end=11)), c(as.list(1:9), as.list(rep(TRUE, 2))))
})

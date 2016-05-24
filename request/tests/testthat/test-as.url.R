context("as.url")

test_that("as.url works", {
  aa <- as.url(5000)
  bb <- as.url("5000")
  cc <- as.url(':9200')
  dd <- as.url('9200')
  ee <- as.url('9200/stuff')

  expect_is(aa, "url")
  expect_is(bb, "url")
  expect_is(cc, "url")
  expect_is(dd, "url")
  expect_is(ee, "url")

  expect_true(grepl("localhost", aa))
  expect_true(grepl("localhost", bb))
  expect_true(grepl("localhost", cc))
  expect_true(grepl("localhost", dd))
  expect_true(grepl("localhost", ee))

  expect_true(grepl("5000", aa))
  expect_true(grepl("5000", bb))
  expect_true(grepl("9200", cc))
  expect_true(grepl("9200", dd))
  expect_true(grepl("9200", ee))
})

test_that("as.url fails well", {
  expect_error(as.url(), "no applicable method")
})

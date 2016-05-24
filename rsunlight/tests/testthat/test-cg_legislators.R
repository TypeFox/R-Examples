context("cg_legislators")

test_that("cg_legislators returns the correct", {
  skip_on_cran()

  aa <- cg_legislators()
  bb <- cg_legislators(latitude = 35.778788, longitude = -78.787805)

  # classes
  expect_is(aa, "data.frame")
  expect_is(aa, "sunlight")

  expect_is(bb, "data.frame")
  expect_is(bb, "sunlight")
  expect_equal(bb$state[1], "NC")

  # dimensions
  expect_equal(NROW(aa), 20)
})

test_that("cg_legislators geographic search works", {
  skip_on_cran()

  aa <- cg_legislators(zip = 97206)
  bb <- cg_legislators(latitude = 45.509013, longitude = -122.660602)

  # classes
  expect_is(aa, "data.frame")
  expect_is(aa, "sunlight")

  expect_is(bb, "data.frame")
  expect_is(bb, "sunlight")

  # right names back
  expect_true(any(grepl("Wyden", aa$last_name)))
  expect_true(any(grepl("Wyden", bb$last_name)))
  expect_true(any(grepl("Merkley", aa$last_name)))
  expect_true(any(grepl("Merkley", bb$last_name)))
})

test_that("cg_legislators vectorizing works", {
  skip_on_cran()

  res <- cg_legislators(last_name = c('Pelosi', 'Merkley'))
  expect_is(res, "data.frame")
  expect_is(res, "sunlight")
  expect_equal(NROW(res), 2)
  expect_equal(res$last_name, c('Pelosi', 'Merkley'))
})

test_that("cg_legislators pagination works", {
  skip_on_cran()

  expect_equal(NROW(cg_legislators(per_page = 1)), 1)
  expect_equal(NROW(cg_legislators(per_page = 10)), 10)
  expect_equal(NROW(cg_legislators(per_page = 1, page = 2)), 1)
})

test_that("cg_legislators curl options work", {
  skip_on_cran()

  library("httr")
  expect_error(cg_legislators(config = httr::timeout(0.001)), "Timeout")
})

test_that("cg_legislators fails well", {
  skip_on_cran()

  expect_error(cg_legislators(latitude = 35.778788, longitude = -78.787805, first_name = "Afd"),
               "If latitude, longitude, or zip are used, all other parameters must be NULL")
})

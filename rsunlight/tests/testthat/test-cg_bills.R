context("cg_bills")

test_that("cg_bills returns the correct", {
  skip_on_cran()

  aa <- cg_bills()
  bb <- cg_bills(sponsor.party = 'R', history.vetoed = TRUE)

  # classes
  expect_is(aa, "data.frame")
  expect_is(aa, "sunlight")

  expect_is(bb, "data.frame")
  expect_is(bb, "sunlight")
  expect_true(bb$history.vetoed[1])

  # dimensions
  expect_equal(NROW(aa), 20)
})

test_that("cg_bills vectorizing works", {
  skip_on_cran()

  bills <- c("hjres131-113", "hjres129-113", "s2921-113")
  res <- cg_bills(bill_id = bills)
  expect_is(res, "data.frame")
  expect_is(res, "sunlight")
  expect_equal(NROW(res), 3)
  expect_equal(res$bill_id, bills)
})

test_that("cg_bills pagination works", {
  skip_on_cran()

  expect_equal(NROW(cg_bills(per_page = 1)), 1)
  expect_equal(NROW(cg_bills(per_page = 10)), 10)
  expect_equal(NROW(cg_bills(per_page = 1, page = 2)), 1)
})

test_that("cg_bills curl options work", {
  skip_on_cran()

  library("httr")
  expect_error(cg_bills(config = httr::timeout(0.001)), "Timeout")
})

test_that("cg_bills fails well", {
  skip_on_cran()

  expect_error(cg_bills(history.active = "asdf"), "history.active must be logical")
  expect_error(cg_bills(history.vetoed = 5555), "history.vetoed must be logical")

  expect_equal(NROW(cg_bills(sponsor.party = 5)), 0)
})

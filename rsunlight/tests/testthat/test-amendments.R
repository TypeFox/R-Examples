context("cg_amendments")

test_that("returns the correct", {
  skip_on_cran()

  a <- cg_amendments()
  b <- cg_amendments(chamber='house', congress=113)
  d <- cg_amendments(sponsor_type='committee', sponsor_id='HSRU')
  e <- cg_amendments(amends_bill_id='hr624-113')

  # classes
  expect_is(a, "sunlight")
  expect_is(a, "data.frame")
  expect_is(b, "sunlight")
  expect_is(d, "sunlight")
  expect_is(e, "sunlight")

  # values
  expect_equal(unique(b$chamber), "house")
  expect_equal(unique(b$congress), 113)
  expect_equal(unique(d$sponsor_type), "committee")
  expect_equal(unique(d$sponsor_id), "HSRU")
  expect_equal(unique(e$amends_bill_id), 'hr624-113')
})

test_that("vectorizing works", {
  skip_on_cran()

  res <- cg_amendments(chamber = c('house', 'senate'))
  expect_is(res, "sunlight")
  expect_equal(unique(res$chamber), c('house', 'senate'))
})

test_that("paging works", {
  skip_on_cran()

  res <- cg_amendments(chamber = 'house', per_page = 2)
  expect_equal(NROW(res), 2)
})

test_that("fails well", {
  # Sunlight APIs don't error when given the wrong type of data for a
  # parameter, so returning nothing, or dropping silently is failing
  # well - or as good as we can
  skip_on_cran()

  expect_equal(NROW(cg_amendments(congress = 23434)), 0)
})

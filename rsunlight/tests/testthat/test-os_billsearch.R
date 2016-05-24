context("os_billsearch")

test_that("returns the correct", {
  skip_on_cran()

  a <- os_billsearch(terms = 'agriculture', state = 'tx')
  b <- os_billsearch(terms = 'taxi', state = 'dc', sort='created_at')
  d <- os_billsearch(terms = 'taxi', state = 'dc', type='resolution')

  # classes
  expect_is(a, "data.frame")
  expect_is(b, "data.frame")
  expect_is(d, "data.frame")
  expect_is(a$state, "character")

  # values
  expect_equal(unique(a$state), "tx")
  expect_equal(unique(b$state), "dc")
  dates <- as.Date(b$created_at)
  expect_more_than(max(dates), min(dates))
})

test_that("paging works", {
  skip_on_cran()

  page1 <- os_billsearch(terms = 'agriculture', state = 'tx', per_page=2)
  page2 <- os_billsearch(terms = 'agriculture', state = 'tx', per_page=4, page=4)

  expect_is(page1, "data.frame")
  expect_is(page2, "data.frame")

  expect_equal(NROW(page1), 2)
  expect_equal(NROW(page2), 4)
})

test_that("vectorizing works", {
  skip_on_cran()

  vec1 <- os_billsearch(terms = 'agriculture', state = c('tx', 'or'))

  expect_is(vec1, "data.frame")
  expect_equal(unique(vec1$state), c('tx', 'or'))
})

test_that("fails well", {
  skip_on_cran()

  expect_error(os_billsearch(), "Bad Request")
  expect_equal(NROW(os_billsearch(state = 3)), 0)
})

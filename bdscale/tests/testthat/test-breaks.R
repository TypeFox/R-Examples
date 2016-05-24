library(bdscale)

data(nyse)

context("breaks")

test_that('week breaks', {
  range <- as.Date(c('2015-01-02', '2015-01-30'))
  
  f <- bd_breaks(nyse)()
  
  mondays <- f(range)
  
  expect_equal(mondays, as.Date(c('2015-01-05', '2015-01-12', '2015-01-20', '2015-01-26')))
})

context("groups") 

test_that('months', {
  dates <- as.Date('2015-01-01') + 0:100
  firsts <- firstInGroup(dates, month_format)
  expect_equal(unname(firsts), as.Date(c('2015-01-01', '2015-02-01', '2015-03-01', '2015-04-01')))
})

test_that('quarters', {
  dates <- as.Date('2015-01-01') + 0:100
  firsts <- firstInGroup(dates, quarter_format)
  expect_equal(unname(firsts), as.Date(c('2015-01-01', '2015-04-01')))
})
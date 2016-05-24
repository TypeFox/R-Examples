library(bdscale)

data(nyse)

context("forwards")

test_that('holidays show up as NA', {
  holidays <- bd2t(as.Date(c('2015-07-04', '2015-01-19')), nyse)
  expect_equal(holidays, as.numeric(c(NA,NA)))
})

test_that('weekends show up as NA', {
  weekends <- bd2t(as.Date(c('2015-06-06', '2015-04-05')), nyse)
  expect_equal(weekends, as.numeric(c(NA,NA)))  
})

test_that('first and last dates', {
  limits <- bd2t(as.Date(c('1950-01-03', '2016-03-15')), nyse)
  expect_equal(limits, c(0, length(nyse) - 1))    
})

context('backwards')

test_that('off the ends', {
  limits <- t2bd(c(-2, length(nyse) + 20), nyse)
  expect_equal(limits, c(nyse[1], nyse[length(nyse)]))
})

test_that('samples', {
  sample <- t2bd(c(100, 200, 300), nyse)
  expect_equal(sample, as.Date(c('1950-05-26', '1950-10-19', '1951-03-16')))
})
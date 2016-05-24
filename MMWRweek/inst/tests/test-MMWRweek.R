context('MMWRweek')

dates = seq(as.Date('1900-01-01'), as.Date('2015-12-31'), by='day')
truth = readRDS('MMWRweek-truth.rds')

test_that('MMWRweek returns a factor vector', {
  week = MMWRweek(dates)
  expect_that(week, is_a('data.frame'))
  expect_that(names(week), equals(c('MMWRyear', 'MMWRweek', 'MMWRday')))
})

test_that('MMWRweek returns correct values', {
  dates_chr = as.Date(dates)
  expect_that(MMWRweek(dates    ), equals(truth))
  expect_that(MMWRweek(dates_chr), equals(truth))
})
 

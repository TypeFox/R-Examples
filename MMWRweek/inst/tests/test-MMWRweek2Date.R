context('MMWRweek2Date')

test_that('MMWRweek2Date returns errors', {
  expect_error(MMWRweek2Date())
  expect_error(MMWRweek2Date(2015))
  expect_error(MMWRweek2Date('a',1))
  expect_error(MMWRweek2Date(2015,'a'))
  expect_error(MMWRweek2Date(2015,0))
  expect_error(MMWRweek2Date(2015,54))
})

n = 20
d = data.frame(
  MMWRyear = sample(1900:2100, n, replace=TRUE),
  MMWRweek = sample(1:52,      n, replace=TRUE),
  MMWRday  = sample(1:7,       n, replace=TRUE)
)

test_that('MMWRweek2Date matches MMWRweek', {
  date = with(d, MMWRweek2Date(MMWRyear,MMWRweek,MMWRday))
  expect_equal(d, MMWRweek(date))
})


test_that('MMWRweek2Date defaults to return first day of the MMWR week', {
  expect_equal(MMWRweek2Date(d$MMWRyear, d$MMWRweek),
               MMWRweek2Date(d$MMWRyear, d$MMWRweek, rep(1,length(d$MMWRyear))))
})
 


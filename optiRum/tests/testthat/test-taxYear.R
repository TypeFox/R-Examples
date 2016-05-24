context("taxYear")

test_that("taxYear correct results, uk default", {
  res<-taxYear(seq(from=as.Date("2015-04-05"),by=1,length=3))
  check<-c(2014,2015,2015)
  expect_equal(res,check)
})

test_that("taxYear correct results, jan 1", {
  res<-taxYear(seq(from=as.Date("2015-01-01"),by=1,length=2),"01-01")
  check<-c(2015,2015)
  expect_equal(res,check)
})


test_that("taxYear feb 29 aka mar 1st!", {
  expect_error(taxYear(seq(from=as.Date("2015-02-27"),by=1,length=3),"02-29"))
  
  res<-taxYear(seq(from=as.Date("2016-02-27"),by=1,length=3),"02-29")
  check<-c(2015,2015,2016)
  expect_equal(res,check)
})
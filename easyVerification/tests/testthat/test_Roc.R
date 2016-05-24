library(easyVerification)
context('ROC scores')

obs <- diag(1, 4)[-4,]
fcst <- obs
fcst2 <- obs*0.4 + 0.2
fcst3 <- array(0.25, dim(obs))

test_that('ROC area values', {
  expect_is(EnsRoca(fcst, obs), 'list')
  expect_equal(length(EnsRoca(fcst, obs)), 4)
  expect_equal(EnsRoca(fcst, obs), list(cat1=1, cat2=1, cat3=1, cat4=as.numeric(NA)))
  expect_equal(EnsRoca(fcst3, obs), list(cat1=0.5, cat2=0.5, cat3=0.5, cat4=as.numeric(NA)))
  expect_equal(EnsRoca(fcst, obs), EnsRoca(fcst2, obs))
})
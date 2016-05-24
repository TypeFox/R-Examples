context('Check scale functions')

load('./tabari-demo.RData')
gold <- make_scale('goldstein', file='./demo-scale.csv')

test_that("make_eventscale", {
  expect_is(gold, 'eventscale')
  expect_equal(length(gold), 109)
  expect_equal(score(gold, '062'), 2.5)
})

test_that("score", {
  scs <- score(gold, tabari.demo$code)
  s25to36 <- c(1.9,2.8,3.4,3.4,3.4,NA,NA,NA,NA,3.4,3.6,4.0)
  expect_true(all.equal(scs[25:36], s25to36)) ## deal with NAs
})

test_that("add_eventscale", {
  dd <- add_eventscale(tabari.demo, gold)
  expect_true('goldstein' %in% names(dd)) 

  gold2 <- make_scale('goldstein2', file='./demo-scale.csv', default='-99')
  dd <- add_eventscale(dd, gold2)
  expect_true('goldstein2' %in% names(dd))
  expect_true(is.na(dd$goldstein[58]))
  expect_true(dd$goldstein2[58]==-99)

})




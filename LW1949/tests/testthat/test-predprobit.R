
test_that("predprobit() throws errors when it should", {

  toxdat <- data.frame(
    dose=c(0.05, 0.0625, 0.125, 0.25, 0.5, 1),
    ntot=rep(8, 6),
    nfx = c(0, 1, 4, 4, 6, 8))
  myfit <- fitprobit(toxdat)

  expect_that(predprobit(pct="A", myfit), throws_error())
  expect_that(predprobit(pct=NA, myfit), throws_error())
  expect_that(predprobit(pct=c(50, 75), myfit), throws_error())
  expect_that(predprobit(pct=-5, myfit), throws_error())
  expect_that(predprobit(pct=0, myfit), throws_error())
  expect_that(predprobit(pct=100, myfit), throws_error())
  expect_that(predprobit(pct=200, myfit), throws_error())

  expect_that(predprobit(50, myfit, alpha="A"), throws_error())
  expect_that(predprobit(50, myfit, alpha=NA), throws_error())
  expect_that(predprobit(50, myfit, alpha=c(0.5, 0.6)), throws_error())
  expect_that(predprobit(50, myfit, alpha=-5), throws_error())
  expect_that(predprobit(50, myfit, alpha=0), throws_error())
  expect_that(predprobit(50, myfit, alpha=100), throws_error())
  expect_that(predprobit(50, myfit, alpha=200), throws_error())

  expect_that(predprobit(50, myfit, logbase="A"), throws_error())
  expect_that(predprobit(50, myfit, logbase=TRUE), throws_error())
  expect_that(predprobit(50, myfit, logbase=NA), throws_error())
  expect_that(predprobit(50, myfit, logbase=c(50, 75)), throws_error())
  expect_that(predprobit(50, myfit, logbase=-5), throws_error())
  expect_that(predprobit(50, myfit, logbase=1), throws_error())
  expect_that(predprobit(50, myfit, logbase=200), throws_error())

})

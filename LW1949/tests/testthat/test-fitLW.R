
test_that("fitLWauto() throws errors when it should", {

  dose <- c(0.0625, 0.125, 0.25, 0.5, 1)
  ntested <- rep(8, 5)
  nalive <- c(1, 4, 4, 7, 8)
  mydat <- dataprep(dose=dose, ntot=ntested, nfx=nalive)

  mymat <- as.matrix(mydat)

  mydat1 <- mydat[, -1]
  mydat2 <- mydat[, -2]
  mydat3 <- mydat[, -3]
  mydat5 <- mydat[, -5]
  mydat6 <- mydat[, -6]
  mydat7 <- mydat[, -7]
  mydat8 <- mydat[, -8]
  mydat9 <- mydat[, -9]

  mydatdup <- rbind(mydat, mydat[1, ])

  expect_that(fitLWauto(mymat), throws_error())

  expect_that(fitLWauto(mydat1), throws_error())
  expect_that(fitLWauto(mydat2), throws_error())
  expect_that(fitLWauto(mydat3), throws_error())
  expect_that(fitLWauto(mydat5), throws_error())
  expect_that(fitLWauto(mydat6), throws_error())
  expect_that(fitLWauto(mydat7), throws_error())
  expect_that(fitLWauto(mydat8), throws_error())
  expect_that(fitLWauto(mydat9), throws_error())

  expect_that(fitLWauto(mydatdup), throws_error())

})

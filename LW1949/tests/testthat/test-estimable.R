
test_that("estimable() throws errors when it should", {

  conc <- c(0.0625, 0.125, 0.25, 0.5, 1)
  pfx <- c(1, 4, 4, 7, 8)/8

  mydat <- data.frame(dose=conc, pfx=pfx)

  mymat <- as.matrix(mydat)

  mydat1 <- mydat
  names(mydat1)[1] <- "bob"
  mydat2 <- mydat
  names(mydat2)[2] <- "dorothy"

  conc.dup <- conc[c(1, 1, 3:5)]
  mydatdup <- data.frame(dose=conc.dup, pfx=pfx)

  expect_that(estimable(mymat), throws_error())
  expect_that(estimable(mydat1), throws_error())
  expect_that(estimable(mydat2), throws_error())
  expect_that(estimable(mydatdup), throws_error())

})

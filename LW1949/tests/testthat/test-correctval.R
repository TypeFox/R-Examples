test_that("correctval() throws errors when it should", {

  gamfit <- gamtable1()

  val1 <- c(0.37, 0.5, 0.63)
  val2 <- 0.37
  val3 <- c(0.37, NA, 0.63)
  val4 <- as.numeric(rep(NA, 3))
  val5 <- as.numeric(NA)

  val6 <- letters[1:3]
  val7 <- rep(NA, 3)
  val8 <- c(-0.2, 0.5, 0.63)
  val9 <- c(0, 0.5, 0.63)
  val10 <- c(0.37, 0.5, 1)
  val11 <- c(0.37, 0.5, 1.2)

  expect_that(is.numeric(correctval(val1, gamfit)), is_true())
  expect_that(is.numeric(correctval(val2, gamfit)), is_true())
  expect_that(is.numeric(correctval(val3, gamfit)), is_true())
  expect_that(all(is.na(correctval(val4, gamfit))), is_true())
  expect_that(all(is.na(correctval(val5, gamfit))), is_true())
  expect_that(correctval(val6, gamfit), throws_error())
  expect_that(correctval(val7, gamfit), throws_error())
  expect_that(correctval(val8, gamfit), throws_error())
  expect_that(correctval(val9, gamfit), throws_error())
  expect_that(correctval(val10, gamfit), throws_error())
  expect_that(correctval(val11, gamfit), throws_error())
})

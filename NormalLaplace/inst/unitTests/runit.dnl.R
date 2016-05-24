# Testing dnl
test.dnl <- function() {
  param <- c(0, 1, 2, 3)
  tol <- 0.001

  testdnl1 <- dnl(1, param = param)
  pdnl1 <- 0.259102

  testdnl2 <- dnl(2, param = param)
  pdnl2 <- 0.093693

  testdnl3 <- dnl(3, param = param)
  pdnl3 <- 0.019365

  checkEquals(testdnl1, pdnl1, tol = tol)
  checkEquals(testdnl2, pdnl2, tol = tol)
  checkEquals(testdnl3, pdnl3, tol = tol)
}


# Testing pnl
test.pnl <- function() {
  param <- c(0, 1, 2, 3)
  tol <- 0.001

  testpnl1 <- pnl(0, param = param)
  ppnl1 <- 0.44774

  testpnl2 <- pnl(1, param = param)
  ppnl2 <- 0.76905

  testpnl3 <- pnl(2, param = param)
  ppnl3 <- 0.94081

  checkEquals(testpnl1, ppnl1, tol = tol)
  checkEquals(testpnl2, ppnl2, tol = tol)
  checkEquals(testpnl3, ppnl3, tol = tol)

  # Checking that the CDF is equal to the inverse
  # of the inverse of the CDF
  checkTrue(pnl(1) == pnl(qnl(pnl(1))))
}


# Testing qnl
test.qnl <- function() {
  param <- c(0, 1, 2, 3)
  tol <- 0.001

  testqnl1 <- qnl(0.25, param = param)
  pqnl1 <- -0.61571

  testqnl2 <- qnl(0.5, param = param)
  pqnl2 <- 0.14971

  testqnl3 <- qnl(0.75, param = param)
  pqnl3 <- 0.92817

  checkEquals(testqnl1, pqnl1, tol = tol)
  checkEquals(testqnl2, pqnl2, tol = tol)
  checkEquals(testqnl3, pqnl3, tol = tol)

  # Checking that the CDF is equal to the inverse
  # of the inverse of the CDF
  checkTrue(pnl(1) == pnl(qnl(pnl(1))))
}


# Testing rnl
test.rnl <- function() {
  param <- c(0, 1, 2, 3)
  names(param) <- c("mu", "sigma", "alpha", "beta")

  # RUnit uses kind = "M-M", normal.kind = "K-R" for RNG. See ?RNGkind
  set.seed(2242, kind = "Marsaglia-Multicarry")
  dataVector <- rnl(5, param = param)

  checkEquals(dataVector, c(0.8048718, 1.6375198, -0.4043393, 0.9219744, -0.6499051))
}

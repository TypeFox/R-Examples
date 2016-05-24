# Testing millsR
test.millsR <- function() {
  ## check some ordinary values
  xVals <- 1:5
  yVals <- millsR(xVals)
  yReqd <- c(0.65568, 0.42137, 0.30459, 0.23665, 0.19281)
  tol <- 0.001
  checkEquals(yReqd, yVals, tol = tol)

  ## check some large values
  tol <- 10^(-14)
  xVals <- (1:10)*10^7
  yVals <- millsR(xVals)

  yReqd <-
    c(1.00090856345226e-07, 4.99371810711756e-08,
      3.43212891632624e-08, 2.51099915574398e-08, 1.95556810878505e-08,
      1.95556810878505e-08, 1.52299797447126e-08, 1.52299797447126e-08,
      1.52299797447126e-08, 1.52299797447126e-08)
  ## Use percentage difference
  checkEquals((yVals - yReqd)/yReqd, rep(0, length(yVals)), tol = tol)

  # Checking that e^log(millsR) = millsR.
  checkEquals(millsR(1:5, log = TRUE), log(millsR(1:5)), tol = tol)
}

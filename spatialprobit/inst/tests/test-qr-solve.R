context("qr-solve")

################################################################################
#
# QR decomposition for sparse Matrix: 
# test qr.solve(S, b) vs. solve(qr(S), b)
# mu <-qr.solve(S,X %*% beta) and 
# mu <-solve(qr(S),X %*% beta)
#
################################################################################

library(spatialprobit)
set.seed(2)
n <- 200
I_n <- sparseMatrix(i=1:n, j=1:n, x=1)
rho <- 0.75
W <- kNearestNeighbors(x=rnorm(n), y=rnorm(n))
rho <- 0.5
S <- (I_n - rho * W)
QR <- qr(S)

test_that("solve(qr(S),b) is faster than qr.solve(S,b) for sparse S", {  
  b <- rep(1, n)
  
  time1 <- system.time(for (i in 1:500) y1 <- as.double(solve(QR, b)))
  # user  system elapsed 
  # 1.83    0.02    1.84
  time2 <- system.time(for (i in 1:500) y2 <- as.double(qr.solve(S, b)))
  # user  system elapsed
  # 5.36    0.14    5.52
  
  expect_that(y1, equals(y2))  # expect same results
  #expect_that(time1["elapsed"] < time2["elapsed"] * 0.7, is_true())  # expect at least 30% performance gain
  cond <- time1["elapsed"] < (time2["elapsed"] * 0.7)
  names(cond) <- NULL                                  # problem with named condition
  expect_true(cond)
})
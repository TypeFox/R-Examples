## Daniel Gerlanc and Kris Kirby (2010-2012)
## Input and regression tests for the scaleLambdasBySide function

test.scaleLambdasBySide <- function() {
  ## Test the functioning of the scaleLambdasBySide function

  ## Simplest Case and its permutation
  x1   = c(-1, 1)
  res1 = bootES:::scaleLambdasBySide(x1)
  checkEquals(x1, res1)
  
  x1a   = c(1, -1)
  res1a = bootES:::scaleLambdasBySide(x1a)
  checkEquals(x1a, res1a)

  ## Simplest Case w/ names and its permutation
  x1na  = c(A=-1, B=1)
  res1na = bootES:::scaleLambdasBySide(x1na)
  checkEquals(x1na, res1na)
  
  x1nb   = c(B=-1, A=1)
  res1nb = bootES:::scaleLambdasBySide(x1nb)
  checkEquals(x1nb, res1nb)
  
  x2   = c(-0.33, 0.33)
  res2 = bootES:::scaleLambdasBySide(x2)
  checkEquals(res2, c(-1, 1))

  x3   = c(-0.25, -0.5, 0.5)
  res3 = bootES:::scaleLambdasBySide(x3)
  truth3 = c(-1/3, -2/3, 1)
  checkEquals(res3, truth3, tol=1e-8)
  
  x4     = c(-5, 0, 6)
  res4   = bootES:::scaleLambdasBySide(x4)
  truth4 = c(-1, 0, 1)
  checkEquals(res4, truth4)

  x5     = c(0, 0)
  res5   = bootES:::scaleLambdasBySide(x5)
  checkEquals(res5, x5)
  
}

test.calcSlopeLambdas <- function() {
  ## Regression test for calcSlopeLambdas

  ## Simplest Case
  x1    <- c(A=30, B=60, C=120)
  truth <- c(A=-0.0095, B=-0.0023, C=0.0119)
  res   <- bootES:::calcSlopeLambdas(x1)
  checkEquals(res, truth, tolerance=1e-2)
}

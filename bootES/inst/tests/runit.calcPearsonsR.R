## Daniel Gerlanc and Kris Kirby (2010-2012)
## Input and regression tests for the calcPearsonsR function

test.calcPearsonsR <- function() {
  ## Test the functioning of calcPearsonsR function

  ## Create 2-group data.frame.
  g1        = c(11, 12, 13, 14, 15)
  g2        = c(26, 27, 28, 29)
  g3        = c(17, 18, 19, 20, 21, 22, 23)
  grpLabels = rep(c("A", "B"), times=c(length(g1), length(g2)))
  twoGpsA   = data.frame(x=c(g1, g2), team=grpLabels)
  twoGpsErr = data.frame(x=c(g1, g2), team=rep("A", length(c(g1, g2))))
  lambdas   = c(A=1, B=-1)

  ## Create 3-group data.frame.
  grpLabels3  = rep(c("A", "B", "C"), times=c(length(g1), length(g2), length(g3)))
  threeGps    = data.frame(grpLabels3, scores=c(g1, g2, g3))
  threeGpsVec = c(g1, g2, g3)
  lambdas3    = c(A=-1, B=2, C=-1)
  
  ## Regression test of calcPearsonsR for a 3-group contrast
  set.seed(1)
  truth = 0.9452
  r.res = bootES:::calcPearsonsR(threeGpsVec, freq=rep(1, length(threeGpsVec)), grps=grpLabels3, lambdas=lambdas3)
  checkEquals(truth, r.res, tol=1e-4)
}

## Daniel Gerlanc and Kris Kirby (2010-2012)
## Input and regression tests for the calcCohensD function

test.calcCohensD <- function() {
  ## Test the functioning of calcCohensD function

  ## Create 2-group data.frame.
  g1        = c(11, 12, 13, 14, 15)
  g2        = c(26, 27, 28, 29)
  twoGpsVec = c(g1, g2)
  grpLabels = rep(c("A", "B"), times=c(length(g1), length(g2)))
  freqVec   = rep(1, length(twoGpsVec))
  lambdas   = c(A=1, B=-1)
  
  ## Regression test of calcCohensD
  truth = (mean(g1) - mean(g2)) / bootES:::pooledSD(twoGpsVec, grpLabels)
  d.res = bootES:::calcCohensD(twoGpsVec, freq=freqVec, grps=grpLabels,
    contrast=lambdas)  
  
  checkEquals(truth, d.res, tol=1e-4)
  
  d.res.switched = bootES:::calcCohensD(twoGpsVec, freq=freqVec,
    grps=grpLabels, contrast=c(A=-1, B=1))
  checkEquals(-1 * truth, d.res.switched, tol=1e-4)
  
  ## Regression test of calcCohensD w/ cohens.d.sigma=TRUE
  truth = (mean(g1) - mean(g2)) / bootES:::pooledSD(twoGpsVec, grpLabels,
                                                    pop.sd=TRUE)
  d.res = bootES:::calcCohensD(twoGpsVec, freq=freqVec, grps=grpLabels,
    contrast=lambdas, cohens.d.sigma=TRUE)
  checkEquals(truth, d.res, tol=1e-4)
  
  ## Regression test of calcCohensD w/ glass.control=TRUE
  truth = (mean(g1) - mean(g2)) / sqrt(sum(((g1 - mean(g1))^2) / length(g1)))
  d.res.glass = bootES:::calcCohensD(twoGpsVec, freq=freqVec,
    grps=grpLabels, contrast=lambdas, cohens.d.sigma=TRUE, glass.control="A")

  checkEquals(truth, d.res.glass, tol=1e-4)
  
  ## Regression test of calcCohensD w/ glass.control=TRUE and
  ## cohens.d.sigma=FALSE
  truth = (mean(g1) - mean(g2)) / sd(g1)
  d.res.glass = bootES:::calcCohensD(twoGpsVec, freq=freqVec,
    grps=grpLabels, contrast=lambdas, cohens.d.sigma=FALSE, glass.control="A")
  checkEquals(truth, d.res.glass, tol=1e-4)
}

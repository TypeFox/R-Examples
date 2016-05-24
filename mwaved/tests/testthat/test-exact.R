context("exact reconstruction tests")

n <- 256
shape <- 0.5
scale <- 0.25
signal <- makeLIDAR(n)
G <- gammaBlur(n, shape = 0.5, scale = 0.25)
X <- blurSignal(signal, G)
# Check multiFunctions consistent
mwvd <- multiWaveD(X, G)
multiEst <- multiEstimate(X,thr=0)
j1 <- mwvd$j1
j0 <- mwvd$j0
coef <- multiCoef(X, G, j1 = j1)
sig <- multiSigma(X)
thresh <- multiThresh(X, G)
scoef <- waveletThresh(coef, thresh)
teta <- theoreticalEta(alpha = 1, blur = detectBlur(G), G = G, sigma = sig)
test_that("Correct objects in mWaveD-list", {
  expect_equal(mwvd$estimate, multiEstimate(X, G))
  expect_equal(mwvd$coef, coef)
  expect_equal(mwvd$sigmaEst, sig)
  expect_equal(mwvd$thresholds, thresh[1:(j1 - j0 +1)])
  expect_equal(mwvd$shrinkCoef, scoef)
  expect_equal(mwvd$eta, teta)
})

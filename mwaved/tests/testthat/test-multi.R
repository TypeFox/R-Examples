context("multiFunction tests")

n <- 256
m <- 3
shape <- seq(from = 0.25, to = 0.75, length = m)
scale <- rep(0.25, m)
signal <- makeLIDAR(n)
G <- gammaBlur(n, shape = c(0.25, 0.5, 0.75), scale = rep(0.25, m))
X <- blurSignal(signal, G) 
E <- multiNoise(n, sigma = sigmaSNR(X, SNR = rep(20, m)))
Y <- X + E
# Check multiFunctions consistent
mwvd <- multiWaveD(Y, G)
j1 <- mwvd$j1
j0 <- mwvd$j0
coef <- multiCoef(Y, G, j1 = j1)
sig <- multiSigma(Y)
thresh <- multiThresh(Y, G)
scoef <- waveletThresh(coef, thresh)
teta <- theoreticalEta(alpha = rep(1, m), blur = detectBlur(G), G = G, sigma = sig)
test_that("Correct objects in mWaveD-list", {
  expect_equal(mwvd$estimate, multiEstimate(Y, G))
  expect_equal(mwvd$coef, coef)
  expect_equal(mwvd$sigmaEst, sig)
  expect_equal(mwvd$thresholds, thresh[1:(j1 - j0 +1)])
  expect_equal(mwvd$shrinkCoef, scoef)
  expect_equal(mwvd$eta, teta)
})

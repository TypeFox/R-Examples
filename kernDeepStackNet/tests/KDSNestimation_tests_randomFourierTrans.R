# Test randomFourierTrans
library(kernDeepStackNet)
# Generate covariate matrix
sampleSize <- 100
X <- matrix(0, nrow=sampleSize, ncol=25)
for(j in 1:25) {
  set.seed (j)
  X [, j] <- rnorm(sampleSize)
}

results <- randomFourierTrans (X=X, Dim=5, sigma=1, seedW=0)
stopifnot(all(dim(results$Z)==c(5*2, sampleSize)))

# Test fourierTransPredict
newX <- matrix(0, nrow=sampleSize, ncol=25)
for(j in 1:25) {
  set.seed (j+25)
  newX [, j] <- rnorm(sampleSize)
}
newResults <- fourierTransPredict (newx=newX, rW=results$rW)
stopifnot(all(dim(newResults)==c(5*2, sampleSize)))

# Check convergence of random Fourier transformation to radial basis function matrix
# 1. Higher dimension should yield better approximation
library(mvtnorm)
set.seed(10)
realDat <- rmvnorm(n=100, mean=rep(0, 10))

# Sigma=1
# Gaussian radial basis function matrix
sigmaSpec <- 1
distances <- dist(realDat)
gaussianRBFmatrix <- exp(-as.matrix(distances)^2 / (2*sigmaSpec))
# Calculate random Fourier transform
rftMat <- randomFourierTrans (X=realDat, Dim=1, sigma=sigmaSpec, seedW=0)
tZZ <- crossprodRcpp(rftMat$Z)[[1]]
# Compare with rbf matrix
approxSigma1Dim1 <- mean(abs(tZZ-gaussianRBFmatrix)) # 0.6277217
rftMat <- randomFourierTrans (X=realDat, Dim=10, sigma=sigmaSpec, seedW=0)
tZZ <- crossprodRcpp(rftMat$Z)[[1]]
approxSigma1Dim10 <- mean(abs(tZZ-gaussianRBFmatrix)) # 0.1776484
stopifnot(approxSigma1Dim1 > approxSigma1Dim10)
rftMat <- randomFourierTrans (X=realDat, Dim=100, sigma=sigmaSpec, seedW=0)
tZZ <- crossprodRcpp(rftMat$Z)[[1]]
approxSigma1Dim100 <- mean(abs(tZZ-gaussianRBFmatrix)) # 0.05640924
stopifnot(approxSigma1Dim10 > approxSigma1Dim100)

# Sigma=10
# Gaussian radial basis function matrix
sigmaSpec <- 10
distances <- dist(realDat)
gaussianRBFmatrix <- exp(-as.matrix(distances)^2 / (2*sigmaSpec))
# Calculate random Fourier transform
rftMat <- randomFourierTrans (X=realDat, Dim=1, sigma=sigmaSpec, seedW=0)
tZZ <- crossprodRcpp(rftMat$Z)[[1]]
# Compare with rbf matrix
approxSigma1Dim1 <- mean(abs(tZZ-gaussianRBFmatrix)) # 0.525662
rftMat <- randomFourierTrans (X=realDat, Dim=10, sigma=sigmaSpec, seedW=0)
tZZ <- crossprodRcpp(rftMat$Z)[[1]]
approxSigma1Dim10 <- mean(abs(tZZ-gaussianRBFmatrix)) # 0.1511192
stopifnot(approxSigma1Dim1 > approxSigma1Dim10)
rftMat <- randomFourierTrans (X=realDat, Dim=100, sigma=sigmaSpec, seedW=0)
tZZ <- crossprodRcpp(rftMat$Z)[[1]]
approxSigma1Dim100 <- mean(abs(tZZ-gaussianRBFmatrix)) # 0.04587967
stopifnot(approxSigma1Dim10 > approxSigma1Dim100)

# Sigma=100
# Gaussian radial basis function matrix
sigmaSpec <- 100
distances <- dist(realDat)
gaussianRBFmatrix <- exp(-as.matrix(distances)^2 / (2*sigmaSpec))
# Calculate random Fourier transform
rftMat <- randomFourierTrans (X=realDat, Dim=1, sigma=sigmaSpec, seedW=0)
tZZ <- crossprodRcpp(rftMat$Z)[[1]]
# Compare with rbf matrix
approxSigma1Dim1 <- mean(abs(tZZ-gaussianRBFmatrix)) # 0.09362182
rftMat <- randomFourierTrans (X=realDat, Dim=10, sigma=sigmaSpec, seedW=0)
tZZ <- crossprodRcpp(rftMat$Z)[[1]]
approxSigma1Dim10 <- mean(abs(tZZ-gaussianRBFmatrix)) # 0.03064888
stopifnot(approxSigma1Dim1 > approxSigma1Dim10)
rftMat <- randomFourierTrans (X=realDat, Dim=100, sigma=sigmaSpec, seedW=0)
tZZ <- crossprodRcpp(rftMat$Z)[[1]]
approxSigma1Dim100 <- mean(abs(tZZ-gaussianRBFmatrix)) # 0.01002129
stopifnot(approxSigma1Dim10 > approxSigma1Dim100)

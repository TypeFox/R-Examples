test.vgL3fitNM <- function () {
   for (i in 1:nrow(testParam)) {
    param <- testParam[i,]
    vgC <- param[1]
    sigma <- param[2]
    theta <- param[3]
    nu <- param[4]
    vgCHat <- vector(length=Nfit)
    sigmaHat <- vector(length=Nfit)
    thetaHat <- vector(length=Nfit)
    nuHat <- vector(length=Nfit)
    
    for (j in 1:Nfit) {
      vgFitReturn <- vgFit(rvg(n = n, param = param))$param
      vgCHat[j] <- vgFitReturn[1]
      sigmaHat[j] <- vgFitReturn[2]
      thetaHat[j] <- vgFitReturn[3]
      nuHat[j] <- vgFitReturn[4]
    }
    vgCHatAvg <- mean(vgCHat)
    sigmaHatAvg <- mean(sigmaHat)
    thetaHatAvg <- mean(thetaHat)
    nuHatAvg <- mean(nuHat)
    
    checkTrue(abs(vgCHatAvg - vgC) < errorThresholdFit,
      msg = paste(param[1], param[2], param[3], param[4], vgCHatAvg))
    checkTrue(abs(sigmaHatAvg - sigma) < errorThresholdFit,
      msg = paste(param[1], param[2], param[3], param[4], sigmaHatAvg))
    checkTrue(abs(thetaHatAvg - theta) < errorThresholdFit,
      msg = paste(param[1], param[2], param[3], param[4], thetaHatAvg))
    checkTrue(abs(nuHatAvg - nu) < errorThresholdFit,
      msg = paste(param[1], param[2], param[3], param[4], nuHatAvg))
    }
  }
  
  test.vgL3fitnlm <- function () {
   for (i in 1:nrow(testParam)) {
    param <- testParam[i,]
    vgC <- param[1]
    sigma <- param[2]
    theta <- param[3]
    nu <- param[4]
    vgCHat <- vector(length=Nfit)
    sigmaHat <- vector(length=Nfit)
    thetaHat <- vector(length=Nfit)
    nuHat <- vector(length=Nfit)

    for (j in 1:Nfit) {
      vgFitReturn <- vgFit(rvg(n = n, param = param), method = "nlm")$param
      vgCHat[j] <- vgFitReturn[1]
      sigmaHat[j] <- vgFitReturn[2]
      thetaHat[j] <- vgFitReturn[3]
      nuHat[j] <- vgFitReturn[4]
    }
    vgCHatAvg <- (cumsum(vgCHat)[Nfit])/Nfit
    sigmaHatAvg <- (cumsum(sigmaHat)[Nfit])/Nfit
    thetaHatAvg <- (cumsum(thetaHat)[Nfit])/Nfit
    nuHatAvg <- (cumsum(nuHat)[Nfit])/Nfit

    checkTrue(abs(vgCHatAvg - vgC) < errorThresholdFit,
      msg = paste(param[1], param[2], param[3], param[4], vgCHatAvg))
    checkTrue(abs(sigmaHatAvg - sigma) < errorThresholdFit,
      msg = paste(param[1], param[2], param[3], param[4], sigmaHatAvg))
    checkTrue(abs(thetaHatAvg - theta) < errorThresholdFit,
      msg = paste(param[1], param[2], param[3], param[4], thetaHatAvg))
    checkTrue(abs(nuHatAvg - nu) < errorThresholdFit,
      msg = paste(param[1], param[2], param[3], param[4], nuHatAvg))
    }
  }
  
  test.vgL3fitBFGS <- function () {
   for (i in 1:nrow(testParam)) {
    param <- testParam[i,]
    vgC <- param[1]
    sigma <- param[2]
    theta <- param[3]
    nu <- param[4]
    vgCHat <- vector(length=Nfit)
    sigmaHat <- vector(length=Nfit)
    thetaHat <- vector(length=Nfit)
    nuHat <- vector(length=Nfit)

    for (j in 1:Nfit) {
      vgFitReturn <- vgFit(rvg(n = n, param = param), method = "BFGS")$param
      vgCHat[j] <- vgFitReturn[1]
      sigmaHat[j] <- vgFitReturn[2]
      thetaHat[j] <- vgFitReturn[3]
      nuHat[j] <- vgFitReturn[4]
    }
    vgCHatAvg <- (cumsum(vgCHat)[Nfit])/Nfit
    sigmaHatAvg <- (cumsum(sigmaHat)[Nfit])/Nfit
    thetaHatAvg <- (cumsum(thetaHat)[Nfit])/Nfit
    nuHatAvg <- (cumsum(nuHat)[Nfit])/Nfit

    checkTrue(abs(vgCHatAvg - vgC) < errorThresholdFit,
      msg = paste(param[1], param[2], param[3], param[4], vgCHatAvg))
    checkTrue(abs(sigmaHatAvg - sigma) < errorThresholdFit,
      msg = paste(param[1], param[2], param[3], param[4], sigmaHatAvg))
    checkTrue(abs(thetaHatAvg - theta) < errorThresholdFit,
      msg = paste(param[1], param[2], param[3], param[4], thetaHatAvg))
    checkTrue(abs(nuHatAvg - nu) < errorThresholdFit,
      msg = paste(param[1], param[2], param[3], param[4], nuHatAvg))
    }
  }
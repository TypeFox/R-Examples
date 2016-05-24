test.vgL3momentsMean <- function () {
     for (i in 1:nrow(testParam)) {
      param <- testParam[i,]
      ### Calculate mean first
      mn <- vgMom(1, param = param, about = 0)

      # Calculate theoretical variance of sample means
      m1 <- vgMom(2, param = param, about = mn)
      theoVar <- m1/N
      # Calculate theoretical standard error of sample means
      theoStaError <- sqrt(theoVar)

      # Get N set of random numbers with n random numbers in
      # each set
      sampleMean <- vector(length = N)
      for (j in 1:N) {
         x <- rvg(n, param = param)
         # Compute mean of each sample data:
         sampleMean[j] <- mean(x)
      }

      # Get mean value from vgMean function:
      funMean <- vgMean(param = param)

      # compute sample error of sample means from the random samples
      sampStaError <- sqrt(var(sampleMean - funMean)/N)

      # Sample precision within the theoretical accuracy value?
      checkTrue(abs(sampStaError - theoStaError) < errorThresholdM, 
        msg = paste(param[1], param[2], param[3], param[4]))
  }
}


test.vgL3momentsVar <- function () {
   for (i in 1:nrow(testParam)) {
      param <- testParam[i,]
      ### Calculate mean first
      mn <- vgMom(1, param = param, about = 0)

      # Calculate theoretical variance of sample variances
      m2 <- (vgMom(4, param = param, about = mn) -
        (vgMom(2, param = param, about = mn))^2)
      theoVar <- m2/N
      # Calculate theoretical standard error of sample means
      theoStaError <- sqrt(theoVar)

      # Get N set of random numbers with n random numbers in
      # each set
      sampleVar <- vector(length = N)
      for (j in 1:N) {
         x <- rvg(n, param = param)
         # Compute variance of each sample data:
         sampleVar[j] <- var(x)
      }

      # Get variance value from vgVar function:
      funVar <- vgVar(param = param)

      # compute sample error of sample variances from the random samples
      sampStaError <- sqrt(var(sampleVar - funVar)/N)

      # Sample precision within the theoretical accuracy value?
      checkTrue(abs(sampStaError - theoStaError) < errorThresholdV, 
        msg = paste(param[1], param[2], param[3], param[4]))
  }
}

test.vgL3momentsSkew <- function () {
   for (i in 1:nrow(testParam)) {
      param <- testParam[i,]
      ### Calculate mean first
      mn <- vgMom(1, param = param, about = 0)

      # Calculate theoretical variance of sample skewness
      m3 <- vgMom(6, param = param, about = mn) -
        (vgMom(3, param = param, about = mn)^2) -
        6*vgMom(4, param = param, about = mn)*
        vgMom(2, param = param, about = mn) +
        9*(vgMom(2, param = param, about = mn)^3)
      m2 <- (vgMom(4, param = param, about = mn) -
        (vgMom(2, param = param, about = mn))^2)

      theoVar <- (m3/m2^(3/2))/N

      # Calculate theoretical standard error of sample skewness
      theoStaError <- sqrt(theoVar)

      # Get N set of random numbers with n random numbers in
      # each set
      sampleSkew <- vector(length = N)
      for (j in 1:N) {
         x <- rvg(n, param = param)
         # Compute variance of each sample data:
         sampleSkew[j] <- skewness(x)
      }

      # Get skewness value from vgSkew function:
      funSkew <- vgSkew(param = param)

      # compute sample error of sample skewnesses from the random samples
      sampStaError <- sqrt(var(sampleSkew - funSkew)/N)

      # Sample precision within the theoretical accuracy value?
      checkTrue(abs(sampStaError - theoStaError) < errorThresholdS, 
        msg = paste(param[1], param[2], param[3], param[4]))
  }
}

test.vgL3momentsKurt <- function () {
   for (i in 1:nrow(testParam)) {
      param <- testParam[i,]
      ### Calculate mean first
      mn <- vgMom(1, param = param, about = 0)

      # Calculate theoretical variance of sample kurtosis
      m4 <- vgMom(8, param = param, about = mn) -
        (vgMom(4, param = param, about = mn)^2) -
        8*vgMom(5, param = param, about = mn)*
        vgMom(3, param = param, about = mn) +
        16*vgMom(2, param = param, about = mn)*
        (vgMom(3, param = param, about = mn)^2)
      m2 <- (vgMom(4, param = param, about = mn) -
        (vgMom(2, param = param, about = mn))^2)

      theoVar <- (m4/m2^2)/N

      # Calculate theoretical standard error of sample kurtosis
      theoStaError <- sqrt(theoVar)

      # Get N set of random numbers with n random numbers in
      # each set
      sampleKurt <- vector(length = N)
      for (j in 1:N) {
         x <- rvg(n, param = param)
         # Compute variance of each sample data:
         sampleKurt[j] <- kurtosis(x)
      }

      # Get skewness value from vgSkew function:
      funKurt <- vgKurt(param = param)

      # compute sample error of sample skewnesses from the random samples
      sampStaError <- sqrt(var(sampleKurt - funKurt)/N)

      # Sample precision within the theoretical accuracy value?
      checkTrue(abs(sampStaError - theoStaError) < errorThresholdK, 
        msg = paste(param[1], param[2], param[3], param[4]))
  }
}

test.vgL3momentsMom <- function () {
  for (i in 1:nrow(testParam)) {
    param <- testParam[i,]
    orderSet <- c(1:8)
    ## raw moments
    for (j in 1:length(orderSet)) {
      momInte <- momIntegrated(densFn ="vg", order = j , param = param,
                               about = 0)
      momFun <- vgMom(order = j, param = param, momType = "raw")
      checkTrue(abs(momInte - momFun) < errorThresholdMom, 
        msg = paste(param[1], param[2], param[3], param[4], j))
    }
    
    ## mu moments
    for (j in 1:length(orderSet)) {
      momInte <- momIntegrated(densFn ="vg", order = j , param = param,
                               about = param[1])
      momFun <- vgMom(order = j, param = param, momType = "mu")
      checkTrue(abs(momInte - momFun) < errorThresholdMom, 
        msg = paste(param[1], param[2], param[3], param[4], j))
    }
    
    ## central moments
    for (j in 1:length(orderSet)) {
      ### Calculate mean first
      mn <- vgMom(1, param = param, about = 0)
      momInte <- momIntegrated(densFn ="vg", order = j , param = param,
                               about = mn)
      momFun <- vgMom(order = j, param = param, momType = "central")
      checkTrue(abs(momInte - momFun) < errorThresholdMom, 
        msg = paste(param[1], param[2], param[3], param[4], j))
    }
  }
}

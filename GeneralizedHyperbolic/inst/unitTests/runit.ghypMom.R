### Unit tests of function ghypMom

### Functions with name test.* are run by R CMD check or by make if
### LEVEL=1 in call to make
### Functions with name levelntest.* are run by make if
### LEVEL=n in call to make
### Functions with name graphicstest.* are run by make if
### LEVEL=graphics in call to make

test.ghypMom <- function()
{
  ## Purpose: Level 1 test of ghypMom
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: David Scott, Date: 11 Nov 2010, 13:02

  ## select a random parameter value for testing
  ## eliminate problem values first
  data(ghypParam)
  testParam <- ghypSmallParam
  paramSampleSize <- 1

  ## sample parameter values
  np <- NROW(testParam)
  paramNum <- sample(1:np, paramSampleSize, replace = FALSE)
  paramVals <- matrix(testParam[paramNum,], ncol = 5)
  nv <- NROW(paramVals)

  ## specify orders of moments
  orders <- c(1,2,3,4,5)

  ## mu moments first
  ## initialize result matrix
  results <- matrix(nrow = nv, ncol = 5 + length(orders)*2)
  differences <- matrix(nrow = nv, ncol = 5 + length(orders))

  ## loop over param values
  for (i in 1:nv){
    param <- as.numeric(paramVals[i,])
    mu <- param[1]
    results[i, 1:5] <- param
    differences[i, 1:5] <- param
    for (j in 1:length(orders)){
      results[i, 5 + 2*j - 1] <-
        ghypMom(orders[j], param = param, momType = "mu")
      results[i, 5 + 2*j] <-
        momIntegrated("ghyp", order = orders[j], param = param, about = mu)
      if (abs(results[i, 5 + 2*j]) < 10^(-10)){
        differences[i, 5 + j] <-
          results[i, 5 + 2*j - 1] - results[i, 5 + 2*j]
      } else {
        differences[i, 5 + j] <-
          (results[i, 5 + 2*j - 1] - results[i, 5 + 2*j])/
          results[i, 5 + 2*j]
      }
    }
    maxDiff <- max(abs(differences[i, 5 + 1:length(orders)]))
    jMax <- which.max(abs(differences[i, 5 + 1:length(orders)]))
    checkTrue(maxDiff < 10^(-10),
            msg = paste("Maximum difference", maxDiff, "for param",
                        param[1], param[2], param[3], param[4], param[5],
                        "order", orders[jMax]))
  }

  ## raw moments next
  ## initialize result matrix
  results <- matrix(nrow = nv, ncol = 5 + length(orders)*2)
  differences <- matrix(nrow = nv, ncol = 5 + length(orders))

  ## loop over param values
  for (i in 1:nv){
    param <- as.numeric(paramVals[i,])
    mu <- param[1]
    results[i, 1:5] <- param
    differences[i, 1:5] <- param
    for (j in 1:length(orders)){
      results[i, 5 + 2*j - 1] <-
        ghypMom(orders[j], param = param, momType = "raw")
      results[i, 5 + 2*j] <-
        momIntegrated("ghyp", order = orders[j], param = param, about = 0)
      if (abs(results[i, 5 + 2*j]) < 10^(-10)){
        differences[i, 5 + j] <-
          results[i, 5 + 2*j - 1] - results[i, 5 + 2*j]
      } else {
        differences[i, 5 + j] <-
          (results[i, 5 + 2*j - 1] - results[i, 5 + 2*j])/
          results[i, 5 + 2*j]
      }
    }
    maxDiff <- max(abs(differences[i, 5 + 1:length(orders)]))
    jMax <- which.max(abs(differences[i, 5 + 1:length(orders)]))
    checkTrue(maxDiff < 10^(-10),
            msg = paste("Maximum difference", maxDiff, "for param",
                        param[1], param[2], param[3], param[4], param[5],
                        "order", orders[jMax]))
  }

  ## now central moments
  ## initialize result matrix
  results <- matrix(nrow = nv, ncol = 5 + length(orders)*2)
  differences <- matrix(nrow = nv, ncol = 5 + length(orders))

  ## loop over param values
  for (i in 1:nv){
    param <- as.numeric(paramVals[i,])
    mu <- param[1]
    m1 <- ghypMom(1, param = param, momType = "raw")
    results[i, 1:5] <- param
    differences[i, 1:5] <- param
    for (j in 1:length(orders)){
      results[i, 5 + 2*j - 1] <-
        ghypMom(orders[j], param = param, momType = "central")
      results[i, 5 + 2*j] <-
        momIntegrated("ghyp", order = orders[j], param = param, about = m1)
      if (abs(results[i, 5 + 2*j]) < 10^(-10)){
        differences[i, 5 + j] <-
          results[i, 5 + 2*j - 1] - results[i, 5 + 2*j]
      } else {
        differences[i, 5 + j] <-
          (results[i, 5 + 2*j - 1] - results[i, 5 + 2*j])/
          results[i, 5 + 2*j]
      }
    }
    maxDiff <- max(abs(differences[i, 5 + 1:length(orders)]))
    jMax <- which.max(abs(differences[i, 5 + 1:length(orders)]))
    checkTrue(maxDiff < 10^(-10),
            msg = paste("Maximum difference", maxDiff, "for param",
                        param[1], param[2], param[3], param[4], param[5],
                        "order", orders[jMax]))
  }

  return()
}

level2test.ghypMom <- function()
{
  ## Purpose: Level 2 test of ghypMom
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: David Scott, Date: 11 Nov 2010, 13:02

  ## select a random parameter value for testing
  ## eliminate problem values first
  data(ghypParam)
  paramVals <- ghypLargeParam

  ## test all parameter values
  nv <- NROW(paramVals)

  ## specify orders of moments
  orders <- c(1,2,3,10,50,51)

  ## mu moments first
  ## initialize result matrix
  results <- matrix(nrow = nv, ncol = 5 + length(orders)*2)
  differences <- matrix(nrow = nv, ncol = 5 + length(orders))

  ## loop over param values
  for (i in 1:nv){
    param <- as.numeric(paramVals[i,])
    mu <- param[1]
    results[i, 1:5] <- param
    differences[i, 1:5] <- param
    for (j in 1:length(orders)){
      results[i, 5 + 2*j - 1] <-
        ghypMom(orders[j], param = param, momType = "mu")
      tryInt <- try(momIntegrated("ghyp", order = orders[j],
                                  param = param, about = mu),
                    silent = TRUE)
      if (class(tryInt) == "try-error"){
        results[i, 5 + 2*j] <- results[i, 5 + 2*j - 1]
      } else {
        results[i, 5 + 2*j] <- tryInt
      }
      if (abs(results[i, 5 + 2*j]) < 10^(-10)){
        differences[i, 5 + j] <-
          results[i, 5 + 2*j - 1] - results[i, 5 + 2*j]
      } else {
        differences[i, 5 + j] <-
          (results[i, 5 + 2*j - 1] - results[i, 5 + 2*j])/
          results[i, 5 + 2*j]
      }
    }
    maxDiff <- max(abs(differences[i, 5 + 1:length(orders)]))
    jMax <- which.max(abs(differences[i, 5 + 1:length(orders)]))
    checkTrue(maxDiff < 10^(-7),
            msg = paste("Maximum difference", maxDiff, "for param",
                        param[1], param[2], param[3], param[4], param[5],
                        "order", orders[jMax]))
  }

  ## raw moments next
  ## initialize result matrix
  results <- matrix(nrow = nv, ncol = 5 + length(orders)*2)
  differences <- matrix(nrow = nv, ncol = 5 + length(orders))

  ## loop over param values
  for (i in 1:nv){
    param <- as.numeric(paramVals[i,])
    mu <- param[1]
    results[i, 1:5] <- param
    differences[i, 1:5] <- param
    for (j in 1:length(orders)){
      results[i, 5 + 2*j - 1] <-
        ghypMom(orders[j], param = param, momType = "raw")
      tryInt <- try(momIntegrated("ghyp", order = orders[j],
                                  param = param, about = 0),
                    silent = TRUE)
      if (class(tryInt) == "try-error"){
        results[i, 5 + 2*j] <- results[i, 5 + 2*j - 1]
      } else {
        results[i, 5 + 2*j] <- tryInt
      }
      if (abs(results[i, 5 + 2*j]) < 10^(-10)){
        differences[i, 5 + j] <-
          results[i, 5 + 2*j - 1] - results[i, 5 + 2*j]
      } else {
        differences[i, 5 + j] <-
          (results[i, 5 + 2*j - 1] - results[i, 5 + 2*j])/
          results[i, 5 + 2*j]
      }
    }
    maxDiff <- max(abs(differences[i, 5 + 1:length(orders)]))
    jMax <- which.max(abs(differences[i, 5 + 1:length(orders)]))
    checkTrue(maxDiff < 10^(-7),
            msg = paste("Maximum difference", maxDiff, "for param",
                        param[1], param[2], param[3], param[4], param[5],
                        "order", orders[jMax]))
  }

  ## now central moments
  ## initialize result matrix
  results <- matrix(nrow = nv, ncol = 5 + length(orders)*2)
  differences <- matrix(nrow = nv, ncol = 5 + length(orders))

  ## loop over param values
  for (i in 1:nv){
    param <- as.numeric(paramVals[i,])
    mu <- param[1]
    m1 <- ghypMom(1, param = param, momType = "raw")
    results[i, 1:5] <- param
    differences[i, 1:5] <- param
    for (j in 1:length(orders)){
      results[i, 5 + 2*j - 1] <-
        ghypMom(orders[j], param = param, momType = "central")
      tryInt <- try(momIntegrated("ghyp", order = orders[j],
                                  param = param, about = m1),
                    silent = TRUE)
      if (class(tryInt) == "try-error"){
        results[i, 5 + 2*j] <- results[i, 5 + 2*j - 1]
      } else {
        results[i, 5 + 2*j] <- tryInt
      }
      if (abs(results[i, 5 + 2*j]) < 10^(-10)){
        differences[i, 5 + j] <-
          results[i, 5 + 2*j - 1] - results[i, 5 + 2*j]
      } else {
        differences[i, 5 + j] <-
          (results[i, 5 + 2*j - 1] - results[i, 5 + 2*j])/
          results[i, 5 + 2*j]
      }
    }
    maxDiff <- max(abs(differences[i, 5 + 1:length(orders)]))
    jMax <- which.max(abs(differences[i, 5 + 1:length(orders)]))
    checkTrue(maxDiff < 10^(-7),
            msg = paste("Maximum difference", maxDiff, "for param",
                        param[1], param[2], param[3], param[4], param[5],
                        "order", orders[jMax]))
  }

  return()
}

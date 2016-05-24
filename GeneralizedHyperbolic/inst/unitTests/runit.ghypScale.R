### Unit tests of function ghypScale

### Functions with name test.* are run by R CMD check or by make if
### LEVEL=1 in call to make
### Functions with name levelntest.* are run by make if
### LEVEL=n in call to make
### Functions with name graphicstest.* are run by make if
### LEVEL=graphics in call to make

test.ghypScale <- function()
{
  ## Purpose: Level 1 test of hypStandPars
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: David Scott, Date:  4 Jan 2011, 06:11

  ## select a random parameter value for testing
  data(ghypParam)
  testParam <- ghypLargeParam
  paramSampleSize <- 1

  ## sample parameter values
  np <- NROW(testParam)
  paramNum <- sample(1:np, paramSampleSize, replace = FALSE)
  param <- testParam[paramNum,]

  ## get a random mean and sd
  newMean <- rnorm(1, 0, 5)
  newSD <- rgamma(1, 1, scale = 5)

  ## rescale
  newParam <- ghypScale(newMean, newSD, param = param)

  ## obtain new mean and sd
  newMeanTest <- ghypMean(param = newParam)
  newSDTest <- sqrt(ghypVar(param = newParam))
  checkTrue(abs(newMean - newMeanTest) < 10^(-13),
            msg = paste("newMean = ", newMean, "newMeanTest = ", newMeanTest,
                        "diff = ", abs(newMean - newMeanTest), "for param",
                        param[1], param[2], param[3], param[4], param[5]))
  checkTrue(abs(newSD - newSDTest) < 10^(-13),
            msg = paste("newSD = ", newSD, "newSDTest = ", newSDTest,
                        "diff = ", abs(newSD - newSDTest), "for param",
                        param[1], param[2], param[3], param[4], param[5]))
}

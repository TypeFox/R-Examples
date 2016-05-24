### Unit tests of function ghypStandPars

### Functions with name test.* are run by R CMD check or by make if
### LEVEL=1 in call to make
### Functions with name levelntest.* are run by make if
### LEVEL=n in call to make
### Functions with name graphicstest.* are run by make if
### LEVEL=graphics in call to make

test.ghypStandPars <- function()
{
  ## Purpose: Level 1 test of hypStandPars
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: David Scott, Date:  4 Jan 2011, 06:11

  ## select a random parameter value for testing
  ## eliminate problem values first
  data(ghypParam)
  testParam <- ghypLargeParam
  paramSampleSize <- 1

  ## sample parameter values
  np <- NROW(testParam)
  paramNum <- sample(1:np, paramSampleSize, replace = FALSE)
  param <- testParam[paramNum,]
  param <- ghypChangePars(1, 2, testParam[paramNum,])

  ## convert to standardized parameters
  paramStar <- ghypStandPars(param[3], param[4], param[5])
  checkTrue(abs(ghypMean(param = paramStar) - 0) < 10^(-12),
            msg = paste("paramStar = ", paramStar,
                        "ghypMean(param = paramStar) = ",
                        ghypMean(param = paramStar)))
  checkTrue(abs(ghypVar(param = paramStar) - 1) < 10^(-12),
            msg = paste("paramStar = ", paramStar,
                        "ghypVar(param = paramStar) = ",
                        ghypVar(param = paramStar)))
}

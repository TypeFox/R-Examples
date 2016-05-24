### Unit tests of function ghypChangePars

### Functions with name test.* are run by R CMD check or by make if
### LEVEL=1 in call to make
### Functions with name levelntest.* are run by make if
### LEVEL=n in call to make
### Functions with name graphicstest.* are run by make if
### LEVEL=graphics in call to make

test.ghypChangePars <- function()
{
  ## Purpose: Level 1 test of hypChangePars
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
  param1 <- testParam[paramNum,]

  ## get other parameterizations
  param2 <- ghypChangePars(1, 2, param1)
  param3 <- ghypChangePars(1, 3, param1)
  param4 <- ghypChangePars(1, 4, param1)
  param5 <- ghypChangePars(1, 5, param1)

  ## check all pairs
  param11 <- ghypChangePars(1, 1, ghypChangePars(1, 1, param1))
  checkTrue(max(abs(param11 - param1)) < 10^(-14),
            msg = paste("param1 =", param1, "param11 =", param11))
  param12 <- ghypChangePars(2, 1, ghypChangePars(1, 2, param1))
  checkTrue(max(abs(param12 - param1)) < 10^(-14),
            msg = paste("param1 =", param1, "param12 =", param12))
  param13 <- ghypChangePars(3, 1, ghypChangePars(1, 3, param1))
  checkTrue(max(abs(param13 - param1)) < 10^(-14),
            msg = paste("param1 =", param1, "param13 =", param13))
  param14 <- ghypChangePars(4, 1, ghypChangePars(1, 4, param1))
  checkTrue(max(abs(param14 - param1)) < 10^(-14),
            msg = paste("param1 =", param1, "param14 =", param14))
  param15 <- ghypChangePars(5, 1, ghypChangePars(1, 5, param1))
  checkTrue(max(abs(param15 - param1)) < 10^(-14),
            msg = paste("param1 =", param1, "param15 =", param15))

  param21 <- ghypChangePars(1, 2, ghypChangePars(2, 1, param2))
  checkTrue(max(abs(param21 - param2)) < 10^(-14),
            msg = paste("param1 =", param1, "param21 =", param21))
  param22 <- ghypChangePars(2, 2, ghypChangePars(2, 2, param2))
  checkTrue(max(abs(param22 - param2)) < 10^(-14),
            msg = paste("param2 =", param2, "param22 =", param22))
  param23 <- ghypChangePars(3, 2, ghypChangePars(2, 3, param2))
  checkTrue(max(abs(param23 - param2)) < 10^(-14),
            msg = paste("param2 =", param2, "param23 =", param23))
  param24 <- ghypChangePars(4, 2, ghypChangePars(2, 4, param2))
  checkTrue(max(abs(param24 - param2)) < 10^(-14),
            msg = paste("param2 =", param2, "param24 =", param24))
  param25 <- ghypChangePars(5, 2, ghypChangePars(2, 5, param2))
  checkTrue(max(abs(param25 - param2)) < 10^(-14),
            msg = paste("param2 =", param2, "param25 =", param25))

  param31 <- ghypChangePars(1, 3, ghypChangePars(3, 1, param3))
  checkTrue(max(abs(param31 - param3)) < 10^(-14),
            msg = paste("param1 =", param1, "param11 =", param11))
  param32 <- ghypChangePars(2, 3, ghypChangePars(3, 2, param3))
  checkTrue(max(abs(param32 - param3)) < 10^(-14),
            msg = paste("param3 =", param3, "param32 =", param32))
  param33 <- ghypChangePars(3, 3, ghypChangePars(3, 3, param3))
  checkTrue(max(abs(param33 - param3)) < 10^(-14),
            msg = paste("param3 =", param3, "param33 =", param33))
  param34 <- ghypChangePars(4, 3, ghypChangePars(3, 4, param3))
  checkTrue(max(abs(param34 - param3)) < 10^(-14),
            msg = paste("param3 =", param3, "param34 =", param34))
  param35 <- ghypChangePars(5, 3, ghypChangePars(3, 5, param3))
  checkTrue(max(abs(param35 - param3)) < 10^(-14),
            msg = paste("param3 =", param3, "param35 =", param35))

  param41 <- ghypChangePars(1, 4, ghypChangePars(4, 1, param4))
  checkTrue(max(abs(param41 - param4)) < 10^(-14),
            msg = paste("param1 =", param1, "param11 =", param11))
  param42 <- ghypChangePars(2, 4, ghypChangePars(4, 2, param4))
  checkTrue(max(abs(param42 - param4)) < 10^(-14),
            msg = paste("param4 =", param4, "param42 =", param42))
  param43 <- ghypChangePars(3, 4, ghypChangePars(4, 3, param4))
  checkTrue(max(abs(param43 - param4)) < 10^(-14),
            msg = paste("param4 =", param4, "param43 =", param43))
  param44 <- ghypChangePars(4, 4, ghypChangePars(4, 4, param4))
  checkTrue(max(abs(param44 - param4)) < 10^(-14),
            msg = paste("param4 =", param4, "param44 =", param44))
  param45 <- ghypChangePars(5, 4, ghypChangePars(4, 5, param4))
  checkTrue(max(abs(param45 - param4)) < 10^(-14),
            msg = paste("param4 =", param4, "param45 =", param45))

  param51 <- ghypChangePars(1, 5, ghypChangePars(5, 1, param5))
  checkTrue(max(abs(param51 - param5)) < 10^(-14),
            msg = paste("param1 =", param1, "param11 =", param11))
  param52 <- ghypChangePars(2, 5, ghypChangePars(5, 2, param5))
  checkTrue(max(abs(param52 - param5)) < 10^(-14),
            msg = paste("param5 =", param5, "param52 =", param52))
  param53 <- ghypChangePars(3, 5, ghypChangePars(5, 3, param5))
  checkTrue(max(abs(param53 - param5)) < 10^(-14),
            msg = paste("param5 =", param5, "param53 =", param53))
  param54 <- ghypChangePars(4, 5, ghypChangePars(5, 4, param5))
  checkTrue(max(abs(param54 - param5)) < 10^(-14),
            msg = paste("param5 =", param5, "param54 =", param54))
  param55 <- ghypChangePars(5, 5, ghypChangePars(5, 5, param5))
  checkTrue(max(abs(param55 - param5)) < 10^(-14),
            msg = paste("param5 =", param5, "param55 =", param55))
}

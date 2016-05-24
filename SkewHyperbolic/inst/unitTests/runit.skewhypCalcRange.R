### Unit tests of function skewhypCalcRange

### Functions with name test.* are run by R CMD check or by make if
### LEVEL=1 in call to make
### Functions with name levelntest.* are run by make if
### LEVEL=n in call to make
### Functions with name graphicstest.* are run by make if
### LEVEL=graphics in call to make

test.skewhypCalcRange <- function()
{
  ## Purpose: Level 1 test of skewhypCalcRange
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: David Scott, Date: 18 Mar 2010, 11:35

  ## Distribution function
  ## symmetric case
  param <- c(0,1,0,10)
  range <- skewhypCalcRange(param = param, density = FALSE)
  int <- c(round(integrate(dskewhyp, -Inf, range[1],
                           param = param)$value,7),
           round(integrate(dskewhyp, range[2], Inf,
                           param = param)$value,7))
  checkEquals(rep(10^(-5), 2), int,
              msg = paste(param[1], param[2], param[3], param[4]))

  ## positive mu
  param <- c(1,1,0,10)
  range <- skewhypCalcRange(param = param, density = FALSE)
  int <- c(round(integrate(dskewhyp, -Inf, range[1],
                           param = param)$value,7),
           round(integrate(dskewhyp, range[2], Inf,
                           param = param)$value,7))
  checkEquals(rep(10^(-5), 2), int)

  ## negative mu
  param <- c(-2,1,0,10)
  range <- skewhypCalcRange(param = param, density = FALSE)
  int <- c(round(integrate(dskewhyp, -Inf, range[1],
                           param = param)$value,7),
           round(integrate(dskewhyp, range[2], Inf,
                           param = param)$value,7))
  checkEquals(rep(10^(-5), 2), int)

  ## large nu
  param <- c(0,1,10,20)
  range <- skewhypCalcRange(param = param, density = FALSE)
  int <- c(round(integrate(dskewhyp, -Inf, range[1],
                           param = param)$value,7),
           round(integrate(dskewhyp, range[2], Inf,
                           param = param)$value,7))
  checkEquals(rep(10^(-5), 2), int,
              msg = paste(param[1], param[2], param[3], param[4]))

  ## moderate nu
  param <- c(0,1,10,10)
  range <- skewhypCalcRange(param = param, density = FALSE)
  int <- c(round(integrate(dskewhyp, -Inf, range[1],
                           param = param)$value,7),
           round(integrate(dskewhyp, range[2], Inf,
                           param = param)$value,7))
  checkEquals(rep(10^(-5), 2), int,
              msg = paste(param[1], param[2], param[3], param[4]))

  ## large nu,  -ve skew
  param <- c(0,1,-10,20)
  range <- skewhypCalcRange(param = param, density = FALSE)
  int <- c(round(integrate(dskewhyp, -Inf, range[1],
                           param = param)$value,7),
           round(integrate(dskewhyp, range[2], Inf,
                           param = param)$value,7))
  checkEquals(rep(10^(-5), 2), int,
              msg = paste(param[1], param[2], param[3], param[4]))

  ## negative skew
  param <- c(0,1,-10,10)
  range <- skewhypCalcRange(param = param, density = FALSE)
  int <- c(round(integrate(dskewhyp, -Inf, range[1],
                           param = param)$value,7),
           round(integrate(dskewhyp, range[2], Inf,
                           param = param)$value,7))
  checkEquals(rep(10^(-5), 2), int,
              msg = paste(param[1], param[2], param[3], param[4]))

  ## non-standard limits
  param <- c(2,1,-10,10)
  range <- skewhypCalcRange(param = param, tol= 10^(-7), density = FALSE)
  int <- c(round(integrate(dskewhyp, -Inf, range[1],
                           param = param)$value,7),
           round(integrate(dskewhyp, range[2], Inf,
                           param = param)$value,7))
  checkEquals(rep(10^(-7), 2), int,
              msg = paste(param[1], param[2], param[3], param[4]))




  ## Density  function
  ## symmetric case
  param <- c(0,1,0,10)
  range <- skewhypCalcRange(param = param)
  dens <- c(round(dskewhyp(range[1], param = param),7),
           round(dskewhyp(range[2], param = param),7))
  checkEquals(rep(10^(-5), 2), dens,
              msg = paste(param[1], param[2], param[3], param[4]))

  ## positive mu
  param <- c(1,1,0,10)
  range <- skewhypCalcRange(param = param)
  dens <- c(round(dskewhyp(range[1], param = param),7),
           round(dskewhyp(range[2], param = param),7))
  checkEquals(rep(10^(-5), 2), dens)

  ## negative mu
  param <- c(-2,1,0,10)
  range <- skewhypCalcRange(param = param)
  dens <- c(round(dskewhyp(range[1], param = param),7),
           round(dskewhyp(range[2], param = param),7))
  checkEquals(rep(10^(-5), 2), dens)

  ## large nu
  param <- c(0,1,10,20)
  range <- skewhypCalcRange(param = param)
  dens <- c(round(dskewhyp(range[1], param = param),7),
           round(dskewhyp(range[2], param = param),7))
  checkEquals(rep(10^(-5), 2), dens,
              msg = paste(param[1], param[2], param[3], param[4]))

  ## moderate nu
  param <- c(0,1,10,10)
  range <- skewhypCalcRange(param = param)
  dens <- c(round(dskewhyp(range[1], param = param),7),
           round(dskewhyp(range[2], param = param),7))
  checkEquals(rep(10^(-5), 2), dens,
              msg = paste(param[1], param[2], param[3], param[4]))

  ## small nu, extreme +ve skew
  param <- c(0,1,20,1)
  range <- skewhypCalcRange(param = param)
  dens <- c(round(dskewhyp(range[1], param = param),7),
           round(dskewhyp(range[2], param = param),7))
  checkEquals(rep(10^(-5), 2), dens,
        msg = paste(param[1], param[2], param[3], param[4]))

  ## large nu,  -ve skew
  param <- c(0,1,-10,20)
  range <- skewhypCalcRange(param = param)
  dens <- c(round(dskewhyp(range[1], param = param),7),
           round(dskewhyp(range[2], param = param),7))
  checkEquals(rep(10^(-5), 2), dens,
              msg = paste(param[1], param[2], param[3], param[4]))

  ## negative skew
  param <- c(0,1,-10,10)
  range <- skewhypCalcRange(param = param)
  dens <- c(round(dskewhyp(range[1], param = param),7),
           round(dskewhyp(range[2], param = param),7))
  checkEquals(rep(10^(-5), 2), dens,
              msg = paste(param[1], param[2], param[3], param[4]))

  ## small nu, extreme -ve skew
  param <- c(0,1,-20,1)
  range <- skewhypCalcRange(param = param)
  dens <- c(round(dskewhyp(range[1], param = param),7),
           round(dskewhyp(range[2], param = param),7))
  checkEquals(rep(10^(-5), 2), dens,
              msg = paste(param[1], param[2], param[3], param[4]))

  ## non-standard limits
  param <- c(2,1,-20,1)
  range <- skewhypCalcRange(param = param, tol= 10^(-7))
  dens <- c(round(dskewhyp(range[1], param = param),7),
           round(dskewhyp(range[2], param = param),7))
  checkEquals(rep(10^(-7), 2), dens,
              msg = paste(param[1], param[2], param[3], param[4]))




  return()
}

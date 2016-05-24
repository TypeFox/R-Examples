### Unit tests of function tsHessian

### Functions with name test.* are run by R CMD check or by make if
### LEVEL=1 in call to make
### Functions with name levelntest.* are run by make if
### LEVEL=n in call to make

test.tsHessian <- function()
{
  ## Purpose: Level 1 test of tsHessian
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: David Scott, Date: 12 May 2010, 19:36

  ## Consider Hessian of log(1 + x + 2y)
  ## Example from Lang: A Second Course in Calculus, p.74
  fun <- function(param){
    x <- param[1]
    y <- param[2]
    return(log(1 + x + 2*y))
  }

  ## True value of Hessian
  trueHessian <- matrix( c(-1,-2,
                           -2,-4), byrow = 2, nrow = 2)
  

  ## Value from tsHessian
  approxHessian <- tsHessian(c(0,0), fun = fun)
  maxDiff <- max(abs(trueHessian - approxHessian))
  checkTrue(maxDiff < 0.1)

  return()
}

### Unit tests of function is.wholenumber

### Functions with name test.* are run by R CMD check or by make if
### LEVEL=1 in call to make
### Functions with name levelntest.* are run by make if
### LEVEL=n in call to make
### Functions with name graphicstest.* are run by make if
### LEVEL=graphics in call to make

test.is.wholenumber <- function()
{
  ## Purpose: Level 1 test of is.wholenumber
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: David Local, Date: 25 Jan 2010, 16:25

  checkTrue(is.wholenumber(-3:5))
  checkTrue(!is.wholenumber(c(0,0.1,1.3,5)))
  checkTrue(is.wholenumber(-3:5 + .Machine$double.eps))
  checkTrue(!is.wholenumber(-3:5 + .Machine$double.eps^0.5))
  checkTrue(is.wholenumber(c(2L,3L)))
  checkTrue(!is.wholenumber(c("2L","3L")))
  checkTrue(!is.wholenumber(0i ^ (-3:3)))
  checkTrue(is.wholenumber(matrix(1:6, nrow = 3)))
  checkTrue(!is.wholenumber(list(-1:3,2:6)))
  checkTrue(!is.numeric(list(-1:3,2:6)))
  checkTrue(is.wholenumber(unlist(list(-1:3,2:6))))

  return()
}


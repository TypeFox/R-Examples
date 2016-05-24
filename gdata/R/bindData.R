### bindData.R
###------------------------------------------------------------------------
### What: Bind two data frames - code
### $Id$
### Time-stamp: <2008-12-30 22:01:00 ggorjan>
###------------------------------------------------------------------------

bindData <- function(x, y, common)
{
  ## --- Setup ---
  if(!is.data.frame(x)) stop("'x' must be a data frame")
  if(!is.data.frame(y)) stop("'y' must be a data frame")

  ## --- New data frame ---

  ## First add common column and a dataset indicator column
  z <- rbind(x[common], y[common])

  ## Other columns
  ## - remove common columns in x and y
  namesz <- names(z)
  otherx <- names(x)
  otherx <- otherx[!(otherx %in% namesz)]
  othery <- names(y)
  othery <- othery[!(othery %in% namesz)]

  ## - add all other columns but as a set for each input data frame
  rx <- nrow(x); cx <- length(otherx)
  ry <- nrow(y); cy <- length(othery)
  
  z <- cbind(z, rbind(x[otherx], matrix(rep(NA, times=(ry * cx)), nrow=ry, ncol=cx, dimnames=list(NULL, otherx))))
  z <- cbind(z, rbind(matrix(rep(NA, times=(rx * cy)), nrow=rx, ncol=cy, dimnames=list(NULL, othery)), y[othery]))

  z
}

###------------------------------------------------------------------------
### bindData.R ends here

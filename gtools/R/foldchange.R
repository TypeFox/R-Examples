# $Id: foldchange.R 625 2005-06-09 14:20:30Z nj7w $

foldchange <- function(num,denom)
  {
    ifelse(num >= denom, num/denom, -denom/num)
  }


# Compute foldchange from log-ratio values
logratio2foldchange <- function(logratio, base=2)
  {
    retval <- base^(logratio)
    retval <- ifelse(retval < 1, -1/retval, retval)
    retval
  }

# vice versa
foldchange2logratio <- function(foldchange, base=2)
  {
    retval <- ifelse( foldchange<0, 1/-foldchange, foldchange)
    retval <- log(retval,base)
    retval
  }


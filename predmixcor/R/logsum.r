# FUNCTIONS FOR OPERATING ON NUMBERS REPRESENTED BY THEIR LOGS.
#
# Written by Radford M. Neal, July 2001.


# ADD NUMBERS REPRESENTED BY THEIR LOGS.  Computes the log of the sum of
# the exponentials of its arguments.  The arguments may be vectors or matrices,
# in which case this operation is carried out separately on each element.
# The arguments must all have the same dimensions, or be scalar.  The 
# computation is done in a manner that guarantees that the result is valid 
# even when directly computing the exponentials would result in overflow or 
# underflow.

LOG.add.exp <- function (...)
{
  m <- pmax(...)
  x <- 0

  for (a in list(...))
  { x <- x + exp(a-m)
  }

  log(x) + m
}


# SUM VECTORS OF NUMBERS REPRESENTED BY THEIR LOGS.  Computes the log of the 
# sum of the exponentials of all the elements in all its arguments.  The 
# computation is done in a manner that guarantees that the result is valid 
# even when directly computing the exponentials would result in overflow or 
# underflow.

LOG.sum.exp <- function (...)
{
  m <- max(...)

  log(sum(exp(c(...)-m))) + m
}


# FIND SQUARED DIFFERENCE OF NUMBERS REPRESENTED BY THEIR LOGS.  Computes
# the log of the squared difference of the exponentials of its two arguments.
# The arguments may be vectors or matrices, in which case this operation is
# done element-by-element.  (One of the arguments may be a scalar in this case.)
# The computation is done in a manner that guarantees that the result is valid 
# even when directly computing the exponentials would result in overflow or 
# underflow.

LOG.sqdiff.exp <- function (a, b)
{
  m <- pmax(a,b)
  
  log((exp(a-m)-exp(b-m))^2) + 2*m
}


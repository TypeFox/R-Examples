
#
# create and test objects of type 'rv'
#

rv <- function(length=0) {
  if (is.numeric(length)) {
    x <- as.list(rep.int(NA, length))
  } else {
    stop("length must be numeric")
  }
  class(x) <- "rv"
  return(x)
}

is.rv <- function(x)
{
  inherits(x, "rv")
}

# DEBUG:1 ?? The permutation of !is.null(dim(sims)) is not done yet
# DEBUG:2 (must take permutation of row indices and then permute.


as.double.rv <- function(x, ...)
{
  simapply(x, as.double, ...)
}

as.logical.rv <- function(x, ...)
{
  simapply(x, as.logical, ...)
}

as.integer.rv <- function (x, ...)
{
  simapply(x, as.integer, ...)
}


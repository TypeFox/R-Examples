# adapted from head() and tail()

some <- function(x, ...) UseMethod("some")

some.default <- function(x, n=10, ...){
    len <- length(x)
    ans <- x[sort(sample(len, min(n, len)))]
    if (length(dim(x)) == 1)
        array(ans, n, list(names(ans)))
    else ans
    }

some.matrix <- function(x, n=10, ...){
  nr <- nrow(x)
  x[sort(sample(nr, min(n, nr))), , drop = FALSE]
  }

some.data.frame <- function(x, n=10, ...){
    nr <- nrow(x)
    x[sort(sample(nr, min(n, nr))), , drop=FALSE]
    }

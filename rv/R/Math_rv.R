
# rv-Math.R - standard math functions for the rv class

Math.rv <- function(x, ...) {
  # Componentwise operation
  X <- x # Preserve class and other attributes
  for (i in seq_along(x)) {
    x <- X[[i]] 
    X[[i]] <- NextMethod()
  }
  return(X)
}

# cumsum, cumprod, cummax, cummin

cumsum.rv <- function (x)
{
  simapply(x, cumsum)
}

cumprod.rv <- function (x)
{
  simapply(x, cumprod)
}

cummin.rv <- function (x)
{
  simapply(x, cummin)
}

cummax.rv <- function (x)
{
  simapply(x, cummax)
}


# ----------------
# end of rv-Math.R
# ----------------

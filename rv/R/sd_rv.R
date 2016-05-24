## ========================================================================
## sd.rv
## ========================================================================
##

sd.rv <- function (x, na.rm = FALSE) {
  if (! is.rvobj(x)) {
    return(stats::sd(x, na.rm=na.rm))
  }
  if (is.vector(x)) {
    sqrt(var.rv(x, na.rm = na.rm))
  } else {
    sqrt(var.rv(as.vector(x), na.rm = na.rm))
  }
}


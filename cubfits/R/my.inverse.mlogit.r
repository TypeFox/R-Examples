my.inverse.mlogit <- function(e, log = FALSE){
  ### e is a matrix with dimension # of observations by # of synonymous codons.
  ret <- .Call("invmlogit", e, nrow(e), ncol(e), PACKAGE = "cubfits")
  if(log){
    ret <- log(ret)
  }
  ret
} # End of my.inverse.mlogit().


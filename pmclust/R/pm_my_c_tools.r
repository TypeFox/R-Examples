W.plus.y <- function(W, y, nrow, ncol){
  if(!is.double(W)){
    stop("W should be in double.")
  }
  if(!is.double(y)){
    stop("y should be in double.")
  }

  ret <- .Call("W_plus_y", W, y, as.integer(nrow), as.integer(ncol),
               PACKAGE = "pmclust")
  dim(ret) <- c(nrow, ncol)
  ret
} # End of W.plus.y().

W.plus.y.k <- function(W, y, nrow, ncol, i.k){
  if(!is.double(W)){
    stop("W should be in double.")
  }
  if(!is.double(y)){
    stop("y should be in double.")
  }

  ret <- .Call("W_plus_y_k", W, y, as.integer(nrow), as.integer(ncol),
               as.integer(i.k), PACKAGE = "pmclust")
#  dim(ret) <- c(nrow, ncol)
  dim(ret) <- c(nrow, 1)
  ret
} # End of W.plus.y.k().


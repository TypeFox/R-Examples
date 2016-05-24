### cbindX.R
###------------------------------------------------------------------------
### What: Column-bind objects with different number of rows - code
### $Id: cbindX.R 1300 2008-08-05 11:47:18Z ggorjan $
### Time-stamp: <2008-08-05 13:39:14 ggorjan>
###------------------------------------------------------------------------

cbindX <- function(...)
{
  ## --- Setup ---

  x <- list(...)

  ## Are all objects matrices or data.frames?
  test <- sapply(x, function(z) is.matrix(z) | is.data.frame(z))
  if(any(!test)) stop("only matrices and data.frames can be used")

  ## Get maximum number of rows
  tmp <- sapply(x, nrow)
  maxi <- which.max(tmp)
  test <- tmp < tmp[maxi]

  ## --- Core ---

  ## Adding additional "empty" rows so that all objects have the same number of rows
  for(i in 1:length(tmp)) {
    if(test[i]) {
      add <- matrix(nrow=tmp[maxi] - tmp[i], ncol=ncol(x[[i]]))
      if(is.data.frame(x[[i]])) {
        add <- as.data.frame(add)
      }
      colnames(add) <- colnames(x[[i]])
      x[[i]] <- rbind(x[[i]], add)
    }
  }

  ## Column-bind all objects
  ret <- x[[1]]
  for(i in 2:length(tmp)) {
    ret <- cbind(ret, x[[i]])
  }

  ## --- Return ---
  ret
}

###------------------------------------------------------------------------
### cbindX.R ends here

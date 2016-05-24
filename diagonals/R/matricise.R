#' Matricise
#'
#' @name matricise
#' @param x a higher-order array (length(dim(x)) >= 3)
#' @param row_dim the input array dimension which should be added to the row dimension of the output matrix, the value has to be 3 or 4.
#' @param col_dim the input array dimension which should be added to the column dimension of the output matrix, the value has to be 3 or 4.
#' @return a matrix (length(dim(x)) == 2 )

matricise <- function(x, row_dim = c(NULL,3,4), col_dim = c(NULL,3,4) ) {

  # save dimensions
  dx <- dim(x)

  # create row and col dims
  if (!is.null(row_dim)) {
    dr <- dx[row_dim]
  } else {
    dr <- 1L
  }
  if (!is.null(col_dim)) {
    dc <- dx[col_dim]
  } else {
    dc <- 1L
  }

  # create intermediate array
  a <- array( dim=c(nrow=(dx[1]*dr), ncol=dx[2], dc ) )

  # sort in lists
  s1 <- split_vector( 1:(dx[1]*dr), steps=dx[dr] )

  if (row_dim == 3) {
    # fill output matrix
    for (i in 1:length(s1) ) {
      a[s1[[i]], ,] <- x[,,i,]
    }
  } else if (row_dim == 4) {
    # fill output matrix
    for (i in 1:length(s1) ) {
      a[s1[[i]], ,] <- x[,,,i]
    }
  } else {
    stop("row_dim is not a valid array dimension")
  }

  # use new dimensions
  dx <- dim(a)

  # create output matrix
  m <- matrix(nrow=dx[1],ncol=(dx[2]*dc) )

  # sort in lists
  s2 <- split_vector( 1:(dx[2]*dc), steps=dc )

  # fill output matrix
  for (i in 1:length(s2) ) {
    m[, s2[[i]] ] <- a[,,i]
  }

  # return output matrix
  return(m)

}

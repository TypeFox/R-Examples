#' Get upperbounds of edits, given the boundaries of all variables
#' @param E \code{editmatrix}
#' @param xlim \code{matrix} with columns lower and upper, and rows are variables (in same order as E)
#' @return matrix with upperbounds per edit and a possible value
#' @keywords internal
getUpperBounds <- function(E, xlim){
  #print(xlim)
  A_p <- A_m <- A <- getA(E)
  A_p[A < 0] <- 0
  A_m[A > 0] <- 0
  ub <- A_m %*% xlim[,1] + A_p %*% xlim[,2]
  lb <- A_m %*% xlim[,2] + A_p %*% xlim[,1]
  
  b <- getb(E)
  b[is.na(b)] <- lb[is.na(b)]
  
  ub <- cbind(ub=ub, b=b, dummy=b-ub)
  #colnames(ub) <- c("lim", "dummy")
  rownames(ub) <- rownames(E)
  ub
}


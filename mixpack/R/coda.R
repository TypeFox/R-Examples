#' Orthonormal basis for the Simplex space
#' 
#' Basis from the simplex space with $D$ components
#'
#' @param D numbre of components
#' @export
ilr_basis <- function(D) {
  lapply(1:(D - 1), function(i) {
    I <- i + 1
    l <- exp(1/sqrt(i * I))
    r <- 1/exp(sqrt(i/I))
    s_i <- i * l + r + D - I
    c(rep(l, i), r, rep(1, D - I))/s_i
  })
}

#' Coordinates for an orthonormal basis
#' 
#' Coordinates respect basis \code{\link{ilr_basis}}
#'
#' @param X compositional sample
#' @export
ilr_coordinates <- function(X) {
  D <- NCOL(X)
  r <- data.frame(do.call("cbind", lapply(1:(D - 1), function(i) {
    if (i == 1) {
      XN <- X[1]
    } else {
      XN <- X[1:i]
    }
    as.numeric(sqrt(i/(i + 1)) * log(apply(XN, 1, prod)^(1/i)/X[[i + 1]]))
  })))
  names(r) <- paste("coord", 1:(D - 1), sep = ".")
  r
} 

#' clr coordinates of compositional data
#'
#' @param X compositional sample
#' @export
clr_coordinates <- function(X) {
  lX = log(X)
  lX - apply(lX, 1, mean)
} 
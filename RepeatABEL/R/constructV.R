#' @title Constructs the (co)variance matrix for y
#'
#' @description
#' Constructs the (co)variance matrix for y.
#'
#' @param Z The incidence matrix for the random effects column binded with the Cholesky of the GRM 
#' @param RandC The number of columns in the two matrices combined in Z. 
#' @param ratio The ratios between random effect variances and the residual variance.
#' 
#' @author Lars Ronnegard
#' 
constructV <-
function(Z, RandC, ratio) {
  V <- diag(nrow(Z))
  indx <- 1
  for (i in 1:length(ratio)) {
    Z.tmp <- Z[,indx:(indx + RandC[i] - 1)]
    indx = indx + RandC[i]
    V <- V + tcrossprod(Z.tmp) * ratio[i]
 }
 return(V)
}

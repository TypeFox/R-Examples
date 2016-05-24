#' Is a logit model constant for all observations?
#'
#' Determines whether the specified logit model is constant for all
#' observations. If it is constant then only one integral needs to be computed.
#'
#' @param xmat data
#' @param g0model logit model
#' @param width transect width
#'
#' @return logical value
#' @author Jeff Laake
is.logistic.constant <- function(xmat,g0model,width){
  xmat$distance <- rep(width, nrow(xmat))
  zlist <- setcov(xmat, g0model)
  beta <- rep(1,ncol(zlist))
  logit1 <- beta %*% t(zlist)

  return(all(logit1[1]==logit1))
}

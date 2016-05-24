#' Plots the joint diagonalization. I.e. if it was successful the matrices should all
#' be diagonal.
#' 
#' @param estConnectivity Estimate for connectivity matrix returned by \code{backShift}.
#' @param X Data matrix
#' @param env Indicator of the experiment or the intervention type an 
#' observation belongs to (a numeric vector of length n).
#' @param whichEnv Indicator for the environment to be plotted.
#' @param main Optional title for plot; defaults to paste("Env.", whichEnv)
#'
plotDiagonalization <- function(estConnectivity, X, env, whichEnv, main = NULL){
  deltas <- computeDelta(X, env)$Delta
  p <- ncol(estConnectivity)
  K <- length(deltas)
  settings <- sort(unique(env))
  deltaIntVar <- vector("list", K)
  
  tmp <- (diag(p) - t(estConnectivity))%*%deltas[[which(settings == whichEnv)]]%*%t(diag(p) - t(estConnectivity))
  deltaIntVar <- abs(tmp)/max(abs(tmp))
  
  if(is.null(main)) paste("Env.", whichEnv) else main
  
  fields::image.plot(z = deltaIntVar, col=fields::tim.colors(), legend.shrink = 0.8,
                     main = main)
}

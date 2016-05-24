#' @title Iterative predictor weighting PLS (IPW-PLS)
#'
#' @description An iterative procedure for variable elimination.
#'
#' @param y vector of response values (\code{numeric} or \code{factor}).
#' @param X numeric predictor \code{matrix}.
#' @param ncomp integer number of components (default = 10).
#' @param no.iter the number of iterations (default = 10).
#' @param IPW.threshold threshold for regression coefficients (default = 0.1).
#'
#' @details This is an iterative elimination procedure where a measure of predictor 
#' importance is computed after fitting a PLSR model (with complexity chosen based
#' on predictive performance). The importance measure is used both to re-scale the 
#' original X-variables and to eliminate the least important variables before
#' subsequent model re-fitting
#'  
#' @return Returns a vector of variable numbers corresponding to the model 
#' having lowest prediction error.
#'
#' @author Tahir Mehmood, Kristian Hovde Liland, Solve Sæbø.
#'
#' @references M. Forina, C. Casolino, C. Pizarro Millan, Iterative predictor weighting
#'  (IPW) PLS: a technique for the elimination of useless predictors in regression problems,
#'  Journal of Chemometrics 13 (1999) 165-184.
#'
#' @seealso \code{\link{VIP}} (SR/sMC/LW/RC), \code{\link{filterPLSR}}, \code{\link{spa_pls}}, 
#' \code{\link{stpls}}, \code{\link{truncation}}, \code{\link{bve_pls}}, \code{\link{mcuve_pls}},
#' \code{\link{ipw_pls}}, \code{\link{ga_pls}}, \code{\link{rep_pls}}.
#'
#' @examples
#' data(gasoline, package = "pls")
#' with( gasoline, ipw_pls(octane, NIR) )
#'
#' @export
ipw_pls <- function(y, X, ncomp=10, no.iter=10, IPW.threshold=0.1){
  
  if(is.factor(y)) {
    modeltype <- "classification"
    tb <- as.numeric(names(table(y)))
  } else {
    modeltype <- "prediction"
    y <- scale(y)
  }
  #X<- scale(X)
  for(i in 1:no.iter){
    pls.object <- plsr(y ~ X, ncomp=ncomp, validation = "CV")
    Press    <- pls.object$valid$PRESS[1,]
    opt.comp <- which.min(Press)
    pls.fit  <- plsr(y ~ X, ncomp=opt.comp)
    RC <- pls.fit$coef[,1,opt.comp]	
    SD <- apply(X, 2, sd)
    X  <- X*RC
  }
  ipw.selection <- which(abs(RC) >= IPW.threshold)
  if(length(ipw.selection)<= (ncomp +1)){
    ipw.selection <- sort(RC,decreasing = TRUE, index.return = T)$ix [1:ncomp]
  }
  return(list(ipw.selection=simplify(ipw.selection)))
}

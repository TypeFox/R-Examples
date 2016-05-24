#' @title predict method for l2boost models.
#'
#' @description \code{\link{predict}} is a generic function for predictions from the results 
#' of various model fitting functions. 
#' 
#'  @details \code{\link{predict.l2boost}} takes the optional \emph{xnew} (equivalent \emph{newdata}) \code{\link{data.frame}}
#' and returns the model estimates from an \code{\link{l2boost}} object. If neither \emph{xnew} or \emph{newdata} are
#' provided, \code{\link{predict}} returns estimates for the \code{\link{l2boost}} training data set.
#'
#' By default, \code{\link{predict.l2boost}} returns the function estimates, unless type="coef" then the 
#' set of regression coefficients (beta) are returned from the \code{\link{l2boost}} object.
#'
#' @param object an l2boost objectect
#' @param xnew a new design matrix to fit with the l2boost object
#' @param newdata a new design matrix to fit with the l2boost object
#' @param type "fit" or "coef" determins the values returned. "fit" returns model estimates, "coef" returns the 
#' model coefficients
#' @param ... other arguments (currently not used)
#'
#' @return function estimates for type=fit, coefficient estimates for type=coef
#' \itemize{
#' \item{yhat}{vector of n function estimates from the final step M}
#' \item{yhat.path}{list of M function estimates, one  at each step m}
#' }
#' or
#' \itemize{
#' \item{coef}{vector of p beta coefficient estimates from final step M}          
#' \item{coef.stand}{vector of p standardized beta coefficient estimates from final step M}     
#' \item{coef.path}{list of vectors of p beta coefficient estimates, one for each step m}  
#' \item{coef.stand.path}{list of vectors of p standardized beta coefficient estimates, one for each step m}  
#' }
#'  
#' @seealso \code{\link{predict}} and \code{\link{l2boost}}, \code{\link{coef.l2boost}},  
#' \code{\link{fitted.l2boost}}, \code{\link{residuals.l2boost}} and \code{\link{cv.l2boost}}
#'
#' @examples
#' #--------------------------------------------------------------------------
#' # Example 1: Diabetes 
#' #  
#' # See Efron B., Hastie T., Johnstone I., and Tibshirani R. 
#' # Least angle regression. Ann. Statist., 32:407-499, 2004.
#' data(diabetes)
#' 
#' object <- l2boost(diabetes$x,diabetes$y, M=1000, nu=.01)
#'
#' # With no arguments returns the estimates at the full M from the training data.
#' prd <- predict(object)
#' prd$yhat
#' 
#' # at step m=600
#' prd$yhat.path[[600]]
#' 
#' # Also can return coefficient estimates. This is equivalent to \code{\link{coef.l2boost}}
#' cf <- predict(object, type="coef")
#' cf$coef
#' 
#' # at step m=600
#' cf$coef.path[[600]]
#' 
#' # Or used to predict new data, in this case a subset of training data
#' cbind(diabetes$y[1:5], predict(object, xnew=diabetes$x[1:5,])$yhat)
#' 
#' @method predict l2boost
#' @S3method predict l2boost
predict.l2boost <- function(object, xnew = NULL, type = c("fit", "coef"), newdata=xnew, ...) {
  if(inherits(object,"cv")) object <- object$fit
  
  type <- match.arg(type)
  # extract necessary items from the stagewise objectect
  ybar <- object$ybar
  betam <- betam.stand <- object$betam
  betam.path <- betam.stand.path <- object$betam.path
  x.na <- object$x.na
  #cat(x.na, "\n")
  M <- length(betam.path)
  p <- length(betam)
  #cat(dim(xnew), "\n")
  
  if (any(object$x.na) & !is.null(xnew)) {
    xnew <- xnew[, !x.na, drop=FALSE]
  }
  if (is.null(xnew)) xnew <- object$x
  x.attr <- object$x.attr
  n <- x.attr$dim[1]
  n.new <- nrow(xnew)
  if (type == "fit") {
    # center and rescale xnew using original data
    # !!! do NOT remove NA columns !!!
    
    #cat(length(x.attr$"scaled:center"), "\t", dim(xnew), "\n")
    xnew <- scale(xnew, center = x.attr$"scaled:center", 
                  scale = x.attr$"scaled:scale")/sqrt(n - 1)
    
    #cat(xnew, "\n")
    yhat.path <- lapply(1:M, function(m) {
      be.m <- betam.stand.path[[m]]
      pt.m <- which(abs(be.m) > .Machine$double.eps)
      if (sum(pt.m) > 0) {
        yhat.m <- c(ybar + as.matrix(xnew[, pt.m]) %*% be.m[pt.m])
      }
      else {
        yhat.m <- rep(ybar, n.new)
      }
      yhat.m
    })
    yhat <- yhat.path[[M]]
    #return the predictor
    return(list(yhat = yhat, yhat.path = yhat.path))
  }
  else if (type == "coef") {
    # rescale coefficients back to original x-variable scale
    sf <- x.attr$"scaled:scale" *  sqrt(n - 1)
    betam <- betam.stand / sf
    for (k in 1:M) {
      betam.path[[k]] <- betam.stand.path[[k]] / sf
    }
    return(list(coef = betam,
                coef.stand = betam.stand,
                coef.path = betam.path,
                coef.stand.path = betam.stand.path))
  }
  else {
    stop("type must be set to 'fit' or 'coef'\n")
  }
}

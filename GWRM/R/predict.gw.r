#' GWRM Predictions
#'
#' Obtains predictions from a fitted GWRM object.
#'
#' @param object	a fitted object of class inheriting from \code{"gw"}.
#' @param newdata	optionally, a data frame in which to look for variables with which to predict. If omitted, the fitted linear predictors are used.
#' @param ... 	further arguments passed to or from other methods.
#'
#' @return A data frame with newdata and their fitted means.
#'
#' @importFrom stats delete.response .checkMFClasses
#'
#' @examples
#' data(goals)
#' fit <- gw(goals ~ position, data = goals)
#' predict(fit)
#'
#' @export

predict.gw <- function(object = NULL, newdata = NULL, ...){
  tt <- terms(object)
  if (!inherits(object, "gw"))
    warning("calling predict.gw(<fake-gw-object>) ...")

  if (missing(newdata) || is.null(newdata)) {
    mm <- X <- model.matrix(object)
    mmDone <- TRUE
    offset <- object$offset
  }
  else {
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, xlev = object$xlevels)
    if (!is.null(cl <- attr(Terms, "dataClasses")))
      .checkMFClasses(cl, m)
    nobs<-nrow(as.matrix(m))
    X <- model.matrix(Terms, m, offset <- rep(0, nobs))
    if (!is.null(off.num <- attr(tt, "offset")))
      for (i in off.num) offset <- offset + eval(attr(tt, "variables")[[i + 1]], newdata)
    if (!is.null(object$call$offset))
      offset <- offset + eval(object$call$offset, newdata)
    mmDone <- FALSE
  }

  ncovars <- ncol(X)
  beta <- object$coefficients[1:(ncovars)]
  if(!object$kBool){
    k <- object$betaIIpars[1]
    ro <- object$betaIIpars[2]
  }
  else{
    k <- object$k
    ro <- object$betaIIpars[2]
  }

  if (is.null(offset))
    fits <- exp(X %*% beta)
  else
    fits <- exp(offset + X %*% beta)
  predictor <- cbind(fits)
  colnames(predictor) <- c("fit")
  ans <- data.frame(predictor)
  return(ans)
}

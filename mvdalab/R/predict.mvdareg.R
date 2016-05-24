predict.mvdareg <- function(object, newdata, ncomp = object$ncomp, na.action = na.pass, ...) 
{
  #Adapted from 'pls' package
  if (missing(newdata) || is.null(newdata)) {
    new.X <- model.matrix(object)
  } else if (is.matrix(newdata)) {
    if (ncol(newdata) != length(object$Xmeans)) 
      stop("'newdata' does not have the correct number of columns")
    new.X <- newdata
  }  else {
    Terms <- delete.response(terms(object))
    options(contrasts = c(object$contrasts, "contr.poly"))
    m <- model.frame(Terms, newdata, na.action = na.action)
    if (!is.null(cl <- attr(Terms, "dataClasses"))) 
      .checkMFClasses(cl, m)
    new.X <- no.intercept(model.matrix(Terms, m))
  }
  nobs <- dim(new.X)[1]
  if (!is.null(object$scale)) 
    new.X <- new.X/rep(object$scale, each = nobs)
    if (missing(newdata)) {
      object$iPreds[, ncomp]
    } else {
      B <- object$coefficients[, ncomp]
      new.X <- new.X - rep(object$Xmeans, each = nobs)
      pred <- new.X %*% B + object$Ymean
      if (missing(newdata) && !is.null(object$na.action)) 
        pred <- napredict(object$na.action, pred)
      pred 
    }

}

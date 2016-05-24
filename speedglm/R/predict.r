predict.speedglm <- function (object, newdata, type = c("link", "response"),
                              na.action = na.pass, ...) 
{
  type <- match.arg(type)
  if (missing(newdata)&is.null(object$linear.predictors))        
    warning("fitted values were not returned from the speedglm object: 
            use the original data by setting argument 'newdata' or refit 
            the model by specifying fitted=TRUE.")
  na.act <- object$na.action
  object$na.action <- NULL
  if (missing(newdata)) {
    pred <- switch(type,
                   link = object$linear.predictors, 
                   response = fitted(object)
    )
    if (!is.null(na.act)) pred <- napredict(na.act, pred)
  } else {
    pred <- predict.speedlm(object, newdata, 
                            type = "response", 
                            na.action = na.action)
    switch(type, response = {
      pred <- family(object)$linkinv(pred)
    }, link = )
  }
  pred
}


predict.speedlm <- function (object, newdata, na.action = na.pass, ...) 
{
  tt <- terms(object)
  if (!inherits(object, c("speedlm","speedglm"))) 
    warning("calling predict.speedlm(<fake-speedlm/speedglm-object>) ...")
  if (missing(newdata) || is.null(newdata)) {
    if(is.null(object$fitted.values)) 
      warning("fitted values were not returned from the speedglm object: 
              use the original data by setting argument 'newdata' or refit 
              the model by specifying fitted=TRUE.")
    return(object$fitted.values)
  }
  else {
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, na.action = na.action)#, xlev = object$xlevels)
    if (!is.null(cl <- attr(Terms, "dataClasses"))) 
      .checkMFClasses(cl, m)
    X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
    offset <- rep(0, nrow(X))
    if (!is.null(off.num <- attr(tt, "offset"))) 
      for (i in off.num) offset <- offset + eval(attr(tt, 
                                                      "variables")[[i + 1]], newdata)
    if (!is.null(object$call$offset)) 
      offset <- offset + eval(object$call$offset, newdata)
  }
  p <- object$rank            
  ord <- colnames(X)
  if (p < ncol(X) && !(missing(newdata) || is.null(newdata))) 
    warning("prediction from a rank-deficient fit may be misleading")
  beta <- object$coefficients 
  beta[is.na(beta)] <- 0
  predictor <- drop(X[, ord, drop = FALSE] %*% beta[ord])              
  if (!is.null(offset)) 
    predictor <- predictor + offset
  if (missing(newdata) && !is.null(na.act <- object$na.action)) 
    predictor <- napredict(na.act, predictor)
  predictor
}

family.speedglm <- function(object,...) {
  object$family
}

fitted.speedglm <- function(object,...) {
  return(family(object)$linkinv(object$linear.predictors))
}

fitted.speedlm <- function(object,...) {
  object$fitted.values
}






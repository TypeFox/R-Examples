`estimable.mlm` <-
  function (obj, cm, beta0, conf.int=NULL,  show.beta0, ...)
{
  coef <- coef(obj)
  ny <- ncol(coef)
  effects <- obj$effects
  resid <- obj$residuals
  fitted <- obj$fitted
  ynames <- colnames(coef)
  if (is.null(ynames)) {
    lhs <- obj$terms[[2]]
    if (mode(lhs) == "call" && lhs[[1]] == "cbind") 
      ynames <- as.character(lhs)[-1]
    else ynames <- paste("Y", seq(ny), sep = "")
  }
  value <- vector("list", ny)
  names(value) <- paste("Response", ynames)
  cl <- oldClass(obj)
  class(obj) <- cl[match("mlm", cl):length(cl)][-1]
  for (i in seq(ny)) {
    obj$coefficients <- coef[, i]
    obj$residuals <- resid[, i]
    obj$fitted.values <- fitted[, i]
    obj$effects <- effects[, i]
    obj$call$formula[[2]] <- obj$terms[[2]] <- as.name(ynames[i])
    value[[i]] <- estimable(obj, cm, beta0, conf.int=NULL,  show.beta0, ...)
  }
  class(value) <- "listof"
  value
}

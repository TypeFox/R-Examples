# Last modified 25 Nov 2009 for point marking
# 18 January 2012 added robust estimation from Pendergast and Sheather

inverseResponsePlot <- function(model, lambda=c(-1, 0, 1), robust=FALSE,
   xlab=NULL, ...)
    UseMethod("inverseResponsePlot")

inverseResponsePlot.lm <- function(model, lambda=c(-1, 0, 1), xlab=NULL, 
       labels = names(residuals(model)), ...) {
  mf <- model$model
  if (is.null(mf)) mf <- update(model, model=TRUE, method="model.frame")
  xlab <- if(is.null(xlab)) names(mf)[1]
  y <- mf[, 1]
  yhat <- predict(model)
  invTranPlot(y, yhat, lambda=lambda, xlab=xlab, labels=labels, ...)
}

invResPlot <- function(model, ...) UseMethod("inverseResponsePlot")


##########NEW
inverseResponsePlot.lm <- function(model, lambda=c(-1, 0, 1), robust=FALSE, 
       xlab=NULL, labels = names(residuals(model)), ...) {
  if(robust == TRUE){
    m <- model$call
    m[[1L]] <- as.name("rlm")
    model <- eval(m, parent.frame())
  }
  mf <- model$model
  if (is.null(mf)) mf <- update(model, model=TRUE, method="model.frame")
  
  xlab <- if(is.null(xlab)) names(mf)[1]
  y <- mf[, 1]
  yhat <- predict(model)
  invTranPlot(y, yhat, lambda=lambda, xlab=xlab, labels=labels, 
       robust=robust, ...)
}
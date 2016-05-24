## plot x and y,
## with straight line fit and display of squared residuals
regr1.plot <- function(x, y,
                       model=lm(y~x),
                       coef.model,
                       alt.function,
                       main="put a useful title here",
                       xlab=deparse(substitute(x)),
                       ylab=deparse(substitute(y)),
                       jitter.x=FALSE,
                       resid.plot=FALSE,
                       points.yhat=TRUE,
                       pch=16,
                       ..., length.x.set=51,
                       x.name,
                       pch.yhat=16,
                       cex.yhat=par()$cex*.7,
                       err=-1) {
  if (sum(!missing(model), !missing(coef.model), !missing(alt.function)) > 1)
    stop("Use only one of the arguments: model, coef.model, alt.function.")
  old.err <- par(err=err)
  if (jitter.x) x <- jitter(x)
  plot(x, y, xlab=xlab, ylab=ylab, main=main, pch=pch, ...)

  x.set <- seq(par()$usr[1], par()$usr[2], length=length.x.set)
  newdata <- data.frame(x=x.set)
  if (!missing(model)) {
    if (missing(x.name)) {
      tmp <- substitute(x)
      x.name <- switch(length(tmp),
                       tmp[[1]],
                       stop("Please use the 'x.name' argument."),
                       if.R(s=if (is.character(tmp[[3]])) tmp[[3]],
                            r=if (is.name(tmp[[3]])) as.character(tmp[[3]])),
                       stop("Please use the 'x.name' argument."))
    }
    names(newdata) <- x.name
  }

  y.hat <- predict(model)
  y.hat.set <- predict(model, newdata)

  if (!missing(alt.function)) {
    y.hat <- alt.function(x)
    y.hat.set <- alt.function(newdata$x)
  }

  if (!missing(coef.model)) {
    y.hat <- coef.model[1] + coef.model[2] * x
    y.hat.set <- coef.model[1] + coef.model[2] * newdata$x
  }

  lines(x.set, y.hat.set)
  
  if (points.yhat) points(y=y.hat, x=x, pch=pch.yhat, cex=cex.yhat)
  if (resid.plot != FALSE) ## can't simplify.  The options are character strings
    resid.squares(x, y, y.hat, resid.plot)
  par(old.err)
  invisible(NULL)
}

calplot <- function(object, 
  xlim = c("auto", "auto"), ylim = c("auto", "auto"), 
  xlab = "Concentration", ylab = "Response", alpha = 0.05,
  varfunc = NULL)
{
  UseMethod("calplot")
}

calplot.default <- function(object, 
  xlim = c("auto","auto"), ylim = c("auto","auto"), 
  xlab = "Concentration", ylab = "Response",
  alpha=0.05, varfunc = NULL)
{
  stop("Calibration plots only implemented for univariate lm objects.")
}

calplot.lm <- function(object,
  xlim = c("auto","auto"), ylim = c("auto","auto"), 
  xlab = "Concentration", ylab = "Response", alpha=0.05,
  varfunc = NULL)
{
  if (length(object$coef) > 2)
    stop("More than one independent variable in your model - not implemented")

  if (alpha <= 0 | alpha >= 1)
    stop("Alpha should be between 0 and 1 (exclusive)")

  m <- object
  level <- 1 - alpha
  y <- m$model[[1]]
  x <- m$model[[2]]
  if (xlim[1] == "auto") xlim[1] <- 0
  if (xlim[2] == "auto") xlim[2] <- max(x)
  xlim <- as.numeric(xlim)
  newdata <- list(
    x = seq(from = xlim[[1]], to = xlim[[2]], length=250))
  names(newdata) <- names(m$model)[[2]]
  if (is.null(varfunc)) {
    varfunc <- if (length(m$weights)) {
        function(variable) mean(m$weights)
      } else function(variable) rep(1,250)
  }
  pred.lim <- predict(m, newdata, interval = "prediction",
    level=level, weights.newdata = varfunc(m))
  conf.lim <- predict(m, newdata, interval = "confidence",
    level=level)
  yrange.auto <- range(c(0,pred.lim))
  if (ylim[1] == "auto") ylim[1] <- yrange.auto[1]
  if (ylim[2] == "auto") ylim[2] <- yrange.auto[2]
  plot(1,
    type = "n", 
    xlab = xlab,
    ylab = ylab,
    xlim = as.numeric(xlim),
    ylim = as.numeric(ylim)
  )
  points(x,y, pch = 21, bg = "yellow")
  matlines(newdata[[1]], pred.lim, lty = c(1, 4, 4), 
    col = c("black", "red", "red"))
  if (length(object$weights) > 0) {
    legend(min(x), 
      max(pred.lim, na.rm = TRUE), 
      legend = c("Fitted Line", "Confidence Bands"),
      lty = c(1, 3), 
      lwd = 2, 
      col = c("black", "green4"), 
      horiz = FALSE, cex = 0.9, bg = "gray95")
  } else {
  matlines(newdata[[1]], conf.lim, lty = c(1, 3, 3), 
    col = c("black", "green4", "green4"))
  legend(min(x), 
    max(pred.lim, na.rm = TRUE), 
    legend = c("Fitted Line", "Confidence Bands", 
        "Prediction Bands"), 
    lty = c(1, 3, 4), 
    lwd = 2, 
    col = c("black", "green4", "red"), 
    horiz = FALSE, cex = 0.9, bg = "gray95")
  }
}

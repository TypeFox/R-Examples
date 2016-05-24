lod <- function(object, ..., alpha = 0.05, beta = 0.05, method = "default", tol = "default")
{
  UseMethod("lod")
}

lod.default <- function(object, ..., alpha = 0.05, beta = 0.05, method = "default", tol = "default")
{
  stop("lod is only implemented for univariate lm objects.")
}

lod.lm <- function(object, ..., alpha = 0.05, beta = 0.05, method = "default", tol = "default")
{
  if (length(object$weights) > 0) {
    stop(paste(
      "\nThe detemination of a lod from calibration models obtained by",
      "weighted linear regression requires confidence intervals for",
      "predicted y values taking into account weights for the x values",
      "from which the predictions are to be generated.",
      "This is not supported by the internally used predict.lm method.",
      sep = "\n"
    ))
  }
  xname <- names(object$model)[[2]]
  xvalues <- object$model[[2]]
  yname <- names(object$model)[[1]]
  newdata <- data.frame(0)
  names(newdata) <- xname
  y0 <- predict(object, newdata, interval = "prediction", 
    level = 1 - 2 * alpha)
  yc <- y0[[1,"upr"]]
  if (method == "din") {
    y0.d <- predict(object, newdata, interval = "prediction",
      level = 1 - 2 * beta)
    deltay <- y0.d[[1, "upr"]] - y0.d[[1, "fit"]]
    lod.y <- yc + deltay 
    lod.x <- inverse.predict(object, lod.y)$Prediction
  } else {
    f <- function(x) {
      newdata <- data.frame(x)
      names(newdata) <- xname
      pi.y <- predict(object, newdata, interval = "prediction",
        level = 1 - 2 * beta)
      yd <- pi.y[[1,"lwr"]]
      (yd - yc)^2
    }
    if (tol == "default") tol = min(xvalues[xvalues !=0]) / 1000
    lod.x <- optimize(f, interval = c(0, max(xvalues) * 10), tol = tol)$minimum
    newdata <- data.frame(x = lod.x)
    names(newdata) <- xname
    lod.y <-  predict(object, newdata)
  }
  lod <- list(lod.x, lod.y)
  names(lod) <- c(xname, yname)
  return(lod)
}

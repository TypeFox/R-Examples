# This is an implementation of Equation (8.28) in the Handbook of Chemometrics
# and Qualimetrics, Part A, Massart et al (1997), page 200, validated with
# Example 8 on the same page, extended as specified in the package vignette

inverse.predict <- function(object, newdata, ...,
  ws = "auto", alpha = 0.05, var.s = "auto")
{
  UseMethod("inverse.predict")
}

inverse.predict.default <- function(object, newdata, ...,
  ws = "auto", alpha = 0.05, var.s = "auto")
{
  stop("Inverse prediction only implemented for univariate lm objects.")
}

inverse.predict.lm <- function(object, newdata, ...,
  ws = "auto", alpha = 0.05, var.s = "auto")
{
  yname = names(object$model)[[1]]
  xname = names(object$model)[[2]]
  if (ws == "auto") {
    ws <- ifelse(length(object$weights) > 0, mean(object$weights), 1)
  }
  if (length(object$weights) > 0) {
    wx <- split(object$weights,object$model[[xname]])
    w <- sapply(wx,mean)
  } else {
    w <- rep(1,length(split(object$model[[yname]],object$model[[xname]])))
  }
  .inverse.predict(object = object, newdata = newdata, 
    ws = ws, alpha = alpha, var.s = var.s, w = w, xname = xname, yname = yname)
}

inverse.predict.rlm <- function(object, newdata, ...,
  ws = "auto", alpha = 0.05, var.s = "auto")
{
  yname = names(object$model)[[1]]
  xname = names(object$model)[[2]]
  if (ws == "auto") {
    ws <- mean(object$w)
  }
  wx <- split(object$weights,object$model[[xname]])
  w <- sapply(wx,mean)
  .inverse.predict(object = object, newdata = newdata, 
    ws = ws, alpha = alpha, var.s = var.s, w = w, xname = xname, yname = yname)
}

.inverse.predict <- function(object, newdata, ws, alpha, var.s, w, xname, yname)
{
  if (length(object$coef) > 2)
    stop("More than one independent variable in your model - not implemented")

  if (alpha <= 0 | alpha >= 1)
    stop("Alpha should be between 0 and 1 (exclusive)")

  ybars <- mean(newdata)
  m <- length(newdata)

  yx <- split(object$model[[yname]], object$model[[xname]])
  n <- length(yx)
  df <- n - length(object$coef)
  x <- as.numeric(names(yx))
  ybar <- sapply(yx,mean)
  yhatx <- split(object$fitted.values, object$model[[xname]])
  yhat <- sapply(yhatx, mean)
  se <- sqrt(sum(w * (ybar - yhat)^2)/df)

  if (var.s == "auto") {
    var.s <- se^2/ws
  }

  b1 <- object$coef[[xname]]

  ybarw <- sum(w * ybar)/sum(w)

# This is the adapted form of equation 8.28 (see package vignette)
  sxhats <- 1/b1 * sqrt(
    (var.s / m) + 
    se^2 * (1/sum(w) + 
      (ybars - ybarw)^2 * sum(w) /
        (b1^2 * (sum(w) * sum(w * x^2) - sum(w * x)^2)))
  )

  if (names(object$coef)[1] == "(Intercept)") {
    b0 <- object$coef[["(Intercept)"]]
  } else {
    b0 <- 0
  }

  xs <- (ybars - b0) / b1
  t <- qt(1-0.5*alpha, n - 2)
  conf <- t * sxhats
  result <- list("Prediction"=xs,"Standard Error"=sxhats,
    "Confidence"=conf, "Confidence Limits"=c(xs - conf, xs + conf))
  return(result)
}

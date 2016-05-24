# Prediction data frame
# Get predictions with standard errors into data frame
#
# @keyword internal
# @alias predictdf2d.default
# @alias predictdf2d.glm
# @alias predictdf2d.loess
# @alias predictdf2d.locfit
predictdf2d <- function(model, xseq, yseq) UseMethod("predictdf2d")

#' @export
predictdf2d.default <- function(model, xseq, yseq ) {
  newdata = expand.grid(x=xseq,y=yseq)
  pred    = stats::predict(model, newdata = newdata, se.fit = FALSE, interval = "none")
  data.frame(newdata, z = as.vector(pred))
}

#' @export
predictdf2d.glm <- function(model, xseq, yseq) {
  newdata = expand.grid(x=xseq,y=yseq)
  pred    = stats::predict(model, newdata = newdata, se.fit = FALSE, type = "link")
  data.frame(newdata, z = model$family$linkinv(as.vector(pred)))
}

#' @export
predictdf2d.loess <- function(model, xseq, yseq ) {
  newdata = expand.grid(x=xseq,y=yseq)
  pred    = stats::predict(model,newdata, se = FALSE)
  data.frame(newdata, z = as.vector(pred))
}

#' @export
predictdf2d.locfit <- function(model, xseq, yseq ) {
  newdata = expand.grid(x=xseq,y=yseq)
  pred    = stats::predict(model, newdata = newdata, se.fit = FALSE)
  data.frame(newdata, z = as.vector(pred))
}



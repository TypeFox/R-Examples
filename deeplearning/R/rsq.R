#' Calculate the RSQ of a regression model
#' Utilitiy function that calcualtes RSQ of a model. It measures the goodness-of-
#' fit of a regression model.
#'
#' @param  x Regression Model
#' @param  ... Additional Input
#'
#' @import  futile.logger
#' @export

rsq <- function(x, ...) {
  UseMethod("rsq", x)
}

#' Utilitiy function that calcualtes RSQ of a DArch instance
#'
#' Calcualte a regression model's RSQ of a deep neural network
#'
#' @param  x DArch Model
#' @param  input Input data
#' @param  target Target data
#' @param ... addtional inputs
#' @import futile.logger
#' @importFrom stats predict
#' @importFrom graphics plot
#' @export

rsq.DArch <- function(x,
                      input = x@dataSet@data,
                      target = x@dataSet@targets, ...) {
  y <- target
  pred <- predict(x, newdata = input)
  plot(y, pred, xlab = "target", ylab = "prediction")
  RSQ <- 1 - sum((pred-y)^2)/sum((y-mean(y))^2)
  flog.info(paste0("RSQ = ", RSQ))
}

#' Utilitiy function that calcualtes RSQ of a linear model
#'
#' Calcualte a regression model's RSQ
#'
#' @param x linear Model
#' @param input Input data
#' @param target Target data
#' @param ... additional inputs
#' @importFrom  stats predict
#' @importFrom graphics plot
#' @import  futile.logger
#' @export

rsq.lm <- function(x, input, target, ...) {
  y <- target
  pred <- predict(x, newdata = data.frame(input))
  plot(y, pred)
  plot(y, pred, xlab = "target", ylab = "prediction")
  RSQ <- 1 - sum((pred-y)^2)/sum((y-mean(y))^2)
  flog.info(paste0("RSQ = ", RSQ))
}

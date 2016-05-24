#' Calculates the Accuracy Ratio of a classifier
#'
#' This function calculates the Accuracy Ratio of a binary classification
#'  model
#'
#'
#' @param x model
#' @param ... additional inputs
#'
#' @export

AR <- function(x, ...) {
  UseMethod("AR")
}


#' Calculates the Accruacy Ratio of a given set of probability
#'
#' This function calculates the Accuracy Ratio of a binary classification model
#'  output against its targets
#'
#' @param x a list of model output in the form of probabilities
#' @param target binary response
#' @param ... additional inputs
#' @export

AR.numeric <- function(x, target, ...) {
  AR.default(x, target)
}


#' Calculates the Accruacy Ratio of a given set of probability
#'
#' This function calculates the Accuracy Ratio of a binary classification model
#'  output against its targets
#'
#' @param x a list of model output in the form of probabilities
#' @param target binary response
#' @param ... additional inputs
#' @importFrom graphics plot
#'
#' @export

AR.default <- function(x, target, ...) {
  N <- length(x)
  seq = order(x, decreasing = T)
  target <- target[seq]
  auc <- 0
  totTarget <- sum(target)
  y <- c()
  for (i in 1:N) {
    lorenzeCurve <- sum(target[1:i]) / totTarget
    auc <- auc + lorenzeCurve * 1 / N
    y <- cbind(y, lorenzeCurve)
  }
  auc <- auc
  pd <- sum(target) / N
  ar <- (2 * auc - 1) / (1 - pd)
  plot(as.vector(y), xlab = "Population", ylab = "Fraction of Positive")
  if(ar > 1) ar <- 1
  return (ar)
}

#' Calculates the Accruacy Ratio of a given set of probability
#'
#' This function calculates the Accuracy Ratio of a trained darch instance
#'
#' @param x a DArch instance
#' @param input the input matrix
#' @param target binary response
#' @param ... additional inputs
#'
#' @importFrom stats predict
#'
#' @export



AR.DArch <- function(x, input = x@dataSet@data,
                     target = x@dataSet@targets, ...) {
  pred <- predict(x, newdata = input)
  AR.default(pred, target)
}



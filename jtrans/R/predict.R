#' Predict functions for Johnson Transformation
#' 
#' Generic functions to apply the fitted Johnson transformation to transform new
#' data.
#' 
#' After the johnson transformation is used, if you want to use the same
#' transformation on different data, you can use these functions. This is
#' designed to be the same functionality as the \code{\link{predict}} functions.
#' 
#' @param object a \code{jtrans} object with a specific type
#' @param newdata new data to be fitted
#' @param ... further arguments to match \code{\link{predict}}
#' @return Numeric vector of the transformed values 
#' 
#' #' @examples 
#' #' # if you want to predict based on a fitted distribution, you must set the
#' # parameters in the qtls() function using the fitted model object jt.
#' 
#' jt <- jtrans(rexp(300, .4))
#' 
#' # good prediction
#' predict(jt, rexp(10, .4))
#' 
#' # will generate NaN because newx is from different distribution 
#' predict(jt, rexp(10, .1))
#' 
#' 
#' 
#' @rdname predict_functions
#' @export
predict.sb <- function(object, newdata, ...) {
  gamma <- object$gamma
  eta <- object$eta
  epsilon <- object$epsilon
  lambda <- object$lambda
  return(gamma + eta * log((newdata - epsilon) / (lambda + epsilon - newdata)))
}

#' @rdname predict_functions
#' @export
predict.su <- function(object, newdata, ...) {
  gamma <- object$gamma
  eta <- object$eta
  epsilon <- object$epsilon
  lambda <- object$lambda
  return(gamma + eta * asinh((newdata - epsilon) / lambda))
}

#' @rdname predict_functions
#' @export
predict.sl <- function(object, newdata, ...) {
  gamma <- object$gamma
  eta <- object$eta
  epsilon <- object$epsilon
  return(gamma + eta * log(newdata - epsilon))
}


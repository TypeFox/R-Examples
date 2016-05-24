#' @title Parameter estimates in a PPSfit Object
#' 
#' @description \code{coef} returns the parameter estimates in a \code{PPSfit} Object
#' 
#' @param object A \code{PPSfit} object, typically from \code{PPS.fit()}.
#' @param \dots Other arguments.
#' 
#' @return A list with the parameter estimates.
#' 
#' @seealso \code{\link{PPS.fit}}.
#' 
#' @examples
#' x <- rPPS(50, 1.2, 100, 2.3)
#' fit <- PPS.fit(x)
#' coef(fit)

#' @export
coef.PPSfit <-
  function (object, ...) 
  {
    if (!class(object) == "PPSfit") stop("Object must belong to class PPS")            
    print.default(format(object$estimate), print.gap = 2, quote = FALSE)
    invisible(object)
  }
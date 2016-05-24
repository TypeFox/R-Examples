#' @title Log-likelihood value of a PPSfit Object
#' @description It returns the log-likelihood value of a \code{PPSfit} Object
#' 
#' @param object A \code{PPSfit} Object.
#' @param \dots Other arguments.
#' 
#' @return The log-likelihood.
#' 
#' @references Sarabia, J.M and Prieto, F. (2009). The Pareto-positive stable distribution: A new descriptive model for city size data, \emph{Physica A: Statistical Mechanics and its Applications}, \bold{388}(19), 4179-4191.
#'
#' @seealso \code{\link{PPS.fit}}
#' @examples
#' x <- rPPS(50, 1.2, 100, 2.3)
#' fit <- PPS.fit(x)
#' logLik(fit)

#' @export
logLik.PPSfit <-
  function (object, ...) 
  {
    val <- object$loglik
    attr(val, "nobs") <- object$n
    attr(val, "df") <- length(object$estimate)
    class(val) <- "logLik"
    val
  }
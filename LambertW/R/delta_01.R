#' @title Input parameters to get zero mean, unit variance output given delta
#' 
#' @description
#' 
#' Computes the input mean \eqn{\mu_x(\delta)} and standard deviation
#'     \eqn{\sigma_x(\delta)} for input \eqn{X \sim F(x \mid \boldsymbol \beta)}
#'     such that the resulting heavy-tail Lambert W x F RV \eqn{Y} with
#'     \eqn{\delta} has zero-mean and unit-variance.  So far works only for
#'     Gaussian input and scalar \eqn{\delta}.
#' 
#' The function works for any output mean and standard deviation, but default
#'     values are \eqn{\mu_y = 0} and \eqn{\sigma_y = 1} since they are the most
#'     useful, e.g., to generate a standardized Lambert W white noise sequence.
#'  
#' @param delta scalar; heavy-tail parameter.
#' @param mu.y output mean; default: \code{0}.
#' @param sigma.y output standard deviation; default: \code{1}.
#' @param distname string; distribution name.  Currently this function only supports
#' \code{"normal"}.
#' @return 
#' 5-dimensional vector (\eqn{\mu_x(\delta)}, \eqn{\sigma_x(\delta)}, 0, \eqn{\delta}, 1), 
#' where \eqn{\gamma = 0}  and \eqn{\alpha = 1} are set for the sake of compatiblity with other functions.
#' @keywords math univar
#' @export
#' @examples
#' 
#' delta_01(0) # for delta = 0, input == output, therefore (0,1,0,0,1)
#' # delta > 0 (heavy-tails): 
#' #   Y is symmetric for all delta: 
#' #   mean = 0; however, sd must be smaller 
#' delta_01(0.1) 
#' delta_01(1/3)  # only moments up to order 2 exist
#' delta_01(1)  # neither mean nor variance exist, thus NA

delta_01 <- function(delta, mu.y = 0, sigma.y = 1, distname = "normal") {
  stopifnot(is.numeric(delta), 
            length(delta) == 1,
            delta >= 0,  # for moments we only consider positive moments
            is.numeric(mu.y),
            length(mu.y) == 1,
            length(sigma.y) == 1,
            sigma.y >= 0)
  
  distname <- match.arg(distname)
  if (distname == "normal") {
    .moments_N01_input <- function(n = 1, delta = 0) {
      # n: moments
      if (n >= 1 / delta) {
        # n-th moment exists only for n < 1 / delta
        return(NA)
      }
      if (n %% 2) {
        return(0)
      } else {
        return(factorial(n) * (1 - n * delta)^(-(n + 1)/2)/(2^(n/2) * factorial(n/2)))
      }
    }
  } else {
    stop("Distribution", distname, "is not implemented yet.")
  }
  # moments for a standard Normal(0,1) input and 'delta' heavy tail parameter
  out <- c(mu_x = ifelse(delta < 1, 0, NA) + mu.y,
           sigma_x = sqrt(1/(.moments_N01_input(n = 2, delta))) * sigma.y, 
           gamma = 0, 
           delta = delta,
           alpha = 1)
  return(out)  
} 


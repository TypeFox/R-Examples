#' @title Required Sample Size
#' 
#' @description Function to determine required sample size to be within a given margin of error
#' 
#' @details Answer is based on a normal approximation when using type \code{"pi"}.
#' 
#' @param b the desired bound
#' @param sigma population standard deviation; not required if using type \code{"pi"}
#' @param p estimate for the population proportion of successes; not required if using type \code{"mu"}
#' @param conf.level confidence level for the problem, restricted to lie between zero and one
#' @param type character string, one of \code{"mu"} or \code{"pi"}, or just the initial letter of each, indicating the appropriate parameter; default value is \code{"mu"}
#' 
#' @author Alan T. Arnholt <arnholtat@@appstate.edu> 
#' 
#' @export
#' 
#' @examples
#' nsize(b = 0.015, p = 0.5, conf.level = 0.95, type = "pi")
#' # Returns the required sample size (n) to estimate the population 
#' # proportion of successes with a 0.95 confidence interval
#' # so that the margin of error is no more than 0.015 when the
#' # estimate of the population propotion of successes is 0.5.
#' nsize(b = 0.02, sigma = 0.1, conf.level = 0.95, type = "mu")
#' # Returns the required sample size (n) to estimate the population 
#' # mean with a 0.95 confidence interval so that the margin
#' # of error is no more than 0.02.
#' 
#' @keywords programming
#####################################################################################
nsize <- function(b, sigma = NULL, p = 0.5, conf.level = 0.95, type = c("mu", "pi"))
{
  type <- match.arg(type)
  if(length(type) > 1 || is.na(type))
    stop("type must be one \"mu\", \"pi\"")
  if(type == "pi" && b > 1)
    stop("b must be less than 1")
  if(!missing(b))
    if(length(b) != 1 || is.na(b))
      stop("b must be a single number")
  if(type == "mu") {
    z <- qnorm(1 - (1 - conf.level)/2)
    n <- ((z * sigma)/b)^2
    n <- ceiling(n)
    cat("\n")
    cat("The required sample size (n) to estimate the population",
        "\n")
    cat("mean with a", conf.level,
        "confidence interval so that the margin", "\n")
    cat("of error is no more than", b, "is", n,"\b.", "\n")
    cat("\n\n")
  }
  else if(type == "pi") {
    z <- qnorm(1 - (1 - conf.level)/2)
    n <- p * (1 - p) * (z/b)^2
    n <- ceiling(n)
    cat("\n")
    cat("The required sample size (n) to estimate the population",
        "\n")
    cat("proportion of successes with a", conf.level,
        "confidence interval", "\n")
    cat("so that the margin of error is no more than", b, "is",
        n,"\b.", "\n")
    cat("\n\n")
  }
}
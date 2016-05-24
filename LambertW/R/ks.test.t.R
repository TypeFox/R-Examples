#' @title One-sample Kolmogorov-Smirnov test for student-t distribution
#' 
#' @description
#' Performs a two-sided KS test for \eqn{H_0: X \sim t_{\nu}} with \eqn{c},
#'     scale \eqn{s}, and degrees of freedom \eqn{\nu}. If parameters are not
#'     specified, the MLE given the data will be used (see
#'     \code{\link[MASS]{fitdistr}}).
#' 
#' For estimated parameters of the t-distribution the p-values are incorrect and
#'     should be adjusted. See \code{\link[stats]{ks.test}} and the references
#'     therein (Durbin (1973)).  As a more practical approach consider
#'     bootstrapping and estimating the p-value empirically.
#' 
#' @param x a numeric vector of data values.
#' @param param 3-dimensional named vector \code{('location', 'scale', 'df')} 
#' which parametrizes the student t distribution. Default: \code{param = NULL}, 
#' in which case it will be estimated from \code{x}.
#' @return 
#' A list of class \code{"htest"} containing:
#' \item{statistic}{the value of the Kolomogorv-Smirnov statistic.}
#' \item{p.value }{the p-value for the test.} 
#' \item{alternative }{a character string describing the alternative hypothesis.} 
#' \item{method}{the character string "One-sample Kolmogorov-Smirnov test
#'               student-t" plus rounded parameter values.} 
#' \item{data.name}{a character string giving the name(s) of the data.}
#' 
#' @export
#' @seealso \code{\link[MASS]{fitdistr}}, \code{\link[stats]{ks.test}}
#' @keywords htest
#' @examples
#' set.seed(1021)
#' beta.true <- c(location = 0, scale = 1, df = 4)
#' xx <- rt(n = 1000, df = beta.true['df'])
#' ks.test.t(xx)
#' ks.test.t(xx, beta.true)
#' 
ks.test.t <- function(x, param = NULL) {
  
  if (is.null(param)) {
    param <- suppressWarnings(fitdistr(x, "t")$est)
    names(param) <- c("location", "scale", "df")
  }
  stopifnot(length(param) == 3,
            is.numeric(param),
            all(match(names(param), c("location", "scale", "df"))),
            param['scale'] > 0, # standard deviation must be positive
            param['df'] > 2) # degrees of freedom must be 
  
  standardized.data <- (x - param['location']) / param['scale']
  test.result <- suppressWarnings(ks.test(standardized.data, "pt",
                                          df = param['df']))
  test.result$data.name <- deparse(substitute(x))
  test.result$method <- 
    paste0("One-sample Kolmogorov-Smirnov test \n student-t with \n", 
           " df=", round(param['df'], 2), 
           ", location=", round(param['location'], 2), 
           ", scale=", round(param['scale'], 2))
  return(test.result)
} 

#' Predefined functions
#'
#' \code{getfunctions} lists the predefined functions in \code{quickpsy}.
#' @seealso \code{\link{cum_normal_fun}},
#' \code{\link{logistic_fun}},
#' \code{\link{weibull_fun}}
#' @export
get_functions <- function() {
  list(cum_normal_fun = cum_normal_fun,
       logistic_fun = logistic_fun,
       weibull_fun = weibull_fun)
}

#' @title Means and Confidence Intervals
#'
#' @description A function for calculating and formatting means and 
#' confidence interval.
#'
#' @details
#' Given a numeric vector, \code{mean_ci} will return a vector with the mean,
#' LCL, and UCL.  Using \code{frmtci} will be helpfull for reporting the results
#' in print.
#'
#' The \code{transform} arguement allows the use to transform the data.  A
#' common occurance of using \code{mean_ci(log(x), transform = exp)} will return
#' the geometric mean and confidence interval for x.
#'
#' @param x a numeric vector
#' @param na_rm if true, omit NA values
#' @param transform function transform to the mean and the confidence limits.
#' See Details.
#' @param alpha defaults to \code{getOption('qwraps2_alpha', 0.05)}.  The
#' symmetric 100(1-alpha)\% CI will be determined.
#' @param qdist defaults to \code{qnorm}.  use \code{qt} for a Student t
#' intervals.
#' @param ... args passed to \code{qdist}
#'
#' @return a vector with the mean, lower confidence limit (LCL), and the upper
#' confidence limit (UCL).
#'
#' @examples
#' # using the standard normal for the CI
#' mean_ci(mtcars$mpg)
#' 
#' # print it nicely
#' qwraps2::frmtci(mean_ci(mtcars$mpg))
#' qwraps2::frmtci(mean_ci(mtcars$mpg), show_level = TRUE)
#' qwraps2::frmtci(mean_ci(mtcars$mpg, alpha = 0.01), show_level = TRUE)
#' 
#' # Compare to the ci that comes form t.test
#' t.test(mtcars$mpg)
#' t.test(mtcars$mpg)$conf.int
#' mean_ci(mtcars$mpg, qdist = stats::qt, df = 31)
#' 
#' # geometric version
#' mean_ci(log(mtcars$mpg), transform = exp, qdist = stats::qt, df = 31)
#'
#' @export   
mean_ci <- function(x, 
                    na_rm = FALSE, 
                    transform,
                    alpha = getOption("qwraps2_alpha", 0.05), 
                    qdist = stats::qnorm, 
                    ...) { 

  m <- mean(x, na.rm = na_rm)
  s <- stats::sd(x,   na.rm = na_rm)
  n <- if (na_rm) { sum(!is.na(x)) } else { length(x) }

  qd <- match.fun(qdist)

  scores <- qd(c(alpha/2, 1 - alpha/2), ...)

  out <- c("mean" = m, "lcl"  = m + scores[1] * s / sqrt(n),  "ucl"  = m + scores[2] * s / sqrt(n))

  if (!missing(transform)) {
    trans <- match.fun(transform)
    out <- trans(out)
  }

  attr(out, 'alpha') <- alpha
  class(out) <- c("qwraps2_mean_ci", class(out)) 
  out 
}


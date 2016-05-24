#' Confidence Intervals for Proportions
#' 
#' Alternatives to \code{prop.test} and \code{binom.test}.
#' 
#' \code{wald.ci} produces Wald confidence intervals. \code{wilson.ci}
#' produces Wilson confidence intervals (also called ``plus-4'' confidence
#' intervals) which are Wald intervals computed from data formed by adding 2
#' successes and 2 failures.  The Wilson confidence intervals have better
#' coverage rates for small samples.
#' 
#' @aliases wilson.ci wald.ci
#' @param x number of 'successes'
#' @param n number of trials
#' @param conf.level confidence level
#' @return Lower and upper bounds of a two-sided confidence interval.
#' @author Randall Pruim
#' @references A. Agresti and B. A. Coull, Approximate is better then `exact'
#' for interval estimation of binomial proportions, \emph{American
#' Statistician} 52 (1998), 119--126.
#' @export
#' @examples
#' 
#' prop.test(12,30)
#' prop.test(12,30, correct=FALSE)
#' wald.ci(12,30)
#' wilson.ci(12,30)
#' wald.ci(12+2,30+4)
#' 
wilson.ci <-
function (x, n = 100, conf.level = 0.95) 
{
    alpha = 1 - conf.level
    p = (x + 2)/(n + 4)
    zstar <- -qnorm(alpha/2)
    interval <- p + c(-1, 1) * zstar * sqrt(p * (1 - p)/(n+4))
    attr(interval, "conf.level") <- conf.level
    return(interval)
}

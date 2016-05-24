#' @encoding UTF-8
#' @title An S4 Class to Confidence Intervals
#'
#' @slot lower Lower bound of interval.
#' @slot mean Estimated mean.
#' @slot upper Upper bound of interval.
#' @slot stderr Standard Error of the mean.
#'
#' @export
`ci`<-setClass(Class = "ci",
               slots = list(lower = "numeric",
                            mean = "numeric",
                            upper = "numeric",
                            stderr="numeric"))



#' @encoding UTF-8
#' @title Confidence Intervals
#' @description Calculates the confidence intervals for a vector of data values.
#' @keywords Exploratory
#' @param x A vector of data values.
#' @param level The confidence level. Default is \code{0.95}.
#' @param alpha The significance level. Default is \code{1-level}. If alpha equals 0.05, then your confidence level is 0.95.
#' @param na.rm A logical value, default is \code{FALSE}
#' @param \dots Additional arguements (currently ignored)
#'
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}.
#' @return
#' \item{CI lower}{}
#' \item{Est. Mean}{Mean of data.}
#' \item{CI upper}{Upper bound of interval.}
#' \item{Std. Error}{Standard Error of the mean.}
#'
#' @examples
#' x <- c(1, 2.3, 2, 3, 4, 8, 12, 43, -1,-4)
#'
#' ci(x, level=.90)
#' @export
#'
`ci` <- function(x, level=0.95, alpha=1-level,na.rm=FALSE,...){
  estimate <- mean(x, na.rm = na.rm);
  stderr <- stats::sd(x, na.rm=na.rm)/sqrt(length(x));
  ci.low <- estimate + stats::qt(alpha/2,length(x)-1)*stderr;
  ci.high <- estimate - stats::qt(alpha/2,length(x)-1)*stderr;
  retval <- new("ci",
                lower=ci.low,
                mean=estimate,
                upper=ci.high,
                stderr=stderr
  );
  retval;
}
NULL


setMethod("show", signature(object="ci"),
          definition=function(object) {
            retval <- c(
              "CI Lower"=object@lower,
              "Est. Mean"=object@mean,
              "CI Upper"=object@upper,
              "Std. Error"=object@stderr
  )
print(retval, digits = max(3, getOption("digits") - 3))
}
)

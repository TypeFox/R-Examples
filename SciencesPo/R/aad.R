#' @encoding UTF-8
#' @title An S4 Class to Average Absolute Deviation
#'
#' @slot estimate Estimated value.
#'
#' @export
setClass(Class = "aad",
               slots = list(estimate = "numeric"))


#' @encoding UTF-8
#' @title Average Absolute Deviation
#'
#' @description Calculates the average (mean) absolute deviation from the sample mean.
#' @param x	A numeric vector containing the observations.
#' @param na.rm A logical value for \code{na.rm}, default is \code{na.rm=TRUE}.
#' @param \dots Additional arguements (currently ignored)

#' @details The statistical literature has not yet adopted a standard notation for  the "Mean Absolute Deviation" and the "Median Absolute Deviation". As a result, both statistics have been denoted as "MAD", which may lead to confusion once they may produce different values.
#' The R \code{\link[stats]{mad}} by default computes the "Median Absolute Deviation"; to obtain the "Mean Absolute Deviation" one has to use \code{mad(x, constant = 1)}.
#' Thus, the function \code{\link[SciencesPo]{aad}} will calculate the "Mean Absolute Deviation"--or "Average Deviation (AD)" as proposed by Garrett, who defines it as "the mean of the deviation of all the separate scores in the series taken from their mean (occasionally from the median or mode)", (1971, p. 481).
#'
#' @references
#' Garrett, Henry (1982) \emph{Statistics in Psychology and Education}. 6th, Paragon.
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}
#' @keywords Exploratory
#' @seealso \code{\link[stats]{mad}}
#' @examples
#' x <- c(15, 10, 6, 8, 11)
#' aad(x)
#'
#' @export
`aad`<-function(x, na.rm = TRUE, ...){
  if (!is(x, "numeric") & !is(x, "integer")) {
    stop("\"x\" must be numeric")
  }
  if (!is(na.rm, "logical") | length(na.rm) != 1) {
    stop("\"na.rm\" must be a single logical value")
  }
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  ans <- mean(abs(x - mean(x)))
  retval <- new("aad", estimate=ans);
  retval;
}## -- end of aad
NULL


setMethod("show", signature(object="aad"),
          definition=function(object) {
            retval <- c(object@estimate)
            print(retval, digits = max(3, getOption("digits") - 3))
})

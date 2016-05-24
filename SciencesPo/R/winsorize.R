#' @encoding UTF-8
#' @title Winsorized Mean
#'
#' @description Compute the winsorized mean, which consists of recoding the top k values in a vector.
#'
#' @param x The vector to be winsorized
#' @param k An integer for the quantity of outlier elements that to be replaced in the calculation process
#' @param na.rm a logical value for \code{na.rm}, default is \code{na.rm=TRUE}.
#'
#' @details Winsorizing a vector will produce different results than trimming it. While by trimming a vector causes extreme values to be discarded, by winsorizing it in the other hand, causes extreme values to be replaced by certain percentiles.
#'
#' @return An object of the same type as \code{x}
#'
#' @references  Dixon, W. J., and Yuen, K. K. (1999) Trimming and winsorization: A review. \emph{The American Statistician,} \bold{53(3),} 267--269.
#' @references  Dixon, W. J., and Yuen, K. K. (1960) Simplified Estimation from Censored Normal Samples, \emph{The Annals of Mathematical Statistics,} \bold{31,} 385--391.
#'  @references  Wilcox, R. R. (2012) \emph{Introduction to robust estimation and hypothesis testing.} Academic Press, 30-32. Statistics Canada (2010) \emph{Survey Methods and Practices.}
#'
#'  @note One may want to winsorize estimators, however, winsorization tends to be used for one-variable situations.
#'
#'  @author Daniel Marcelino, \email{dmarcelino@@live.com}
#'
#'
#'  @examples
#' set.seed(51)  # for reproducibility
#' x <- rnorm(50)
#' ## introduce outliers
#' x[1] <- x[1] * 10
#' # Compare to mean:
#'  mean(x)
#'  winsorize(x)
#' @keywords Descriptive
#'
#' @export
#'
`winsorize` <-
  function (x, k = 1, na.rm=TRUE) {
    if (any(is.na <- is.na(x))) {
      if (na.rm)
        x <- x[!is.na]
      else return(NA)
    }
    n <- length(x)
    if (!(k %in% (0:n)))
      stop("'k' should be > 0 and less than half the number of non-missing observations.")
    else {
      x <- sort(x)
      x[1:k] <- x[k+1] # Here I solve the lower values
      x[(n-k+1):n] <- x[n-k] #Then I go over the higher ones
      return(mean(x))
    }
  }### end -- winsorize function
NULL

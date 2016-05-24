################################################################################
#' Validates if frequencies are Fourier frequencies from
#' \eqn{[0,\pi]}{[0,pi]}.
#'
#' Validation of the parameter \code{freq} is perfomed in six steps:
#' \enumerate{
#' \item Throw an error if parameter is not a vector or not numeric.
#' \item Transform each element \eqn{\omega}{w} of the vector to
#'            \eqn{[0,2\pi)}{[0,2pi)}, by replacing it with
#'            \eqn{\omega \, \mbox{mod} \, 2\pi}{w mod 2pi}.
#' \item Check whether all elements \eqn{\omega}{w} of the vector are
#'            Fourier frequency \eqn{2 \pi j / T}{2 pi j / T}, \eqn{j \in Z}.
#'            If this is not
#'            the case issue a warning and round each frequency to the next
#'            Fourier frequency of the mentioned type; the smaller one, if
#'            there are two.
#' \item Transform each element \eqn{\omega}{w} with
#'            \eqn{\pi < \omega < 2\pi}{pi < w < 2pi} of the vector to
#'            \eqn{[0,\pi]}{[0,pi]}, by replacing it with
#'            \eqn{2\pi - \omega}{2pi - w}.
#' \item Check for doubles and remove all but the first appearance.
#' \item Sort in ascending order.
#' }
#' Any subset of the six steps can be chosen, but 1 should almost always be
#' among the steps to be performed.
#'
#' @name frequenciesValidator
#' @export
#'
#' @keywords Validator-functions
#'
#' @param freq the vector of frequencies to be validated.
#' @param N    the base of the Fourier frequencies against which the values in
#'             \code{freq} will be compared.
#' @param steps a vector containing a subset of {1,2,3,4,5,6}, indicating
#'              which of the steps are to be performed.
#'
#' @return Returns a vector of Fourier frequencies that is yield by the
#'         transformations described above.
#'
#' @examples
#' freq <- 2*pi*c(3,2,5,8,9)/10
#'
#' res <- frequenciesValidator(freq, N=10, steps=1:3)
#' res * 10 / (2*pi) # Returns: [1] 3 2 5 8 9
#'
#' res <- frequenciesValidator(freq, N=10, steps=1:4)
#' res * 10 / (2*pi) # Returns: [1] 3 2 5 2 1
#'
#' res <- frequenciesValidator(freq, N=10, steps=1:5)
#' res * 10 / (2*pi) # Returns: [1] 3 2 5 1
#'
#' res <- frequenciesValidator(freq, N=10, steps=1:6)
#' res * 10 / (2*pi) # Returns: [1] 1 2 3 5
################################################################################

frequenciesValidator <- function(freq, N, steps=1:6) {

  f <- function(s) {Vectorize(toString)(s)}
  steps <- match.arg(f(steps), choices = f(1:6), several.ok=TRUE)

  if (is.element("1", steps)) {
    if (!(is.vector(freq)  && is.numeric(freq))) {
      stop("'frequencies' needs to be specified as a vector of real numbers")
    }
  }

  if (is.element("2", steps)) {
    # Transform all frequencies to [0,2pi)
    freq <- freq %% (2*pi)
  }

  setequal.approx <- function(x,y) {
    X <- round(x, .Machine$double.exponent-2)
    Y <- round(y, .Machine$double.exponent-2)
    return(setequal(X,Y))
  }
  if (is.element("3", steps)) {
    # Transform all frequencies to Fourier frequencies
    # and issue a warning if that is necessary.
    if (!setequal.approx(round(N*freq/(2*pi)), N*freq/(2*pi))) {
      warning("Some 'frequencies' were rounded to the next Fourier frequencies")
      freq <- (2*pi/N) * round(N/(2*pi)*freq)
    }
  }

  if (is.element("4", steps)) {
    # Transform all frequencies from (pi,2pi) to [0,pi]
    pos <- which(freq > pi & freq <= 2*pi)
    freq[pos] <- 2*pi - freq[pos]
  }

  if (is.element("5", steps)) {
    pos.unique <- !duplicated(round(freq, .Machine$double.exponent-2))
    freq <- freq[pos.unique]
  }

  if (is.element("6", steps)) {
    freq <- sort(freq)
  }

  return(freq)
}


################################################################################
#' Positions of elements which are closest to some reference elements.
#'
#' For two vectors \code{X} and \code{Y} a vector of indices \code{I} is returned,
#' such that \code{length(Y)} and \code{length(I)} coincide and \code{X[I[j]]}
#' is an element of \code{X} which has minimal distance to \code{Y[j]}, for all
#' \code{j=1,\ldots,length(Y)}.
#' In case that there are multiple elements with minimal distance, the smallest
#' index (the index of the first element with minimal distance) is returned.
#'
#' @name   closest.pos
#' @rdname closest.pos
#' @export
#' @param X Vector of elements among which to find the closest one for each
#'          element in \code{Y}.
#' @param Y Vector of elements for which to find the clostest element in \code{X}.
#'
#' @return Returns a vector of same length as \code{X}, with indices indicating
#'         which element in \code{Y} is closest.
#' @examples
#' X1 <- c(1,2,3)
#' closest.pos(X1, 1.7)
#' closest.pos(X1, c(1.3,2.2))
#'
#' X2 <- c(2,1,3)
#' closest.pos(X2, 1.5)

closest.pos <- function(X, Y) {
  g.ptw <- function(u) {which.min(abs(X-u))}
  return(Vectorize(g.ptw)(Y))
}



################################################################################
#' Validates if \code{Y} is of an appropriate type and converts to a numeric.
#'
#' Checks whether \code{Y} is either
#' \itemize{
#'   \item \code{numeric},
#'   \item a \code{ts} object, or
#'   \item a \code{zoo} object.
#' }
#' If not, an error is returned. If it is one of the three the data is returned
#' as a numeric.
#'
#' @name timeSeriesValidator
#' @export
#' 
#' @importFrom zoo is.zoo coredata
#' @importFrom stats is.ts
#' @importFrom utils head
#'
#' @keywords Validator-functions
#'
#' @param Y the time series to be validated.
#'
#' @return Returns the time series as a matrix.
#'
#' @examples
#' Y <- timeSeriesValidator(sp500)
#' Y <- timeSeriesValidator(wheatprices)
#' Y <- timeSeriesValidator(rnorm(10))
#' \dontrun{Y <- timeSeriesValidator("Not a valid input")}
################################################################################

timeSeriesValidator <- function(Y) {
  if (!(is.matrix(Y)  && is.numeric(Y)) & !(is.vector(Y)  && is.numeric(Y)) & !is.ts(Y) & !is.zoo(Y)) {
    stop("'Y' needs to be specified as a vector or matrix of real numbers, a ts or a zoo object")
  }

  if (is.ts(Y)) {
    Y <- head(Y)
  }
  
  if (is.vector(Y)) {
    Y <- matrix(Y, ncol=1)
  }

  if (is.zoo(Y)) {
    Y <- as.matrix(Y)
  }
  return(Y)
}



################################################################################
#' Validates if \code{Y} is of an appropriate type for a time series and
#' returns the length of the time series.
#'
#' Runs \code{\link{timeSeriesValidator}} and returns the number of rows of the
#' returned matrix.
#'
#' @name lenTS
#' @export
#' @importFrom zoo is.zoo coredata
#'
#' @keywords Validator-functions
#'
#' @param Y the time series to be validated and of which the length is to
#' 					be returned.
#'
#' @return Returns the length of the time series after validating it's valid.
#'
#' @examples
#' Y <- lenTS(sp500)
#' Y <- lenTS(wheatprices)
#' Y <- lenTS(rnorm(10))
#' \dontrun{Y <- lenTS("Not a valid input")}
################################################################################

lenTS <- function(Y) {
  
  Y <- timeSeriesValidator(Y)
  
  return(dim(Y)[1])
}


################################################################################
#' Checks whether \code{x} contains integer numbers.
#'
#' Borrowed from the example in \code{\link[base]{integer}}.
#'
#' @name is.wholenumber
#'
#' @keywords internals
#'
#' @param x   a vector to be checked for integers
#' @param tol an optional parameter specifying to which precision the check is
#'            to be performed.
#'
#' @return Returns a vector of logicals with the same length as \code{x}; each
#'         element \code{i} is \code{TRUE} iff \code{x[i]} is an integer.
#'
#' @examples
#' \dontrun{
#' is.wholenumber(1) # is TRUE
#' (x <- seq(1, 5, by = 0.5) )
#' is.wholenumber( x ) #-->  TRUE FALSE TRUE ...
#' }
################################################################################
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  return(abs(x - round(x)) < tol)
}

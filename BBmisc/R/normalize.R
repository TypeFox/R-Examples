#' @title Normalizes numeric data to a given scale.
#'
#' @description
#' Currently implemented for numeric vectors, numeric matrices and data.frame.
#' For matrixes one can operate on rows or columns
#' For data.frames, only the numeric columns are touched, all others are left unchanged.
#' For constant vectors / rows / columns most methods fail, special behaviour for this
#' case is implemented.
#'
#' The method also handles NAs in in \code{x} and leaves them untouched.
#'
#' @param x [\code{numeric} | \code{matrix} | \code{data.frame}]\cr
#'   Input vector.
#' @param method [\code{character(1)}]\cr
#'   Normalizing method. Available are:\cr
#'   \dQuote{center}: Subtract mean.\cr
#'   \dQuote{scale}: Divide by standard deviation.\cr
#'   \dQuote{standardize}: Center and scale.\cr
#'   \dQuote{range}: Scale to a given range.\cr
#' @param range [\code{numeric(2)}]\cr
#'   Range for method \dQuote{range}.
#'   Default is \code{c(0,1)}.
#' @param margin [\code{integer(1)}]\cr
#'   1 = rows, 2 = cols.
#'   Same is in \code{\link{apply}}
#'   Default is 1.
#' @param on.constant [\code{character(1)}]\cr
#'   How should constant vectors be treated? Only used, of \dQuote{method != center},
#'   since this methods does not fail for constant vectors. Possible actions are:\cr
#'   \dQuote{quiet}: Depending on the method, treat them quietly:\cr
#'     \dQuote{scale}: No division by standard deviation is done, input values.
#'        will be returned untouched.\cr
#'     \dQuote{standardize}: Only the mean is subtracted, no division is done.\cr
#'     \dQuote{range}: All values are mapped to the mean of the given range.\cr
#'   \dQuote{warn}: Same behaviour as \dQuote{quiet}, but print a warning message.\cr
#'   \dQuote{stop}: Stop with an error.\cr
#' @return [\code{numeric} | \code{matrix} | \code{data.frame}].
#' @seealso \code{\link{scale}}
#' @export
normalize = function(x, method = "standardize", range = c(0, 1), margin = 1L, on.constant = "quiet") {
  assertChoice(method, c("range", "standardize", "center", "scale"))
  assertNumeric(range, len = 2L, any.missing = FALSE)
  assertChoice(on.constant, c("quiet", "warn", "stop"))
  UseMethod("normalize")
}

#' @export
normalize.numeric = function(x, method = "standardize", range = c(0, 1), margin = 1L, on.constant = "quiet") {
  y = normalize2(x, method, range, on.constant = on.constant)
  # scale call below returns matrices
  if (is.matrix(y))
    y = y[,1L]
  return(y)
}

#' @export
normalize.matrix = function(x, method = "standardize", range = c(0, 1), margin = 1L, on.constant = "quiet") {
  x = apply(x, margin, normalize2, method = method, range = range, on.constant = on.constant)
  if (margin == 1L)
    x = t(x)
  return(x)
}

#' @export
normalize.data.frame = function(x, method = "standardize", range = c(0, 1), margin = 1L, on.constant = "quiet") {
  isnum = sapply(x, is.numeric)
  if (any(isnum))
    x = as.data.frame(lapply(x[, isnum, drop = FALSE], normalize2, method = method,
      range = range, on.constant = on.constant))
  return(x)
}

normalize2 = function(x, method, range, on.constant) {
  # is x a constant vector?
  if (length(unique(x[!is.na(x)])) == 1L) {
    switch(on.constant,
      warn = warning("Constant vector in normalization."),
      stop = stop("Constant vector in normalization."))
    switch(method,
      center = scale(x, center = TRUE, scale = FALSE),
      range = ifelse(is.na(x), NA, mean(range)),
      standardize = scale(x, center = TRUE, scale = FALSE),
      scale = x
    )
  } else {
    switch(method,
      range = (x - min(x, na.rm = TRUE)) / diff(range(x, na.rm = TRUE)) * diff(range) + range[1L],
      standardize = scale(x, center = TRUE, scale = TRUE),
      center = scale(x, center = TRUE, scale = FALSE),
      scale = scale(x, center = FALSE, scale = sd(x, na.rm = TRUE))
    )
  }
}



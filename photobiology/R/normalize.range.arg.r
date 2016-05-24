#' Normalize a range argument into a true numeric range
#'
#' Several functions in this package and the suite accept a range argument
#' with a flexible syntax. To ensure that all functions and methods behave
#' in the same way this code has been factored out into a separate function.
#'
#' @param range a numeric vector of length two, or any other object for which
#'   function range() will return two
#' @param default.range a numeric vector of length two, missing values are
#'   not allowed, but Inf amd -Inf are.
#'
#' @return a numeric vector of length two, guaranteed not to have missing
#'   values.
#'
#' @details The \code{range} argument can contain NAs which are replaced by
#'   the value at the same position in \code{default.range}. In addition
#'   a NULL argument for \code{range} is converted into \code{default.range}.
#'   The \code{default.range} is also the limit to which the returned value
#'   is trimmed. The idea is that the value supplied as default is the whole
#'   valid range, and as we use range only for wavelength, the default is
#'   0 to Inf.
#'
#' @family auxiliary functions
#'
#' @export
#' @examples
#' normalize_range_arg(sun.spct)
#' normalize_range_arg(c(NA, 500))
#' normalize_range_arg(c(NA, 500), c(100, Inf))
#' normalize_range_arg(c(-100, 500), c(-Inf, Inf))
#'
normalize_range_arg <- function(range, default.range = c(0, Inf)) {
  stopifnot(is.numeric(default.range) &&
              length(default.range) == 2 &&
                 default.range[1] < default.range[2])
  if (is.null(range) || all(is.na(range))) {
    return(default.range)
  }
  if (!is.numeric(range) || (is.numeric(range) && length(range) != 2)) {
    range <- range(range, na.rm = TRUE)
  }
  stopifnot(is.numeric(range) && length(range) == 2)

  if (is.na(range[1]))
    range[1] <- default.range[1]
  else
      range[1] <- max(range[1], default.range[1])

  if (is.na(range[2]))
    range[2] <- default.range[2]
  else
    range[2] <- min(range[2], default.range[2])
  # simplest here as all NAs have been replaced above
  stopifnot(range[1] < range[2])
  range
}

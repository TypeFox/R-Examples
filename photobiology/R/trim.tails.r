#' Trim (or expand) head and/or tail
#'
#' Trimming of tails of a spectrum based on wavelength limits, interpolating the
#' values at the boundaries.Trimming is needed for example to remove short
#' wavelength noise when the measured spectrum extends beyond the known emission
#' spectrum of the measured light source. Occasionally one may want also to
#' expand the wavelength range.
#'
#' @param x numeric array
#' @param y numeric array
#' @param low.limit smallest x-value to be kept (defaults to smallest
#'   x-value in input)
#' @param high.limit largest x-value to be kept (defaults to largets
#'   x-value in input)
#' @param use.hinges logical, if TRUE (the default)
#' @param fill if fill==NULL then tails are deleted, otherwise tails of y
#'   are filled with the value of fill
#'
#' @return a data.frame with variables \code{x} and \code{y}
#'
#' @note When expanding an spectrum, if fill==NULL, then expansion is not
#'   performed with a warning.
#'
#' @family trim functions
#' @export
#' @examples
#' head(sun.data)
#' head(with(sun.data,
#'      trim_tails(w.length, s.e.irrad, low.limit=300)))
#' head(with(sun.data,
#'      trim_tails(w.length, s.e.irrad, low.limit=300, fill=NULL)))
#'
trim_tails <- function(x, y,
                       low.limit=min(x), high.limit=max(x),
                       use.hinges=TRUE, fill=NULL)
{
#  if (!check_spectrum(w.length, s.irrad)) return(NA)

  # check whether we should expand the low end
  low.end <- min(x)
  if (low.end > low.limit) {
    if (!is.null(fill)) {
      # expand short tail
      low.tail.length <- low.end - low.limit
      low.tail <- seq(from = low.limit, to = low.end - 1, length.out = low.tail.length)
      y <- c(rep(fill, length.out = low.tail.length), y)
      x <- c(low.tail, x)
      low.end <- low.limit
    } else {
      warning("Ignoring low.limit as it is too low.")
      low.limit <- low.end
    }
  }

  # check whether we should expand the high end
  high.end <- max(x)
  if (high.end < high.limit) {
    if (!is.null(fill)) {
      # expand short tail
      high.tail.length <- high.limit - high.end
      high.tail <- seq(from = high.end + 1, to = high.limit, length.out = high.tail.length)
      y <- c(y, rep(fill, length.out = high.tail.length))
      x <- c(x, high.tail)
      high.end <- high.limit
    } else {
      warning("Ignoring high.limit as it is too high.")
      high.limit <- high.end
    }
  }

  # insert hinges
  if (use.hinges) {
    new.data <- insert_hinges(x, y, c(low.limit, high.limit))
  } else {
    new.data <- dplyr::data_frame(x = x, y = y)
  }

  trimmed.selector <- with(new.data, (x >= low.limit) & (x <= high.limit))

  if (is.null(fill)) {
    new.data[trimmed.selector, ]
  }
  else {
    new.data[!trimmed.selector, ] <- fill
    new.data
  }
}

#' Insert wavelength values into spectral data.
#'
#' Inserting wavelengths values immediately before and after a discontinuity in
#' the SWF, greatly reduces the errors caused by interpolating the weighted
#' irradiance during integration of the effective spectral irradiance. This is
#' specially true when data have a large wavelength step size.
#'
#' @param x numeric array (sorted in increasing order)
#' @param y numeric array
#' @param h a numeric array giving the wavelengths at which the y values
#'   should be inserted by interpolation, no interpolation is indicated by an
#'   empty array (numeric(0))
#'
#' @return a data.frame with variables \code{x} and \code{y}
#'
#' @export
#' @examples
#' with(sun.data,
#'     insert_hinges(w.length, s.e.irrad,
#'        c(399.99, 400.00, 699.99, 700.00)))
#' with(sun.data,
#'     insert_hinges(w.length, s.e.irrad,
#'       c(100, 399.50, 399.99, 400.00, 699.99, 700.00, 799.99, 1000)))
#'
#' @note Insertion is a costly operation but I have tried to optimize this
#' function as much as possible by avoiding loops. Earlier this function was
#' implemented in C++, but a bug was discovered and I have now rewritten it
#' using R.
#'
insert_hinges <- function(x, y, h) {
  # sanitize 'hinges'
  h <- setdiff(h, x)
  h <- h[h > x[1] & h < x[length(x)]]
  h <- sort(h)
  # save lengths
  j <- length(h)
  k <- length(x)
  # allocate vectors for the output
  x.out <- numeric(j + k)
  y.out <- numeric(j + k)
  # compute the new idx values after shift caused by insertions
  idxs.in <- findInterval(h, x)
  idxs.out <- idxs.in + 1:j
  idxs.diff <- diff(c(0,idxs.in,k))
  idxs.map <- 1:k + rep(0:j, idxs.diff)
  # we copy everything that does not require interpolation
  x.out[idxs.map] <- x
  y.out[idxs.map] <- y
  x.out[idxs.out] <- h
  # we use recycling to interpolate all values and insert them into the gaps
  if (is.numeric(y)) {
    y.out[idxs.out] <- y[idxs.in + 1] -
      (x[idxs.in + 1] - x.out[idxs.out]) /
      (x[idxs.in + 1] - x[idxs.in]) *
      (y[idxs.in + 1] - y[idxs.in])
  } else {
    y.out[idxs.out] <- NA
  }
  data.frame(x = x.out, y = y.out)
}

#' Insert wavelength values into spectral data.
#'
#' Inserting wavelengths values immediately before and after a discontinuity in
#' the SWF, greatly reduces the errors caused by interpolating the weighted
#' irradiance during integration of the effective spectral irradiance. This is
#' specially true when data have a large wavelength step size.
#'
#' @param x numeric array (sorted in increasing order)
#' @param y numeric array
#' @param h a numeric array giving the wavelengths at which the y values
#'   should be inserted by interpolation, no interpolation is indicated by an
#'   empty array (numeric(0))
#'
#' @return a list with variables \code{x} and \code{y}
#'
#' @keywords internal
#'
l_insert_hinges <- function(x, y, h) {
  # sanitize 'hinges'
  h <- setdiff(h, x)
  h <- h[h > x[1] & h < x[length(x)]]
  h <- sort(h)
  # save lengths
  j <- length(h)
  k <- length(x)
  # allocate vectors for the output
  x.out <- numeric(j + k)
  y.out <- numeric(j + k)
  # compute the new idx values after shift caused by insertions
  idxs.in <- findInterval(h, x)
  idxs.out <- idxs.in + 1:j
  idxs.diff <- diff(c(0,idxs.in,k))
  idxs.map <- 1:k + rep(0:j, idxs.diff)
  # we copy everything that does not require interpolation
  x.out[idxs.map] <- x
  y.out[idxs.map] <- y
  x.out[idxs.out] <- h
  # we use recycling to interpolate all values and insert them into the gaps
  if (is.numeric(y)) {
    y.out[idxs.out] <- y[idxs.in + 1] -
      (x[idxs.in + 1] - x.out[idxs.out]) /
      (x[idxs.in + 1] - x[idxs.in]) *
      (y[idxs.in + 1] - y[idxs.in])
  } else {
    y.out[idxs.out] <- NA
  }
  list(x = x.out, y = y.out)
}

#' Insert wavelength values into spectral data.
#'
#' Inserting wavelengths values immediately before and after a discontinuity in
#' the SWF, greatly reduces the errors caused by interpolating the weighted
#' irradiance during integration of the effective spectral irradiance. This is
#' specially true when data have a large wavelength step size. This fucntion
#' differs from \code{insert_hinges()} in that it returns a vector of \code{y}
#' values instead of a \code{data_frame}.
#'
#' @param x numeric array (sorted in increasing order)
#' @param y numeric array
#' @param h a numeric array giving the wavelengths at which the y values
#'   should be inserted by interpolation, no interpolation is indicated by an
#'   empty array (numeric(0))
#'
#' @return a numeric vector of \code{y} values
#'
#' @keywords internal
#'
v_insert_hinges <- function(x, y, h) {
  # sanitize 'hinges'
  h <- setdiff(h, x)
  h <- h[h > x[1] & h < x[length(x)]]
  h <- sort(h)
  # save lengths
  j <- length(h)
  k <- length(x)
  # allocate vectors for the output
  x.out <- numeric(j + k)
  y.out <- numeric(j + k)
  # compute the new idx values after shift caused by insertions
  idxs.in <- findInterval(h, x)
  idxs.out <- idxs.in + 1:j
  idxs.diff <- diff(c(0,idxs.in,k))
  idxs.map <- 1:k + rep(0:j, idxs.diff)
  # we copy everything that does not require interpolation
  x.out[idxs.map] <- x
  y.out[idxs.map] <- y
  x.out[idxs.out] <- h
  # we use recycling to interpolate all values and insert them into the gaps
  if (is.numeric(y)) {
    y.out[idxs.out] <- y[idxs.in + 1] -
      (x[idxs.in + 1] - x.out[idxs.out]) /
      (x[idxs.in + 1] - x[idxs.in]) *
      (y[idxs.in + 1] - y[idxs.in])
  } else {
    y.out[idxs.out] <- NA
  }
  y.out
}



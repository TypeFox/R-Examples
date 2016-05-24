#' Trim (or expand) head and/or tail
#'
#' Trimming of waveband boundaries can be required needed when the spectral data
#' does not cover the whole waveband.
#'
#' @param w.band an object of class "waveband" or a list of such objects
#' @param range a numeric vector of length two, or any other object for which
#'   function range() will return a numeric vector of two wavelengths (nm)
#' @param low.limit shortest wavelength to be kept (defaults to 0 nm)
#' @param high.limit longest wavelength to be kept (defaults to Inf nm)
#' @param trim logical (default is TRUE which trims the wavebands at the
#'   boundary, while FALSE discards wavebands that are partly off-boundary).
#' @param use.hinges logical, if TRUE (the default) hinges are inserted when
#'   trimming.
#'
#' @return a waveband object or a list of waveband objects trimmed or filtered
#'   depending on whether a single waveband object or a list of waveband
#'   objects was supplied as argument to formal parameter \code{w.band}.
#'
#' @family trim functions
#' @export
#' @examples
#' VIS <- waveband(c(380, 760)) # nanometres
#'
#' trim_waveband(VIS, c(400,700))
#' trim_waveband(VIS, low.limit = 400)
#' trim_waveband(VIS, high.limit = 700)
#'
trim_waveband <-
  function(w.band,
           range = NULL,
           low.limit = 0, high.limit = Inf,
           trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = TRUE)
  {
    input.was.waveband <- is.waveband(w.band)
    if (input.was.waveband) {
      w.band <- list(w.band)
    }

    range <- normalize_range_arg(range)
    low.limit <- range[1]
    high.limit <- range[2]

    w.band.out <- list()
    i <- 0
    for (wb in w.band) {
      if (min(wb) >= high.limit || max(wb) <= low.limit) {
        next
      }
      if (min(wb) >= (low.limit - 5e-12) && max(wb) <= (high.limit + 5e-12)) {
        i <- i + 1L
        w.band.out[i] <- list(wb)
      } else if (trim) {
        trimmed.wb <- wb
        trimmed.high <- trimmed.low <- FALSE
        if (min(wb) < low.limit) {
          trimmed.wb$low <- low.limit
          trimmed.wb$hinges <- unique(sort(c(low.limit - 1e-12, low.limit,
                                             wb$hinges[wb$hinges >= low.limit])))
          trimmed.low <- TRUE
        }
        if (max(wb) > high.limit) {
          trimmed.wb$high <- high.limit
          trimmed.wb$hinges <- unique(sort(c(wb$hinges[wb$hinges <= high.limit],
                                             high.limit - 1e-12, high.limit)))
          trimmed.high <- TRUE
        }
        if (trimmed.low || trimmed.high) {
          trimmed.tag <-  paste("tr", ifelse(trimmed.low, ".lo", ""),
                                ifelse(trimmed.high, ".hi", ""), sep = "")
          trimmed.wb$label <- paste(wb$label, trimmed.tag, sep = " .")
          trimmed.wb$name <- paste(wb$name, trimmed.tag, sep = ".")
          i <- i + 1L
          w.band.out[[i]] <- trimmed.wb
        }
      }
    }
    # w.band.out is always a list, possibly of length zero
    if (input.was.waveband) {
      # we simplify the list if input was a waveband object instead of a list
      if (length(w.band.out) == 1L) {
        w.band.out <- w.band.out[[1]]
      } else {
        w.band.out <- NULL
      }
    }
    w.band.out
  }

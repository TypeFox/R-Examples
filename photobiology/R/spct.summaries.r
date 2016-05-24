# print -------------------------------------------------------------------

#' Print a spectral object
#'
#' Print method for objects of spectral classes.
#'
#' @param x An object of one of the summary classes for spectra
#' @param ... not used in current version
#' @param n	Number of rows to show. If NULL, the default, will print all rows if
#'   less than option dplyr.print_max. Otherwise, will print dplyr.print_min
#' @param width	Width of text output to generate. This defaults to NULL, which
#'   means use getOption("width") and only display the columns that fit on one
#'   screen. You can also set option(dplyr.width = Inf) to override this default
#'   and always print all columns.
#'
#' @return Returns \code{x} invisibly.
#'
#' @export
#'
#' @note At the moment just a modified copy of dplyr:::print.tbl_df.
#'
#' @name print
#'
#' @examples
#'
#' print(sun.spct)
#' print(sun.spct, n = 5)
#'
print.generic_spct <- function(x, ..., n = NULL, width = NULL)
{
  cat("Object: ", class_spct(x)[1], " ", dplyr::dim_desc(x), "\n", sep = "")
  if (nrow(x)) {
    cat("Wavelength (nm): range ",
        paste(signif(range(x), 8), sep = "", collapse = " to "), ", step ",
        paste(unique(signif(stepsize(x), 7)), sep = "", collapse = " to "),
        "\n", sep = "")
  }
  if (!any(is.na(getWhenMeasured(x)))) {
    cat("Measured on: ", as.character(getWhenMeasured(x)), " UTC\n", sep = "")
  }
  if (!any(is.na(getWhereMeasured(x)))) {
    where.measured <- getWhereMeasured(x)
    cat("Measured at: ", where.measured[["lat"]], " N, ",
        where.measured[["lon"]], " E\n", sep = "")
  }
  if (class_spct(x)[1] %in% c("source_spct", "response_spct")) {
    cat("Time unit: ", as.character(getTimeUnit(x, force.duration = TRUE)),
        "\n", sep = "")
  }
  if (is_scaled(x)) {
    scaling <- getScaled(x)[["f"]]
    cat("Rescaled to '", scaling, "' = 1 \n", sep = "")
  }
  if (is_normalized(x)) {
    norm <- getNormalized(x)
    cat("Spectral data normalized to 1 at ", norm,
        ifelse(is.numeric(norm), " nm \n", " \n"), sep = "")
  }
  if (is_effective(x)) {
    BSWF <- getBSWFUsed(x)
    cat("Data weighted using '", BSWF, "' BSWF\n", sep = "")
  }
  cat("\n")
  print(dplyr::trunc_mat(x, n = n, width = width))
  invisible(x)
}

# print method ------------------------------------------------------------

#' @describeIn print
#'
#' @export
#'
print.generic_mspct <- function(x, ..., n = NULL, width = NULL)  {
  cat("Object: ", class(x)[1], " ", dplyr::dim_desc(x), "\n", sep = "")
  member.names <- names(x)
  for (name in member.names) {
    cat("--- Member:", name, "---\n")
    print(x[[name]], n = n, width = width)
  }
  cat("--- END ---")
  invisible(x)
}

# names of all spectral summary classes -----------------------------------

#' Function that returns a vector containing the names of spectral summary
#' classes.
#'
#' @export
#'
#' @return A \code{character} vector of class names.
#'
summary_spct_classes <- function() {
  c("summary_cps_spct",
    "summary_filter_spct", "summary_reflector_spct",
    "summary_source_spct", "summary_object_spct",
    "summary_response_spct", "summary_chroma_spct", "summary_generic_spct")
}
# is functions for spct summary classes --------------------------------------------

#' Query class of spectrum summary objects
#'
#' Functions to check if an object is of a given type of spectrum, or coerce it if
#' possible.
#'
#' @param x an R object.
#'
#' @return These functions return \code{TRUE} if its argument is a of the queried type
#'   of spectrum and \code{FALSE} otherwise.
#'
#' @note Derived types also return TRUE for a query for a base type such as
#' \code{generic_spct}.
#'
#' @export is.summary_generic_spct
#' @rdname is.summary_generic_spct
#' @examples
#' sm <- summary(sun.spct)
#' is.summary_source_spct(sm)
#'
is.summary_generic_spct <- function(x) inherits(x, "summary_generic_spct")

#' @rdname is.summary_generic_spct
#' @export
#'
is.summary_cps_spct <- function(x) inherits(x, "summary_cps_spct")

#' @rdname is.summary_generic_spct
#' @export
#'
is.summary_source_spct <- function(x) inherits(x, "summary_source_spct")

#' @rdname is.summary_generic_spct
#' @export
#'
is.summary_response_spct <- function(x) inherits(x, "summary_response_spct")

#' @rdname is.summary_generic_spct
#' @export
#'
is.summary_filter_spct <- function(x) inherits(x, "summary_filter_spct")

#' @rdname is.summary_generic_spct
#' @export
#'
is.summary_reflector_spct <- function(x) inherits(x, "summary_reflector_spct")

#' @rdname is.summary_generic_spct
#' @export
#'
is.summary_object_spct <- function(x) inherits(x, "summary_object_spct")

#' @rdname is.summary_generic_spct
#' @export
#'
is.summary_chroma_spct <- function(x) inherits(x, "summary_chroma_spct")

#' @rdname is.summary_generic_spct
#'
#' @export
#'
is.any_summary_spct <- function(x) {
  inherits(x, summary_spct_classes())
}

# summary -----------------------------------------------------------------

#' Summary of a spectral object
#'
#' Methods of generic function summary for objects of spectral classes.
#'
#' @param object An object of one of the spectral classes for which a summary is
#'   desired
#' @param maxsum integer Indicates how many levels should be shown for factors.
#' @param digits integer Used for number formatting with \code{\link{format}()}.
#' @param ... additional arguments affecting the summary produced, ignored in
#'   current version
#'
#' @return A summary object matching the class of \code{object}.
#'
#' @export
#' @method summary generic_spct
#'
#' @name summary
#'
#' @examples
#' summary(sun.spct)
#'
summary.generic_spct <- function(object,
                                 maxsum = 7,
                                 digits = max(3, getOption("digits") - 3),
                                 ...) {
  z <- list()
  class(z) <- c("summary_generic_spct", class(z))
  z[["orig.class"]] <- class_spct(object)[1]
  z[["orig.dim_desc"]] <- dplyr::dim_desc(object)
  z[["wl.range"]] <- range(object)
  z[["wl.stepsize"]] <- stepsize(object)
  z[["summary"]] <- summary(as.data.frame(object), maxsum = maxsum, digits = digits, ...)
  comment(z) <- comment(object)
  setNormalized(z, getNormalized(object))
  setScaled(z, getScaled(object))
  setWhenMeasured(z, getWhenMeasured(object))
  setWhereMeasured(z, getWhereMeasured(object))
  if (is.source_spct(object)) {
    class(z) <- c("summary_source_spct", class(z))
    setTimeUnit(z, getTimeUnit(object))
    setBSWFUsed(z, getBSWFUsed(object))
  } else if (is.response_spct(object)) {
    class(z) <- c("summary_response_spct", class(z))
    setTimeUnit(z, getTimeUnit(object))
  } else if (is.filter_spct(object)) {
    class(z) <- c("summary_filter_spct", class(z))
    setTfrType(z, getTfrType(object))
  } else if (is.reflector_spct(object)) {
    class(z) <- c("summary_reflector_spct", class(z))
    setRfrType(z, getRfrType(object))
  } else if (is.object_spct(object)) {
    class(z) <- c("summary_object_spct", class(z))
    setTfrType(z, getTfrType(object))
    setRfrType(z, getRfrType(object))
  } else if (is.chroma_spct(object)) {
    class(z) <- c("summary_chroma_spct", class(z))
  } else if (is.cps_spct(object)) {
    class(z) <- c("summary_cps_spct", class(z))
  }
  z
}

# Print spectral summaries ------------------------------------------------

#' Print spectral summary
#'
#' A function to nicely print objects of classes "summary...spct".
#'
#' @param x An object of one of the summary classes for spectra
#' @param ... not used in current version
#'
#' @export
#' @examples
#' print(summary(sun.spct))
#'
print.summary_generic_spct <- function(x, ...) {
  cat("Summary of object: ", x[["orig.class"]], " ", x[["orig.dim_desc"]], "\n", sep = "")
  cat(
    "Wavelength (nm): range ",
    paste(signif(x[["wl.range"]], 8), sep = "", collapse = " to "),
    ", step ",
    paste(unique(signif(x[["wl.stepsize"]], 7)), sep = "", collapse = " to "),
    "\n",
    sep = ""
  )
  if (!any(is.na(getWhenMeasured(x)))) {
    cat("Measured on: ",
        as.character(getWhenMeasured(x)),
        " UTC\n",
        sep = "")
  }
  if (!any(is.na(getWhereMeasured(x)))) {
    where.measured <- getWhereMeasured(x)
    cat("Measured at: ",
        where.measured[["lat"]],
        " N, ",
        where.measured[["lon"]],
        " E\n",
        sep = "")
  }
  if (class(x)[1] %in% c("summary_source_spct", "summary_response_spct")) {
    cat("Time unit: ", as.character(getTimeUnit(x, force.duration = TRUE)),
        "\n", sep = "")
  }
  if (is_scaled(x)) {
    scaling <- getScaled(x)[["f"]]
    cat("Rescaled to '", scaling, "' = 1 \n", sep = "")
  }
  if (is_normalized(x)) {
    norm <- getNormalized(x)
    cat("Spectral data normalized to 1 at ",
        norm,
        ifelse(is.numeric(norm), " nm \n", " \n"),
        sep = "")
  }
  if (is_effective(x)) {
    BSWF <- getBSWFUsed(x)
    cat("Data weighted using '", BSWF, "' BSWF\n", sep = "")
  }
  cat("\n")
  print(x[["summary"]])
  invisible(x)
}




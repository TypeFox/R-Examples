#' Insert new wavelength values into a spectrum
#'
#' Insert new wavelength values into a spectrum interpolating the corresponding
#' spectral data values.
#'
#' @param spct an object of class "generic_spct"
#' @param hinges numeric vector of wavelengths (nm) at which the
#'   s.irrad should be inserted by interpolation, no interpolation is indicated
#'   by an empty array (numeric(0))
#' @param byref logical indicating if new object will be created by reference or
#'   by copy of spct
#'
#' @return a generic_spct or a derived type with variables \code{w.length} and
#'   other numeric variables.
#'
#' @note Inserting wavelengths values "hinges" immediately before and after a
#'   discontinuity in the SWF, greatly reduces the errors caused by
#'   interpolating the weighted irradiance during integration of the effective
#'   spectral irradiance. This is specially true when data has a large
#'   wavelength step size.
#'
#' @export
#'
#' @examples
#'
#' insert_spct_hinges(sun.spct, c(399.99,400.00,699.99,700.00))
#' insert_spct_hinges(sun.spct,
#'                    c(199.99,200.00,399.50,399.99,400.00,699.99,
#'                          700.00,799.99,1000.00))
insert_spct_hinges <- function(spct, hinges=NULL, byref=FALSE) {
  if (!is.any_spct(spct)) {
    warning("Only objects derived from 'generic_spct' are supported")
    return(spct)
  }
  if (is.null(hinges) || length(hinges) == 0) {
    return(spct)
  }
  old.w.length <- spct[["w.length"]]
  if (length(hinges) > 0) {
    name <- substitute(spct)
    names.spct <- names(spct)
    names.data <- names.spct != "w.length"
    idx.wl <- which(!names.data)
    idx.data <- which(names.data)
    class_spct <- class(spct)
    comment.spct <- comment(spct)
    if (is.source_spct(spct) || is.response_spct(spct)) {
      time.unit <- getTimeUnit(spct)
      bswf.used <- getBSWFUsed(spct)
    }
    if (is.filter_spct(spct) || is.object_spct(spct)) {
      Tfr.type <- getTfrType(spct)
    }
    if (is.reflector_spct(spct) || is.object_spct(spct)) {
      Rfr.type <- getRfrType(spct)
    }
    first.iter <- TRUE
    for (data.col in idx.data) {
      temp.data <- spct[[data.col]]
      if (first.iter) {
        new.spct <- insert_hinges(old.w.length, temp.data, hinges)
        names(new.spct) <- c("w.length", names.spct[data.col])
        first.iter <- FALSE
      } else {
        new.spct[ , names.spct[data.col] ] <-
          v_insert_hinges(old.w.length, temp.data, hinges)
      }
    }
    if(class_spct[1] == "source_spct") {
      setSourceSpct(new.spct, time.unit = time.unit, bswf.used = bswf.used)
    } else if (class_spct[1] == "filter_spct") {
      setFilterSpct(new.spct, Tfr.type = Tfr.type)
    } else if (class_spct[1] == "reflector_spct") {
      setReflectorSpct(new.spct, Rfr.type = Rfr.type)
    } else if (class_spct[1] == "object_spct") {
      setObjectSpct(new.spct, Tfr.type = Tfr.type, Rfr.type = Rfr.type)
    } else if (class_spct[1] == "response_spct") {
      setResponseSpct(new.spct, time.unit = time.unit)
    } else if (class_spct[1] == "chroma_spct") {
      setChromaSpct(new.spct)
    } else if (class_spct[1] == "generic_spct") {
      setGenericSpct(new.spct)
    }
    comment(new.spct) <- comment.spct
    if (byref && is.name(name)) {
      name <- as.character(name)
      assign(name, new.spct, parent.frame(), inherits = TRUE)
    }
    return(new.spct)
  } else {
    return(spct)
  }
}

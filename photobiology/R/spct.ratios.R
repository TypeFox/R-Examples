#' Photon:photon ratio
#'
#' This function returns the photon ratio for a given pair of wavebands of a
#' light source spectrum.
#'
#' @param spct an object of class "source_spct"
#' @param w.band.num waveband definition created with new_waveband()
#' @param w.band.denom waveband definition created with new_waveband()
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded
#' @param use.cached.mult logical indicating whether multiplier values should be
#'   cached between calls
#' @param use.hinges logical indicating whether to use hinges to reduce
#'   interpolation errors
#' @param ... other arguments (possibly ignored)
#'
#' @return a single numeric nondimensional value giving a photon ratio between
#'   pairs of wavebands, with name attribute set to the name of the wavebands
#'   unless a named list of wavebands is supplied in which case the names of the
#'   list elements are used, with "(q:q)" appended.
#'
#' @export
#' @examples
#' q_ratio(sun.spct, new_waveband(400,500), new_waveband(400,700))
#'
#' @note Recycling for wavebans takes place when the number of denominator and
#'   denominator wavebands differ. The last two parameters control speed
#'   optimizations. The defaults should be suitable in mosts cases. If you will
#'   use repeatedly the same SWFs on many spectra measured at exactly the same
#'   wavelengths you may obtain some speed up by setting
#'   \code{use.cached.mult=TRUE}. However, be aware that you are responsible for
#'   ensuring that the wavelengths are the same in each call, as the only test
#'   done is for the length of the \code{w.length} vector.
#'
#' @family photon and energy ratio functions
#'
q_ratio <- function(spct, w.band.num, w.band.denom, wb.trim,
                  use.cached.mult, use.hinges, ...) UseMethod("q_ratio")

#' @describeIn q_ratio Default for generic function
#'
#' @export
#'
q_ratio.default <- function(spct, w.band.num, w.band.denom, wb.trim,
                            use.cached.mult, use.hinges, ...) {
  warning("'q_ratio' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' @describeIn q_ratio Method for \code{source_spct} objects
#'
#' @export
#'
q_ratio.source_spct <-
  function(spct,
           w.band.num = NULL, w.band.denom = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = FALSE,
           use.hinges = getOption("photobiology.use.hinges", default = NULL), ... ) {
    q.irrad.num <- irrad_spct(spct, w.band = w.band.num,
                              unit.out = "photon", quantity = "total",
                              wb.trim = wb.trim,
                              use.cached.mult = use.cached.mult,
                              use.hinges = use.hinges,
                              allow.scaled = TRUE)
    q.irrad.denom <- irrad_spct(spct, w.band = w.band.denom,
                                unit.out = "photon", quantity = "total",
                                wb.trim = wb.trim,
                                use.cached.mult = use.cached.mult,
                                use.hinges = use.hinges,
                                allow.scaled = TRUE)
    ratio <- q.irrad.num / q.irrad.denom
    names(ratio) <- paste(names(q.irrad.num), ":", names(q.irrad.denom), "(q:q)", sep = "")
    attr(ratio, "time.unit") <- NULL
    attr(ratio, "radiation.unit") <- "q:q ratio"
    return(ratio)
  }

#' Energy:energy ratio
#'
#' This function returns the photon ratio for a given pair of wavebands of a
#' light source spectrum.
#'
#' @param spct asource_spct
#' @param w.band.num waveband or list of waveband objects
#' @param w.band.denom waveband or list of waveband objects
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded
#' @param use.cached.mult logical Flag telling whether multiplier values should be
#'   cached between calls
#' @param use.hinges logical Flag telling whether to use hinges to reduce
#'   interpolation errors
#' @param ... other arguments (possibly ignored)
#'
#' @return A single numeric nondimensional value giving an energy ratio between
#'   pairs of wavebands, with name attribute set to the name of each waveband
#'   unless a named list of wavebands is supplied in which case the names of the
#'   list elements are used, with "(e:e)" appended.
#'
#'
#' @export e_ratio
#' @examples
#' e_ratio(sun.spct, new_waveband(400,500), new_waveband(400,700))
#'
#' @note Recycling for wavebans takes place when the number of denominator and
#'   denominator wavebands differ. The last two parameters control speed
#'   optimizations. The defaults should be suitable in mosts cases. If you will
#'   use repeatedly the same SWFs on many spectra measured at exactly the same
#'   wavelengths you may obtain some speed up by setting
#'   \code{use.cached.mult=TRUE}. However, be aware that you are responsible for
#'   ensuring that the wavelengths are the same in each call, as the only test
#'   done is for the length of the \code{w.length} vector.
#'
#' @family photon and energy ratio functions
#'
e_ratio <- function(spct, w.band.num, w.band.denom, wb.trim,
                    use.cached.mult, use.hinges, ...) UseMethod("q_ratio")

#' @describeIn e_ratio Default for generic function
#'
#' @export
#'
e_ratio.default <- function(spct, w.band.num, w.band.denom, wb.trim,
                            use.cached.mult, use.hinges, ...) {
  warning("'e_ratio' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' @describeIn e_ratio Method for \code{source_spct} objects
#'
#' @export
#'
e_ratio.source_spct <-
  function(spct,
           w.band.num = NULL, w.band.denom = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = FALSE,
           use.hinges=getOption("photobiology.use.hinges", default = NULL), ...) {
    e.irrad.num <- irrad_spct(spct, w.band = w.band.num, unit.out = "energy", quantity = "total",
                              wb.trim = wb.trim,
                              use.cached.mult = use.cached.mult, use.hinges = use.hinges,
                              allow.scaled=TRUE)
    e.irrad.denom <- irrad_spct(spct, w.band = w.band.denom, unit.out = "energy", quantity = "total",
                                wb.trim = wb.trim,
                                use.cached.mult = use.cached.mult, use.hinges = use.hinges,
                                allow.scaled = TRUE)
    ratio <- e.irrad.num / e.irrad.denom
    names(ratio) <- paste(names(e.irrad.num), ":", names(e.irrad.denom), "(e:e)", sep="")
    attr(ratio, "time.unit") <- NULL
    attr(ratio, "radiation.unit") <- "e:e ratio"
    return(ratio)
  }

#' Photon:energy ratio
#'
#' This function returns the photon to energy ratio for each waveband of a light
#' source spectrum.
#'
#' @param spct source_spct
#' @param w.band waveband or list of waveband objects
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded
#' @param use.cached.mult logical Flag telling whether multiplier values should be
#'   cached between calls
#' @param use.hinges logical Flag telling whether to use hinges to reduce
#'   interpolation errors
#' @param ... other arguments (possibly ignored)
#'
#' @return A vector of \code{numeric} values giving number of moles of photons
#'   per Joule for each waveband, with name attribute set to the name of each
#'   waveband unless a named list of wavebands is supplied in which case the
#'   names of the list elements are used, with "q:e" prepended..
#'
#'
#' @export
#' @examples
#' qe_ratio(sun.spct, new_waveband(400,700))
#'
#' @note The last two parameters control speed optimizations. The defaults
#'   should be suitable in mosts cases. If you will use repeatedly the same SWFs
#'   on many spectra measured at exactly the same wavelengths you may obtain
#'   some speed up by setting \code{use.cached.mult=TRUE}. However, be aware
#'   that you are responsible for ensuring that the wavelengths are the same in
#'   each call, as the only test done is for the length of the \code{w.length}
#'   vector.
#'
#' @family photon and energy ratio functions
#'
qe_ratio <- function(spct, w.band, wb.trim,
                     use.cached.mult, use.hinges, ...) UseMethod("qe_ratio")

#' @describeIn qe_ratio Default for generic function
#'
#' @export
#'
qe_ratio.default <- function(spct, w.band, wb.trim,
                             use.cached.mult, use.hinges, ...) {
  warning("'qe_ratio' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' @describeIn qe_ratio Method for \code{source_spct} objects
#'
#' @export
#'
qe_ratio.source_spct <-
  function(spct, w.band = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = FALSE,
           use.hinges = getOption("photobiology.use.hinges", default = NULL), ...) {
    q.irrad <- irrad_spct(spct, w.band=w.band, unit.out = "photon",
                          quantity ="total",
                          wb.trim = wb.trim,
                          use.cached.mult = use.cached.mult,
                          use.hinges = use.hinges,
                          allow.scaled = TRUE)
    e.irrad <- irrad_spct(spct, w.band=w.band, unit.out = "energy",
                          quantity = "total",
                          wb.trim = wb.trim,
                          use.cached.mult = use.cached.mult,
                          use.hinges = use.hinges,
                          allow.scaled = TRUE)
    ratio <- q.irrad / e.irrad
    names(ratio) <- paste("q:e(", names(q.irrad), ")", sep = "")
    attr(ratio, "time.unit") <- NULL
    attr(ratio, "radiation.unit") <- "q:e ratio"
    return(ratio)
  }

#' Energy:photon ratio
#'
#' This function returns the energy to mole of photons ratio for each waveband and a
#' light source spectrum.
#'
#' @param spct source_spct
#' @param w.band waveband or list of waveband objects
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded
#' @param use.cached.mult logical Flag telling whether multiplier values should
#'   be cached between calls
#' @param use.hinges logical Flag telling whether to use hinges to reduce
#'   interpolation errors
#' @param ... other arguments (possibly ignored)
#'
#' @return a numeric value giving number of Joule per mol of photons for each
#'   waveband, with name attribute set to the name of each waveband unless a
#'   named list of wavebands is supplied in which case the names of the list
#'   elements are used, with "e:q" prepended..
#'
#'
#' @export
#' @examples
#' eq_ratio(sun.spct, new_waveband(400,700))
#'
#' @note The last two parameters control speed optimizations. The defaults
#'   should be suitable in mosts cases. If you will use repeatedly the same SWFs
#'   on many spectra measured at exactly the same wavelengths you may obtain
#'   some speed up by setting \code{use.cached.mult=TRUE}. However, be aware
#'   that you are responsible for ensuring that the wavelengths are the same in
#'   each call, as the only test done is for the length of the \code{w.length}
#'   vector.
#'
#' @family photon and energy ratio functions
#'
eq_ratio <- function(spct, w.band, wb.trim,
                     use.cached.mult, use.hinges, ...) UseMethod("eq_ratio")

#' @describeIn eq_ratio Default for generic function
#'
#' @export
#'
eq_ratio.default <- function(spct, w.band, wb.trim,
                             use.cached.mult, use.hinges, ...) {
  warning("'eq_ratio' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' @describeIn eq_ratio Method for \code{source_spct} objects
#'
#' @export
#'
eq_ratio.source_spct <-
  function(spct, w.band=NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult =FALSE,
           use.hinges  = getOption("photobiology.use.hinges", default = NULL), ...) {
    ratio <- 1 / qe_ratio(spct = spct, w.band = w.band, wb.trim = wb.trim,
                          use.cached.mult = use.cached.mult, use.hinges = use.hinges)
    names(ratio) <- gsub("q:e", "e:q", names(ratio), fixed = TRUE )
    attr(ratio, "time.unit") <- NULL
    attr(ratio, "radiation.unit") <- "e:q ratio"
    return(ratio)
  }

# source_mspct methods ----------------------------------------------------

#' @describeIn q_ratio Calculates photon:photon from a \code{source_mspct}
#'   object.
#'
#' @param idx logical whether to add a column with the names of the elements of spct
#'
#' @export
#'
q_ratio.source_mspct <-
  function(spct,
           w.band.num = NULL, w.band.denom = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = FALSE,
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           ..., idx = !is.null(names(spct))) {
    msdply(
      mspct = spct,
      .fun = q_ratio,
      w.band.num = w.band.num,
      w.band.denom = w.band.denom,
      wb.trim = wb.trim,
      use.cached.mult = use.cached.mult,
      use.hinges = use.hinges,
      idx = idx
    )
  }

#' @describeIn e_ratio Calculates energy:energy ratio from a \code{source_mspct}
#'   object.
#'
#' @param idx logical whether to add a column with the names of the elements of spct
#'
#' @export
#'
e_ratio.source_mspct <-
  function(spct,
           w.band.num = NULL, w.band.denom = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = FALSE,
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           ..., idx = !is.null(names(spct))) {
    msdply(
      mspct = spct,
      .fun = e_ratio,
      w.band.num = w.band.num,
      w.band.denom = w.band.denom,
      wb.trim = wb.trim,
      use.cached.mult = use.cached.mult,
      use.hinges = use.hinges,
      idx = idx
    )
  }

#' @describeIn eq_ratio Calculates energy:photon from a \code{source_mspct}
#'   object.
#'
#' @param idx logical whether to add a column with the names of the elements of spct
#'
#' @export
#'
eq_ratio.source_mspct <-
  function(spct, w.band = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = FALSE,
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           ...,
           idx = !is.null(names(spct))) {
    msdply(
      mspct = spct,
      .fun = eq_ratio,
      w.band = w.band,
      wb.trim = wb.trim,
      use.cached.mult = use.cached.mult,
      use.hinges = use.hinges,
      idx = idx,
      col.names = names(w.band)
    )
  }

#' @describeIn qe_ratio Calculates photon:energy ratio from a
#'   \code{source_mspct} object.
#'
#' @param idx logical whether to add a column with the names of the elements of spct
#'
#' @export
#'
qe_ratio.source_mspct <-
  function(spct, w.band=NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.cached.mult = FALSE,
           use.hinges=getOption("photobiology.use.hinges", default = NULL),
           ...,
           idx = !is.null(names(spct))) {
    msdply(
      spct,
      .fun = qe_ratio,
      w.band = w.band,
      wb.trim = wb.trim,
      use.cached.mult = use.cached.mult,
      use.hinges = use.hinges,
      idx = idx,
      col.names = names(w.band)
    )
  }



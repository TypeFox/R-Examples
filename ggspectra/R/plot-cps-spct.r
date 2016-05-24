#' Plot an instrument counts spectrum.
#'
#' This function returns a ggplot object with an annotated plot of a
#' cps_spct object.
#'
#' @note Note that scales are expanded so as to make space for the annotations.
#' The object returned is a ggplot objects, and can be further manipulated.
#'
#' @param spct a cps_spct object
#' @param w.band list of waveband objects
#' @param range an R object on which range() returns a vector of length 2,
#' with min annd max wavelengths (nm)
#' @param pc.out logical, if TRUE use percents instead of fraction of one
#' @param label.qty character string giving the type of summary quantity to use
#'   for labels
#' @param annotations a character vector
#' @param norm numeric normalization wavelength (nm) or character string "max"
#' for normalization at the wavelength of highest peak.
#' @param ... other arguments
#'
#' @return a \code{ggplot} object.
#'
#' @keywords internal
#'
cps_plot <- function(spct,
                     w.band,
                     range,
                     pc.out,
                     label.qty,
                     annotations,
                     norm,
                     ...) {
  if (!is.cps_spct(spct)) {
    stop("cps_plot() can only plot response_spct objects.")
  }
  if (!is.null(range)) {
    trim_spct(spct, range = range, byref = TRUE)
  }
  if (is.null(norm)) {
    # we will use the original data
    scale.factor <- 1
  } else if (!is.null(norm)) {
    if (is.character(norm)) {
      if (norm %in% c("max", "maximum")) {
        idx <- which.max(spct[["cps"]])
      } else {
        warning("Invalid character '", norm, "'value in 'norm'")
        return(ggplot())
      }
      scale.factor <- 1 / spct[idx, "cps"]
      norm <- spct[idx, "w.length"]
    } else if (is.numeric(norm) && norm >= min(spct) && norm <= max(spct)) {
      scale.factor <- 1 / interpolate_spct(spct, norm)[["cps"]]
    } else if (is.numeric(norm)) {
      warning("'norm = ", norm, "' value outside spectral data range of ",
              round(min(spct)), " to ", round(max(spct)), " (nm)")
      return(ggplot())
    } else {
      stop("'norm' should be numeric or character")
    }
  }

  if (scale.factor != 1) {
    if (!pc.out) {
      multiplier.label <- "rel."
      #      scale.factor <- 1 * scale.factor
    } else {
      multiplier.label <- "%"
      scale.factor <- 100 * scale.factor
    }
    if (is.numeric(norm)) {
      norm <- signif(norm, digits = 4)
    }
    s.cps.label <-
      bquote(Pixel~~response~~rate~~N( italic(lambda) )/N( .(norm))~~(.(multiplier.label)))
    cps.label <- ""
  } else {
    s.cps.label <-
      expression(Pixel~~response~~rate~~N(lambda)~~(counts~~s^{-1}))
    cps.label <- ""
  }

  spct[["cps"]] <- spct[["cps"]] * scale.factor
  y.max <- max(spct[["cps"]], na.rm = TRUE)
  y.min <- 0
  plot <- ggplot(spct)  +
    scale_fill_identity() + scale_color_identity()
  plot <- plot + geom_line()
  plot <- plot + labs(x = "Wavelength (nm)", y = s.cps.label)

  plot <- plot + decoration(w.band = w.band,
                            y.max = y.max,
                            y.min = y.min,
                            x.max = max(spct),
                            x.min = min(spct),
                            annotations = annotations,
                            label.qty = label.qty,
                            summary.label = cps.label)

  if (!is.null(annotations) &&
      length(intersect(c("boxes", "segments", "labels", "summaries", "colour.guide"), annotations)) > 0L) {
    y.limits <- c(0, y.max * 1.25)
    x.limits <- c(min(spct) - spread(spct) * 0.025, NA)
  } else {
    y.limits <- c(0, 1)
    x.limits <- range(spct)
  }
  plot <- plot + scale_y_continuous(limits = y.limits)
  plot + scale_x_continuous(limits = x.limits)

}


#' Plot an instrument counts spectrum, especialization of generic plot function.
#'
#' This function returns a ggplot object with an annotated plot of a
#' response_spct object.
#'
#' @note Note that scales are expanded so as to make space for the annotations.
#' The object returned is a ggplot objects, and can be further manipulated.
#'
#' @param x a cps_spct object
#' @param ... other arguments passed along, such as \code{label.qty}
#' @param w.band a single waveband object or a list of waveband objects
#' @param range an R object on which range() returns a vector of length 2,
#' with min annd max wavelengths (nm)
#' @param unit.out character IGNORED
#' @param pc.out logical, if TRUE use percents instead of fraction of one
#' @param label.qty character IGNORED
#' @param annotations a character vector ("summaries" is ignored)
#' @param norm numeric normalization wavelength (nm) or character string "max"
#' for normalization at the wavelength of highest peak.
#'
#' @return a \code{ggplot} object.
#'
#' @export
#'
#' @keywords hplot
#'
#' @family plot functions
#'
plot.cps_spct <-
  function(x, ...,
           w.band = getOption("photobiology.plot.bands", default = list(UVC(), UVB(), UVA(), PAR())),
           range = NULL,
           unit.out = "cps",
           pc.out = FALSE,
           label.qty = "average",
           annotations = getOption("photobiology.plot.annotations",
                                 default = c("boxes", "labels", "colour.guide", "peaks")),
           norm = NULL ) {
    if ("color.guide" %in% annotations) {
      annotations <- c(setdiff(annotations, "color.guide"), "colour.guide")
    }

    out.ggplot <- cps_plot(spct = x, w.band = w.band, range = range,
                           label.qty = label.qty,
                           pc.out = pc.out,
                           annotations = annotations, norm = norm, ...)
    if ("title" %in% annotations) {
      out.ggplot <- out.ggplot + labs(title = deparse(substitute(x)))
    }
    out.ggplot
  }



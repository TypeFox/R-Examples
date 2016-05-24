#' Plot a source spectrum.
#'
#' This function returns a ggplot object with an annotated plot of a source_spct
#' object.
#'
#' @note Note that scales are expanded so as to make space for the annotations.
#'   The object returned is a ggplot objects, and can be further manipulated.
#'
#' @param spct a source_spct object
#' @param w.band list of waveband objects
#' @param range an R object on which range() returns a vector of length 2, with
#'   min annd max wavelengths (nm)
#' @param label.qty character string giving the type of summary quantity to use
#'   for labels
#' @param annotations a character vector
#' @param ... other arguments passed to annotate_plot()
#'
#' @return a \code{ggplot} object.
#'
#' @keywords internal
#'
e_plot <- function(spct,
                   w.band,
                   range,
                   label.qty,
                   annotations,
                   ...) {
  if (!is.source_spct(spct)) {
    stop("e_plot() can only plot source_spct objects.")
  }
  q2e(spct, byref=TRUE)
  if (!is.null(range)) {
    trim_spct(spct, range = range, byref = TRUE)
  }

  exposure.label <- NA
  if (is_normalized(spct) || is_scaled(spct)) {
    s.irrad.label <- "Spectral~~energy~~exposure~~E(lambda)~~(relative~~units)"
    irrad.label.total <- "atop(E, (relative~~units))"
    irrad.label.avg <- "atop(bar(E(lambda)), (relative~~units))"
    scale.factor <- 1
  } else {
    time.unit <- getTimeUnit(spct)
    if (!length(time.unit)) {
      time.unit <- "unkonwn"
    }
    if (time.unit == "second" || time.unit == lubridate::duration(1, "seconds"))  {
      s.irrad.label <- "Spectral~~energy~~irradiance~~E(lambda)~~(W~m^{-2}~nm^{-1})"
      irrad.label.total  <- "atop(E, (W~m^{-2}))"
      irrad.label.avg  <- "atop(bar(E(lambda)), (W~m^{-2}~nm^{-1}))"
      scale.factor <- 1
    } else if (time.unit == "day" || time.unit == lubridate::duration(1, "days")) {
      s.irrad.label <- "Spectral~~energy~~exposure~~E(lambda)~~(MJ~d^{-1}~m^{-2}~nm^{-1})"
      irrad.label.total <- "atop(E, (MJ~d^{-1}~m^{-2}))"
      irrad.label.avg <- "atop(bar(E(lambda)), (MJ~d^{-1}~m^{-2}~nm^{-1}))"
      scale.factor <- 1e-6
    } else if (time.unit == "hour" || time.unit == lubridate::duration(1, "hours")) {
      s.irrad.label <- "Spectral~~energy~~exposure~~E(lambda)~~(kJ~h^{-1}~m^{-2}~nm^{-1})"
      irrad.label.total <- "atop(E, (kJ~h^{-1}~m^{-2}))"
      irrad.label.avg <- "atop(bar(E(lambda)), (kJ~h^{-1}~m^{-2}~nm^{-1}))"
      scale.factor <- 1e-3
    } else if (time.unit == "exposure" || lubridate::is.duration(time.unit)) {
      s.irrad.label <- "Spectral~~energy~~fluence~~E(lambda)~~(kJ~m^{-2}~nm^{-1})"
      irrad.label.total <- "atop(E, (kJ~m^{-2}))"
      irrad.label.avg <- "atop(bar(E(lambda)), (kJ~m^{-2}~nm^{-1}))"
      exposure.label <- paste("Length of exposure:",
                              ifelse(lubridate::is.duration(time.unit),
                              as.character(time.unit), "unknown"))
      scale.factor <- 1e-3
    } else {
      s.irrad.label <- "Spectral~~energy~~exposure~~E(lambda)~~(arbitrary~~units)"
      irrad.label.total <- "atop(E, (arbitrary~~units))"
      irrad.label.avg <- "atop(bar(E(lambda)), (arbitrary~~units))"
      scale.factor <- 1
    }
  }
  if (label.qty == "total") {
    irrad.label <- irrad.label.total
  } else if (label.qty %in% c("average", "mean")) {
    irrad.label <- irrad.label.avg
  }  else if (label.qty == "contribution") {
    irrad.label <- "atop(Contribution~~to~~total, E~~(fraction))"
  } else if (label.qty == "contribution.pc") {
    irrad.label <- "atop(Contribution~~to~~total, E~~(percent))"
  } else if (label.qty == "relative") {
    irrad.label <- "atop(Relative~~to~~sum, E~~(fraction))"
  } else if (label.qty == "relative.pc") {
    irrad.label <- "atop(Relative~~to~~sum, E~~(percent))"
  } else {
    irrad.label <- ""
  }
  if (is_effective(spct)) {
    s.irrad.label <- sub("E", "E[eff]", s.irrad.label, fixed = TRUE)
    irrad.label <- sub("E", "E[eff]", irrad.label, fixed = TRUE)
    irrad.label.total <- sub("E", "E[eff]", irrad.label.total, fixed = TRUE)
    irrad.label.avg <- sub("E", "E[eff]", irrad.label.avg, fixed = TRUE)
  }
  s.irrad.label <- parse(text = s.irrad.label)
  spct[["s.e.irrad"]] <- spct[["s.e.irrad"]] * scale.factor
  y.max <- max(spct[["s.e.irrad"]], na.rm = TRUE)
  y.min <- 0
  plot <- ggplot(spct, aes_(~w.length, ~s.e.irrad)) +
    scale_fill_identity() + scale_color_identity()
  plot <- plot + geom_line()
  plot <- plot + labs(x = "Wavelength (nm)", y = s.irrad.label)

  if (label.qty == "total") {
    label.qty <- "irrad"
  } else if (label.qty %in% c("mean", "average")) {
    label.qty <- "sirrad"
  }
  plot <- plot + decoration(w.band = w.band,
                            unit.out = "energy",
                            time.unit = getTimeUnit(spct),
                            y.max = y.max,
                            y.min = y.min,
                            x.max = max(spct),
                            x.min = min(spct),
                            annotations = annotations,
                            label.qty = label.qty,
                            summary.label = irrad.label)

  if (is_effective(spct)) {
    plot <- plot +  annotate("text",
                             x = midpoint(spct),
                             y = y.max,
                             label = paste("BSWF:", getBSWFUsed(spct)),
                             vjust = -0.5, size = rel(3) )
  }

  if (!is.na(exposure.label)) {
    plot <- plot +  annotate("text",
                             x = min(spct),
                             y = y.max,
                             label = exposure.label,
                             vjust = -0.5,
                             hjust = 0,
                             size = rel(3) )
  }

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

#' Plot a source spectrum.
#'
#' This function returns a ggplot object with an annotated plot of a source_spct
#' object.
#'
#' @note Note that scales are expanded so as to make space for the annotations.
#'   The object returned is a ggplot objects, and can be further manipulated.
#'
#' @param spct a source_spct object
#' @param w.band list of waveband objects
#' @param range an R object on which range() returns a vector of length 2, with
#'   min annd max wavelengths (nm)
#' @param label.qty character string giving the type of summary quantity to use
#'   for labels
#' @param annotations a character vector
#' @param ... other arguments passed to annotate_plot()
#'
#' @return a \code{ggplot} object.
#'
#' @keywords internal
#'
q_plot <- function(spct,
                   w.band,
                   range,
                   label.qty,
                   annotations,
                   ...) {
  if (!is.source_spct(spct)) {
    stop("q_plot() can only plot source_spct objects.")
  }
  e2q(spct, byref = TRUE)
  if (!is.null(range)) {
    trim_spct(spct, range = range, byref = TRUE)
  }

  exposure.label <- NA
  if (is_normalized(spct) || is_scaled(spct)) {
    s.irrad.label <- "Spectral~~photon~~exposure~~Q(lambda)~~(relative~~units)"
    irrad.label.total <- "atop(Q, (relative~~units))"
    irrad.label.avg <- "atop(bar(Q(lambda)), (relative~~units))"
    scale.factor <- 1
  } else {
    time.unit <- getTimeUnit(spct)
    if (!length(time.unit)) {
      time.unit <- "unkonwn"
    }
    if (time.unit=="second" || time.unit == lubridate::duration(1, "seconds")) {
      s.irrad.label <- "Spectral~~photon~~irradiance~~Q(lambda)~~(mu*mol~s^{-1}~m^{-2}~nm^{-1})"
      irrad.label.total  <- "atop(Q, (mu*mol~s^{-1}~m^{-2}))"
      irrad.label.avg  <- "atop(bar(Q(lambda)), (mu*mol~s^{-1}~m^{-2}~nm^{-1}))"
      scale.factor <- 1e6
    } else if (time.unit=="day" || time.unit == lubridate::duration(1, "days")) {
      s.irrad.label <- "Spectral~~photon~~exposure~~Q(lambda)~~(mol~d^{-1}~m^{-2}~nm^{-1})"
      irrad.label.total <- "atop(Q, (mol~d^{-1}~m^{-2}))"
      irrad.label.avg <- "atop(bar(Q(lambda)), (mol~d^{-1}~m^{-2}~nm^{-1}))"
      scale.factor <- 1
    } else if (time.unit=="hour" || time.unit == lubridate::duration(1, "hours")) {
      s.irrad.label <- "Spectral~~photon~~exposure~~Q(lambda)~~(mmol~h^{-1}~m^{-2}~nm^{-1})"
      irrad.label.total <- "atop(Q, (mmol~h^{-1}~m^{-2}))"
      irrad.label.avg <- "atop(bar(Q(lambda)), (mmol~h^{-1}~m^{-2}~nm^{-1}))"
      scale.factor <- 1e3
    } else if (time.unit=="exposure" || lubridate::is.duration(time.unit)) {
      s.irrad.label <- "Spectral~~photon~~fluence~~Q(lambda)~~(mol~m^{-2}~nm^{-1})"
      irrad.label.total <- "atop(Q, (mol~m^{-2}))"
      irrad.label.avg <- "atop(bar(Q(lambda)), (mol~m^{-2}~nm^{-1}))"
      exposure.label <- paste("Length of exposure:",
                              ifelse(lubridate::is.duration(time.unit),
                                     as.character(time.unit), "unknown"))
      scale.factor <- 1
    } else {
      s.irrad.label <- "Spectral~~photon~~exposure~~Q(lambda)~~(arbitrary~~units)"
      irrad.label.total <- "atop(Q, (arbitrary~~units))"
      irrad.label.avg <- "atop(bar(Q(lambda)), (arbitrary~~units))"
      scale.factor <- 1
    }
  }
  if (label.qty == "total") {
    irrad.label <- irrad.label.total
  } else if (label.qty %in% c("average", "mean")) {
    irrad.label <- irrad.label.avg
  } else if (label.qty == "contribution") {
    irrad.label <- "atop(Contribution~~to~~total, Q~~(fraction))"
  } else if (label.qty == "contribution.pc") {
    irrad.label <- "atop(Contribution~~to~~total, Q~~(percent))"
  } else if (label.qty == "relative") {
    irrad.label <- "atop(Relative~~to~~sum, Q~~(fraction))"
  } else if (label.qty == "relative.pc") {
    irrad.label <- "atop(Relative~~to~~sum, Q~~(percent))"
  } else {
    irrad.label <- ""
  }
  if (is_effective(spct)) {
    s.irrad.label <- sub("Q", "Q[eff]", s.irrad.label, fixed = TRUE)
    irrad.label <- sub("Q", "Q[eff]", irrad.label, fixed = TRUE)
    irrad.label.total <- sub("Q", "Q[eff]", irrad.label.total, fixed = TRUE)
    irrad.label.avg <- sub("Q", "Q[eff]", irrad.label.avg, fixed = TRUE)
  }
  s.irrad.label <- parse(text = s.irrad.label)
  spct[["s.q.irrad"]] <- spct[["s.q.irrad"]] * scale.factor
  y.max <- max(spct[["s.q.irrad"]], na.rm = TRUE)
  y.min <- 0
  plot <- ggplot(spct, aes_(~w.length, ~s.q.irrad)) +
    scale_fill_identity() + scale_color_identity()
  plot <- plot + geom_line()
  plot <- plot + labs(x = "Wavelength (nm)", y = s.irrad.label)

  if (label.qty == "total") {
    label.qty <- "irrad"
  } else if (label.qty %in% c("mean", "average")) {
    label.qty <- "sirrad"
  }
  plot <- plot + decoration(w.band = w.band,
                            unit.out = "photon",
                            time.unit = getTimeUnit(spct),
                            y.max = y.max,
                            y.min = y.min,
                            x.max = max(spct),
                            x.min = min(spct),
                            annotations = annotations,
                            label.qty = label.qty,
                            summary.label = irrad.label)

  if (is_effective(spct)) {
    plot <- plot +  annotate("text",
                             x = midpoint(spct),
                             y = y.max,
                             label = paste("BSWF:", getBSWFUsed(spct)),
                             vjust = -0.5, size = rel(3) )
  }

  if (!is.na(exposure.label)) {
    plot <- plot +  annotate("text",
                             x = min(spct),
                             y = y.max,
                             label = exposure.label,
                             vjust = -0.5,
                             hjust = 0,
                             size = rel(3) )
  }

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

#' Plot a source spectrum, especialization of generic plot function.
#'
#' This function returns a ggplot object with an annotated plot of a source_spct
#' object.
#'
#' @note Note that scales are expanded so as to make space for the annotations.
#'   The object returned is a ggplot objects, and can be further manipulated.
#'
#' @param x a source_spct object
#' @param ... other arguments passed along, such as \code{label.qty}
#' @param w.band a single waveband object or a list of waveband objects
#' @param range an R object on which range() returns a vector of length 2, with
#'   min annd max wavelengths (nm)
#' @param unit.out character string indicating type of radiation units to use
#'   for plotting: "photon" or its synomin "quantum", or "energy"
#' @param label.qty character string giving the type of summary quantity to use
#'   for labels
#' @param annotations a character vector
#'
#' @return a \code{ggplot} object.
#'
#' @method plot source_spct
#' @export
#'
#' @keywords hplot
#'
#' @examples
#' library(photobiology)
#' plot(sun.spct)
#' plot(sun.spct, unit.out = "photon")
#' plot(sun.daily.spct)
#' plot(sun.daily.spct, unit.out = "photon")
#'
#' @family plot functions
#'
plot.source_spct <-
  function(x, ...,
           w.band=getOption("photobiology.plot.bands", default = list(UVC(), UVB(), UVA(), PAR())),
           range=NULL,
           unit.out=getOption("photobiology.radiation.unit", default = "energy"),
           label.qty = "total",
           annotations=getOption("photobiology.plot.annotations",
                                 default = c("boxes", "labels", "summaries", "colour.guide", "peaks")) ) {
    if ("color.guide" %in% annotations) {
      annotations <- c(setdiff(annotations, "color.guide"), "colour.guide")
    }
    if (unit.out %in% c("photon", "quantum")) {
      out.ggplot <- q_plot(spct = x, w.band = w.band, range = range,
                           label.qty = label.qty, annotations = annotations, ...)
    } else if (unit.out == "energy") {
      out.ggplot <- e_plot(spct = x, w.band = w.band, range = range,
                           label.qty = label.qty, annotations = annotations, ...)
    } else {
      stop("Invalid 'radiation.unit' argument value: '", unit.out, "'")
    }
    if ("title" %in% annotations) {
      out.ggplot <- out.ggplot + labs(title = deparse(substitute(x)))
    }
    out.ggplot
  }

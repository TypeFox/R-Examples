#' Plot a filter spectrum.
#'
#' This function returns a ggplot object with an annotated plot of a source_spct
#' object showing transmittance.
#'
#' @note Note that scales are expanded so as to make space for the annotations.
#'   The object returned is a ggplot objects, and can be further manipulated.
#'
#' @param spct a filter_spct object
#' @param w.band list of waveband objects
#' @param range an R object on which range() returns a vector of length 2, with
#'   min annd max wavelengths (nm)
#' @param pc.out logical, if TRUE use percents instead of fraction of one
#' @param label.qty character string giving the type of summary quantity to use
#'   for labels
#' @param annotations a character vector
#' @param ... other arguments passed to transmittance()
#'
#' @return a \code{ggplot} object.
#'
#' @keywords internal
#'
T_plot <- function(spct,
                   w.band,
                   range,
                   pc.out,
                   label.qty,
                   annotations,
                   ...) {
  if (!is.filter_spct(spct)) {
    stop("T_plot() can only plot filter_spct objects.")
  }
  A2T(spct, byref=TRUE)
  Tfr.type <- getTfrType(spct)
  if (!is.null(range)) {
    trim_spct(spct, range = range, byref = TRUE)
  }
  if (!length(Tfr.type)) {
    Tfr.type <- "unknown"
  }
  if (!pc.out) {
    scale.factor <- 1
    if (Tfr.type == "internal") {
      s.Tfr.label <- expression(Internal~~spectral~~transmittance~~T[int](lambda)~~(fraction))
      Tfr.label.total  <- "atop(T[int], (fraction))"
      Tfr.label.avg  <- "atop(bar(T[int](lambda)), (fraction))"
    } else if (Tfr.type == "total") {
      s.Tfr.label <- expression(Total~~spectral~~transmittance~~T[tot](lambda)~~(fraction))
      Tfr.label.total  <- "atop(T[tot], (total))"
      Tfr.label.avg  <- "atop(bar(T[tot](lambda)), (fraction))"
    }  else {
      s.Tfr.label <- expression(Spectral~~transmittance~~T(lambda)~~(fraction))
      Tfr.label.total  <- "atop(T, (total))"
      Tfr.label.avg  <- "atop(bar(T(lambda)), (fraction))"
    }
  } else if (pc.out) {
    scale.factor <- 100
    if (Tfr.type == "internal") {
      s.Tfr.label <- expression(Internal~~spectral~~transmittance~~T[int](lambda)~~(percent))
      Tfr.label.total  <- "atop(T[int], (total %*% 100))"
      Tfr.label.avg  <- "atop(bar(T[int](lambda)), (percent))"
    } else if (Tfr.type == "total") {
      s.Tfr.label <- expression(Total~~spectral~~transmittance~~T[tot](lambda)~~(percent))
      Tfr.label.total  <- "atop(T[tot], (total %*% 100))"
      Tfr.label.avg  <- "atop(bar(T[tot](lambda)), (percent))"
    }  else {
      s.Tfr.label <- expression(Spectral~~transmittance~~T(lambda)~~(percent))
      Tfr.label.total  <- "atop(T, (total  %*% 100))"
      Tfr.label.avg  <- "atop(bar(T(lambda)), (percent))"
    }
  }
  if (label.qty == "total") {
    Tfr.label <- Tfr.label.total
  } else if (label.qty %in% c("average", "mean")) {
    Tfr.label <- Tfr.label.avg
  } else if (label.qty == "contribution") {
    Tfr.label <- "atop(Contribution~~to~~total, T~~(fraction))"
  } else if (label.qty == "contribution.pc") {
    Tfr.label <- "atop(Contribution~~to~~total, T~~(percent))"
  } else if (label.qty == "relative") {
    Tfr.label <- "atop(Relative~~to~~sum, T~~(fraction))"
  } else if (label.qty == "relative.pc") {
    Tfr.label <- "atop(Relative~~to~~sum, T~~(percent))"
  } else {
    Tfr.label <- ""
  }
  y.max <- 1
  y.min <- 0

  plot <- ggplot(spct, aes_(~w.length, ~Tfr)) +
    scale_fill_identity() + scale_color_identity()
  plot <- plot + geom_line()
  plot <- plot + labs(x="Wavelength (nm)", y=s.Tfr.label)

  plot <- plot + decoration(w.band = w.band,
                            label.mult = scale.factor,
                            y.max = y.max,
                            y.min = y.min,
                            x.max = max(spct),
                            x.min = min(spct),
                            annotations = annotations,
                            label.qty = label.qty,
                            summary.label = Tfr.label)

  if (!is.null(annotations) &&
      length(intersect(c("labels", "summaries", "colour.guide"), annotations)) > 0L) {
    y.limits <- c(0, y.max * 1.25)
    x.limits <- c(min(spct) - spread(spct) * 0.025, NA)
  } else {
    y.limits <- c(0, 1)
    x.limits <- range(spct)
  }
  if (pc.out) {
    plot <- plot + scale_y_continuous(labels = scales::percent, breaks = c(0, 0.25, 0.5, 0.75, 1),
                                      limits = y.limits)
  } else {
    plot <- plot + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                                      limits = y.limits)
  }

  plot + scale_x_continuous(limits = x.limits)

}

#' Plot a filter spectrum.
#'
#' This function returns a ggplot object with an annotated plot of a source_spct
#' object showing spectral absorbance.
#'
#' @note Note that scales are expanded so as to make space for the annotations.
#'   The object returned is a ggplot objects, and can be further manipulated.
#'
#' @param spct a filter_spct object
#' @param w.band list of waveband objects
#' @param range an R object on which range() returns a vector of length 2, with
#'   min annd max wavelengths (nm)
#' @param label.qty character string giving the type of summary quantity to use
#'   for labels
#' @param annotations a character vector
#' @param ... other arguments passed to absorbance()
#'
#' @return a \code{ggplot} object.
#'
#' @keywords internal
#'
A_plot <- function(spct,
                   w.band,
                   range,
                   label.qty,
                   annotations,
                   ...) {
  if (!is.filter_spct(spct)) {
    stop("T_plot() can only plot filter_spct objects.")
  }
  T2A(spct, action="replace", byref=TRUE)
  if (!is.null(range)) {
    trim_spct(spct, range = range, byref = TRUE)
  }
  Tfr.type <- getTfrType(spct)
  if (!length(Tfr.type)) {
    Tfr.type <- "unknown"
  }
  if (Tfr.type=="internal") {
    s.A.label <- expression(Internal~~spectral~~absorbance~~A[int](lambda)~~(AU))
    A.label.total  <- "atop(A[int], (AU %*% nm))"
    A.label.avg  <- "atop(bar(A[int](lambda)), (AU))"
  } else if (Tfr.type=="total") {
    s.A.label <- expression(Total~~spectral~~absorbance~~A[tot](lambda)~~(AU))
    A.label.total  <- "atop(A[tot], (AU %*% nm))"
    A.label.avg  <- "atop(bar(A[tot](lambda)), (AU))"
  }  else {
    s.A.label <- expression(Spectral~~absorbance~~A(lambda)~~(AU))
    A.label.total  <- "atop(A, (AU %*% nm))"
    A.label.avg  <- "atop(bar(A(lambda)), (AU))"
  }
  if (label.qty == "total") {
    A.label <- A.label.total
  } else if (label.qty %in% c("average", "mean")) {
    A.label <- A.label.avg
  } else if (label.qty == "contribution") {
    A.label <- "atop(Contribution~~to~~total, A~~(fraction))"
  } else if (label.qty == "contribution.pc") {
    A.label <- "atop(Contribution~~to~~total, A~~(percent))"
  } else if (label.qty == "relative") {
    A.label <- "atop(Relative~~to~~sum, A~~(fraction))"
  } else if (label.qty == "relative.pc") {
    A.label <- "atop(Relative~~to~~sum, A~~(percent))"
  } else {
    A.label <- ""
  }

  y.max <- max(spct[["A"]], na.rm = TRUE)
  y.min <- 0
  plot <- ggplot(spct, aes_(~w.length, ~A)) +
    scale_fill_identity() + scale_color_identity()
  plot <- plot + geom_line()
  plot <- plot + labs(x="Wavelength (nm)", y=s.A.label)

  plot <- plot + decoration(w.band = w.band,
                            y.max = y.max,
                            y.min = y.min,
                            x.max = max(spct),
                            x.min = min(spct),
                            annotations = annotations,
                            label.qty = label.qty,
                            summary.label = A.label)

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

#' Plot a reflector spectrum
#'
#' This function returns a ggplot object with an annotated plot of a
#' reflector_spct object.
#'
#' @note Note that scales are expanded so as to make space for the annotations.
#'   The object returned is a ggplot objects, and can be further manipulated.
#'
#' @param spct a filter_spct object
#' @param w.band list of waveband objects
#' @param range an R object on which range() returns a vector of length 2, with
#'   min annd max wavelengths (nm)
#' @param pc.out logical, if TRUE use percents instead of fraction of one
#' @param label.qty character string giving the type of summary quantity to use
#'   for labels
#' @param annotations a character vector
#' @param ... other arguments passed to reflectance()
#'
#' @return a \code{ggplot} object.
#'
#' @keywords internal
#'
R_plot <- function(spct,
                   w.band,
                   range,
                   pc.out,
                   label.qty,
                   annotations,
                   ...) {
  if (!is.reflector_spct(spct)) {
    stop("R_plot() can only plot reflector_spct objects.")
  }
  if (!is.null(range)) {
    trim_spct(spct, range = range, byref = TRUE)
  }
  Rfr.type <- getRfrType(spct)
  if (length(Rfr.type)==0) {
    Rfr.type <- "unknown"
  }
  if (!pc.out) {
    scale.factor <- 1
    if (Rfr.type == "specular") {
      s.Rfr.label <- expression(Specular~~spectral~~reflectance~~R[spc](lambda)~~(fraction))
      Rfr.label.total  <- "atop(R[spc], (fraction))"
      Rfr.label.avg  <- "atop(bar(R[spc](lambda)), (fraction))"
    } else if (Rfr.type == "total") {
      s.Rfr.label <- expression(Total~~spectral~~reflectance~~R[tot](lambda)~~(fraction))
      Rfr.label.total  <- "atop(R[tot], (total))"
      Rfr.label.avg  <- "atop(bar(R[tot](lambda)), (fraction))"
    }  else {
      s.Rfr.label <- expression(Total~~spectral~~reflectance~~R(lambda)~~(fraction))
      Rfr.label.total  <- "atop(R, (total))"
      Rfr.label.avg  <- "atop(bar(R(lambda)), (fraction))"
    }
  } else if (pc.out) {
    scale.factor <- 100
    if (Rfr.type == "specular") {
      s.Rfr.label <- expression(Specular~~spectral~~reflectance~~R[spc](lambda)~~(percent))
      Rfr.label.total  <- "atop(R[spc], (total %*% 100))"
      Rfr.label.avg  <- "atop(bar(R[spc](lambda)), (percent))"
    } else if (Rfr.type == "total") {
      s.Rfr.label <- expression(Total~~spectral~~reflectance~~R[tot](lambda)~~(percent))
      Rfr.label.total  <- "atop(R[tot], (total %*% 100))"
      Rfr.label.avg  <- "atop(bar(R[tot](lambda)), (percent))"
    }  else {
      s.Rfr.label <- expression(Total~~spectral~~reflectance~~R(lambda)~~(percent))
      Rfr.label.total  <- "atop(R, (total  %*% 100))"
      Rfr.label.avg  <- "atop(bar(R(lambda)), (percent))"
    }
  }
  if (label.qty == "total") {
    Rfr.label <- Rfr.label.total
  } else if (label.qty %in% c("average", "mean")) {
    Rfr.label <- Rfr.label.avg
  } else if (label.qty == "contribution") {
    Rfr.label <- "atop(Contribution~~to~~total, R~~(fraction))"
  } else if (label.qty == "contribution.pc") {
    Rfr.label <- "atop(Contribution~~to~~total, R~~(percent))"
  } else if (label.qty == "relative") {
    Rfr.label <- "atop(Relative~~to~~sum, R~~(fraction))"
  } else if (label.qty == "relative.pc") {
    Rfr.label <- "atop(Relative~~to~~sum, R~~(percent))"
  } else {
    Rfr.label <- ""
  }
  y.max <- 1
  y.min <- 0
  plot <- ggplot(spct, aes_(~w.length, ~Rfr)) +
    scale_fill_identity() + scale_color_identity()
  plot <- plot + geom_line()
  plot <- plot + labs(x="Wavelength (nm)", y=s.Rfr.label)

  plot <- plot + decoration(w.band = w.band,
                            y.max = y.max,
                            y.min = y.min,
                            x.max = max(spct),
                            x.min = min(spct),
                            annotations = annotations,
                            label.qty = label.qty,
                            summary.label = Rfr.label)

  if (!is.null(annotations) &&
      length(intersect(c("labels", "summaries", "colour.guide"), annotations)) > 0L) {
    y.limits <- c(0, y.max * 1.25)
    x.limits <- c(min(spct) - spread(spct) * 0.025, NA)
  } else {
    y.limits <- c(0, 1)
    x.limits <- range(spct)
  }
  if (pc.out) {
    plot <- plot + scale_y_continuous(labels = scales::percent, breaks = c(0, 0.25, 0.5, 0.75, 1),
                                      limits = y.limits)
  } else {
    plot <- plot + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                                      limits = y.limits)
  }
  plot + scale_x_continuous(limits = x.limits)

}

#' Plot an object spectrum
#'
#' This function returns a ggplot object with an annotated plot of a object_spct
#' object.
#'
#' @note Note that scales are expanded so as to make space for the annotations.
#'   The object returned is a ggplot object, and can be further manipulated.
#'
#' @param spct an object_spct object
#' @param w.band list of waveband objects
#' @param range an R object on which range() returns a vector of length 2, with
#'   min annd max wavelengths (nm)
#' @param pc.out logical, if TRUE use percents instead of fraction of one
#' @param label.qty character string giving the type of summary quantity to use
#'   for labels
#' @param annotations a character vector
#' @param stacked logical
#' @param ... other arguments passed to reflectance()
#'
#' @return a \code{ggplot} object.
#'
#' @keywords internal
#'
O_plot <- function(spct,
                   w.band,
                   range,
                   pc.out,
                   label.qty,
                   annotations,
                   stacked,
                   ...) {
  if (!is.object_spct(spct)) {
    stop("O_plot() can only plot object_spct objects.")
  }
  if (!is.null(range)) {
    trim_spct(spct, range = range, byref = TRUE)
  }
   Rfr.type <- getRfrType(spct)
    if (length(Rfr.type)==0) {
    Rfr.type <- "unknown"
  }
  Tfr.type <- getTfrType(spct)
  if (length(Tfr.type)==0) {
    Tfr.type <- "unknown"
  }
  if (Rfr.type == "specular") {
    stop("Only 'total' reflectance can be plotted in a combined plot")
  }
  if (Tfr.type == "internal") {
    stop("Only 'total' transmittance can be plotted in a combined plot")
  }
  s.Rfr.label <- expression(atop(Spectral~~reflectance~R(lambda)~~spectral~~absorptance~~A(lambda), and~~spectral~~transmittance~T(lambda)))
  y.max <- 1
  y.min <- 0
  spct[["Afr"]] <- 1.0 - spct[["Tfr"]] - spct[["Rfr"]]
  tmp.spct <- dplyr::as_data_frame(spct)[c("w.length", "Tfr", "Afr", "Rfr")]
  molten.spct <- tidyr::gather_(tmp.spct, "variable", "value", c("Tfr", "Afr", "Rfr"))
  setGenericSpct(molten.spct, multiple.wl = 3L)
  plot <- ggplot(molten.spct, aes_(~w.length, ~value)) +
    scale_fill_identity()
  if (stacked) {
    plot <- plot + geom_area(aes_(alpha = ~variable), fill = "black", colour = NA)
    plot <- plot + scale_alpha_discrete(range = c(0.55, 0.25),
                                        breaks = c("Rfr", "Afr", "Tfr"),
                                        labels = c(Tfr = expression(T(lambda)),
                                                   Afr = expression(A(lambda)),
                                                   Rfr = expression(R(lambda))),
                                        guide = guide_legend(title = NULL))
  } else {
    plot <- plot + geom_line(aes_(colour = ~variable))
    plot <- plot + scale_colour_hue(labels = c(Tfr = expression(T(lambda)),
                                               Afr = expression(A(lambda)),
                                               Rfr = expression(R(lambda))),
                                    guide = guide_legend(title = NULL))
  }
  plot <- plot + labs(x="Wavelength (nm)", y=s.Rfr.label)

  annotations <- intersect(annotations, c("labels", "boxes", "segments", "colour.guide"))

  plot <- plot + decoration(w.band = w.band,
                            y.max = y.max,
                            y.min = y.min,
                            x.max = max(spct),
                            x.min = min(spct),
                            annotations = annotations,
                            label.qty = label.qty,
                            summary.label = "")
  if (!is.null(annotations) &&
      length(intersect(c("boxes", "segments", "labels", "summaries", "colour.guide"), annotations)) > 0L) {
    y.limits <- c(0, y.max * 1.25)
    x.limits <- c(min(spct) - spread(spct) * 0.025, NA)
  } else {
    y.limits <- c(0, 1)
    x.limits <- range(spct)
  }
  if (pc.out) {
    plot <- plot + scale_y_continuous(labels = scales::percent, breaks = c(0, 0.25, 0.5, 0.75, 1), limits = y.limits)
  } else {
    plot <- plot + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = y.limits)
  }

 plot + scale_x_continuous(limits = x.limits)

}

#' Plot a filter spectrum, especialization of generic plot function
#'
#' This function returns a ggplot object with an annotated plot of a filter_spct
#' object.
#'
#' @note Note that scales are expanded so as to make space for the annotations.
#'   The object returned is a ggplot objects, and can be further manipulated.
#'
#' @param x a filter_spct object
#' @param ... other arguments passed along, such as \code{label.qty}
#' @param w.band a single waveband object or a list of waveband objects
#' @param range an R object on which range() returns a vector of length 2, with
#'   min annd max wavelengths (nm)
#' @param plot.qty character string one of "transmittance" or "absorbance"
#' @param pc.out logical, if TRUE use percents instead of fraction of one
#' @param label.qty character string giving the type of summary quantity to use
#'   for labels
#' @param annotations a character vector
#'
#' @return a \code{ggplot} object.
#'
#' @export
#'
#' @keywords hplot
#'
#' @examples
#' library(photobiology)
#' plot(yellow_gel.spct)
#' plot(yellow_gel.spct, pc.out = TRUE)
#'
#' plot(polyester.spct)
#' plot(polyester.spct, plot.qty = "absorbance")
#'
#' @family plot functions
#'
plot.filter_spct <-
  function(x, ...,
           w.band = getOption("photobiology.plot.bands",
                            default = list(UVC(), UVB(), UVA(), PAR())),
           range = NULL,
           plot.qty = getOption("photobiology.filter.qty", default = "transmittance"),
           pc.out = FALSE,
           label.qty = "average",
           annotations = getOption("photobiology.plot.annotations",
                                 default = c("boxes", "labels", "summaries", "colour.guide", "peaks")) ) {
    if ("color.guide" %in% annotations) {
      annotations <- c(setdiff(annotations, "color.guide"), "colour.guide")
    }
    if (plot.qty == "transmittance") {
      out.ggplot <- T_plot(spct = x, w.band = w.band, range = range,
                           pc.out = pc.out, label.qty = label.qty,
                           annotations = annotations, ...)
    } else if (plot.qty == "absorbance") {
      out.ggplot <- A_plot(spct = x, w.band = w.band, range = range,
                           label.qty = label.qty, annotations = annotations, ...)
    } else {
      stop("Invalid 'plot.qty' argument value: '", plot.qty, "'")
    }
    if ("title" %in% annotations) {
      out.ggplot <- out.ggplot + labs(title = deparse(substitute(x)))
    }
    out.ggplot
  }

#' Plot a reflector spectrum, especialization of generic plot function
#'
#' This function returns a ggplot object with an annotated plot of a
#' reflector_spct object.
#'
#' @note Note that scales are expanded so as to make space for the annotations.
#'   The object returned is a ggplot objects, and can be further manipulated.
#'
#' @param x a reflector_spct object
#' @param ... other arguments passed along, such as \code{label.qty}
#' @param w.band a single waveband object or a list of waveband objects
#' @param range an R object on which range() returns a vector of length 2, with
#'   min annd max wavelengths (nm)
#' @param plot.qty character string (currently ignored)
#' @param pc.out logical, if TRUE use percents instead of fraction of one
#' @param label.qty character string giving the type of summary quantity to use
#'   for labels
#' @param annotations a character vector
#'
#' @return a \code{ggplot} object.
#'
#' @export
#'
#' @keywords hplot
#'
#' @examples
#'
#' library(photobiology)
#' plot(as.reflector_spct(white_body.spct))
#'
#' @family plot functions
#'
plot.reflector_spct <-
  function(x, ...,
           w.band=getOption("photobiology.plot.bands", default=list(UVC(), UVB(), UVA(), PAR())),
           range=NULL,
           plot.qty=getOption("photobiology.reflector.qty", default="reflectance"),
           pc.out=FALSE,
           label.qty="average",
           annotations=getOption("photobiology.plot.annotations",
                                 default = c("boxes", "labels", "summaries", "colour.guide", "peaks")) ) {
    if ("color.guide" %in% annotations) {
      annotations <- c(setdiff(annotations, "color.guide"), "colour.guide")
    }
    out.ggplot <- R_plot(spct = x, w.band = w.band,  range = range,
                         pc.out = pc.out, label.qty = label.qty,
                         annotations = annotations, ...)
    if ("title" %in% annotations) {
      out.ggplot <- out.ggplot + labs(title = deparse(substitute(x)))
    }
    out.ggplot
  }


#' Plot an object spectrum, especialization of generic plot function
#'
#' This function returns a ggplot object with an annotated plot of an
#' object_spct object.
#'
#' @note Note that scales are expanded so as to make space for the annotations.
#'   The object returned is a ggplot objects, and can be further manipulated.
#'
#' @param x an object_spct object
#' @param ... other arguments passed along, such as \code{label.qty}
#' @param w.band a single waveband object or a list of waveband objects
#' @param range an R object on which range() returns a vector of length 2, with
#'   min annd max wavelengths (nm)
#' @param plot.qty character string (currently ignored)
#' @param pc.out logical, if TRUE use percents instead of fraction of one
#' @param label.qty character string giving the type of summary quantity to use
#'   for labels
#' @param annotations a character vector
#' @param stacked logical
#'
#' @return a \code{ggplot} object.
#'
#' @export
#'
#' @keywords hplot
#'
#' @examples
#'
#' library(photobiology)
#' plot(black_body.spct)
#' plot(white_body.spct)
#' plot(clear_body.spct)
#'
#' @family plot functions
#'
plot.object_spct <-
  function(x, ...,
           w.band=getOption("photobiology.plot.bands", default=list(UVC(), UVB(), UVA(), PAR())),
           range=NULL,
           plot.qty=NULL,
           pc.out=FALSE,
           label.qty="average",
           annotations=getOption("photobiology.plot.annotations",
                                 default = c("boxes", "labels", "colour.guide", "peaks")),
           stacked = TRUE) {
    if ("color.guide" %in% annotations) {
      annotations <- c(setdiff(annotations, "color.guide"), "colour.guide")
    }
    out.ggplot <- O_plot(spct = x, w.band = w.band,  range = range,
                         pc.out = pc.out, label.qty = label.qty,
                         annotations = annotations, stacked = stacked, ...)
    if ("title" %in% annotations) {
      out.ggplot <- out.ggplot + labs(title = deparse(substitute(x)))
    }
    out.ggplot
  }

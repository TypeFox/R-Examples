#' Plot a response spectrum.
#'
#' This function returns a ggplot object with an annotated plot of a
#' response_spct object.
#'
#' @note Note that scales are expanded so as to make space for the annotations.
#'   The object returned is a ggplot objects, and can be further manipulated.
#'
#' @param spct a response_spct object
#' @param w.band list of waveband objects
#' @param range an R object on which range() returns a vector of length 2, with
#'   min annd max wavelengths (nm)
#' @param pc.out logical, if TRUE use percents instead of fraction of one
#' @param label.qty character string giving the type of summary quantity to use
#'   for labels
#' @param annotations a character vector
#' @param norm numeric normalization wavelength (nm) or character string "max"
#'   for normalization at the wavelength of highest peak.
#' @param ... other arguments passed to e_response()
#'
#' @return a \code{ggplot} object.
#'
#' @keywords internal
#'
e_rsp_plot <- function(spct,
                       w.band,
                       range,
                       pc.out,
                       label.qty,
                       annotations,
                       norm,
                       ...) {
  if (!is.response_spct(spct)) {
    stop("e_Rsp_plot() can only plot response_spct objects.")
  }
  q2e(spct, action="replace", byref=TRUE)
  if (!is.null(range)) {
    trim_spct(spct, range = range, byref = TRUE)
  }

  exposure.label <- NA
  if (is.null(pc.out) || is.null(norm)) {
    # no rescaling needed
    if (is_normalized(spct) || is_scaled(spct)) {
      s.rsp.label <- expression(Spectral~~energy~~response~~R[E](lambda)~~(relative~~units))
      rsp.label.total <- "atop(R[E], (relative~~units))"
      rsp.label.avg <- "atop(bar(R[E](lambda)), (relative~~units))"
      scale.factor <- 1
    } else {
      time.unit <- getTimeUnit(spct)
      if (!length(time.unit)) {
        time.unit <- "unkonwn"
      }
      if (time.unit=="second" || time.unit == lubridate::duration(1, "seconds")) {
        s.rsp.label <- expression(Spectral~~energy~~response~~R[E](lambda)~~(resp.~~unit~~s^{-1}~nm^{-1}))
        rsp.label.total  <- "atop(R[E], (resp.~~unit~~s^{-1}))"
        rsp.label.avg  <- "atop(bar(R[E](lambda)), (resp.~~unit~~s^{-1}~nm^{-1}))"
        scale.factor <- 1
      } else if (time.unit=="day" || time.unit == lubridate::duration(1, "days")) {
        s.rsp.label <- expression(Spectral~~energy~~response~~R[E](lambda)~~(resp.~~unit~~d^{-1}~nm^{-1}))
        rsp.label.total  <- "atop(R[E], (resp.~~unit~~d^{-1}))"
        rsp.label.avg  <- "atop(bar(R[E](lambda)), (resp.~~unit~~d^{-1}~nm^{-1}))"
        scale.factor <- 1
      } else if (time.unit=="hour" || time.unit == lubridate::duration(1, "hours")) {
        s.rsp.label <- expression(Spectral~~energy~~response~~R[E](lambda)~~(resp.~~unit~~h^{-1}~nm^{-1}))
        rsp.label.total  <- "atop(R[E], (resp.~~unit~~h^{-1}))"
        rsp.label.avg  <- "atop(bar(R[E](lambda)), (resp.~~unit~~h^{-1}~nm^{-1}))"
        scale.factor <- 1
      } else if (time.unit=="exposure" || lubridate::is.duration(time.unit)) {
        s.rsp.label <- expression(Spectral~~energy~~response~~R[E](lambda)~~(resp.~~unit~nm^{-1}))
        rsp.label.total  <- "atop(R[E], (resp.~~unit))"
        rsp.label.avg  <- "atop(bar(R[E](lambda)), (resp.~~unit~nm^{-1}))"
        exposure.label <- paste("Length of time:",
                                ifelse(lubridate::is.duration(time.unit),
                                       as.character(time.unit), "unknown"))
        scale.factor <- 1
      } else {
        s.rsp.label <- expression(Spectral~~energy~~response~~R[E](lambda)~~(arbitrary~~units))
        rsp.label.total <- "atop(R[E], (arbitrary~~units))"
        rsp.label.avg <- "atop(bar(R[E](lambda)), (arbitrary~~units))"
        scale.factor <- 1
      }
    }
  } else {
    # rescaling needed
    if (!is.null(norm)) {
      if (is.character(norm)) {
        if (norm %in% c("max", "maximum")) {
          idx <- which.max(spct[["s.e.response"]])
        } else if (norm %in% c("min", "minimum")) {
          idx <- which.min(spct[["s.e.response"]])
        } else if (norm %in% c("mean", "average")) {
          summary.value <- average_spct(spct)
          idx <- "summary"
        } else if (norm %in% c("total", "integral")) {
          summary.value <- integrate_spct(spct)
          idx <- "summary"
        } else {
          warning("Invalid character '", norm, "'value in 'norm'")
          return(ggplot())
        }
        if (idx == "summary") {
          scale.factor <- 1 / summary.value
        } else {
          scale.factor <- 1 / spct[idx, "s.e.response"]
          norm <- spct[idx, "w.length"]
        }
      } else if (is.numeric(norm) && norm >= min(spct) && norm <= max(spct)) {
        scale.factor <- 1 / interpolate_spct(spct, norm)[["s.e.response"]]
      } else if (is.numeric(norm)) {
        warning("'norm = ", norm, "' value outside spectral data range of ",
                min(spct), " to ", round(max(spct)), " (nm)")
        return(ggplot())
      } else {
        stop("'norm' should be numeric or character")
      }
    }
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
    s.rsp.label <-
      bquote(Spectral~~energy~~response~~R[E]( italic(lambda) )/R[E]( .(norm))~~(.(multiplier.label)))
    rsp.label.total  <- bquote(atop(integral(R[E](lambda), min, max), (.(multiplier.label))))
    rsp.label.avg  <- bquote(atop(bar(R[E](lambda)/R[E](lambda = norm)), (.(multiplier.label))))
  }
  spct[["s.e.response"]] <- spct[["s.e.response"]] * scale.factor
  y.max <- max(spct[["s.e.response"]], na.rm = TRUE)
  y.min <- 0

  if (label.qty == "total") {
    rsp.label <- "integral(R[E](lambda))"
  } else if (label.qty %in% c("average", "mean")) {
    rsp.label <- "bar(R[E](lambda))"
  } else if (label.qty == "contribution") {
    rsp.label <- "atop(Contribution~~to~~total, R[E]~~(fraction))"
  } else if (label.qty == "contribution.pc") {
    rsp.label <- "atop(Contribution~~to~~total, R[E]~~(percent))"
  } else if (label.qty == "relative") {
    rsp.label <- "atop(Relative~~to~~sum, R[E]~~(fraction))"
  } else if (label.qty == "relative.pc") {
    rsp.label <- "atop(Relative~~to~~sum, R[E]~~(percent))"
  } else {
    rsp.label <- ""
  }

  plot <- ggplot(spct, aes_(~w.length, ~s.e.response)) +
    scale_fill_identity() + scale_color_identity()
  plot <- plot + geom_line()
  plot <- plot + labs(x = "Wavelength (nm)", y = s.rsp.label)

  plot <- plot + decoration(w.band = w.band,
                            y.max = y.max,
                            y.min = y.min,
                            x.max = max(spct),
                            x.min = min(spct),
                            annotations = annotations,
                            label.qty = label.qty,
                            summary.label = rsp.label)

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

#' Plot a response spectrum.
#'
#' This function returns a ggplot object with an annotated plot of a
#' response_spct object.
#'
#' @note Note that scales are expanded so as to make space for the annotations.
#'   The object returned is a ggplot objects, and can be further manipulated.
#'
#' @param spct a response_spct object
#' @param w.band list of waveband objects
#' @param range an R object on which range() returns a vector of length 2, with
#'   min annd max wavelengths (nm)
#' @param pc.out logical, if TRUE use percents instead of fraction of one
#' @param label.qty character string giving the type of summary quantity to use
#'   for labels
#' @param annotations a character vector
#' @param norm numeric normalization wavelength (nm) or character string "max"
#'   for normalization at the wavelength of highest peak.
#' @param ... other arguments passed to q_response()
#'
#' @return a \code{ggplot} object.
#'
#' @keywords internal
#'
q_rsp_plot <- function(spct,
                       w.band,
                       range,
                       pc.out,
                       label.qty,
                       annotations,
                       norm,
                       ...) {
  if (!is.response_spct(spct)) {
    stop("q_Rsp_plot() can only plot response_spct objects.")
  }
  e2q(spct, action="replace", byref=TRUE)
  if (!is.null(range)) {
    trim_spct(spct, range = range, byref = TRUE)
  }

  exposure.label <- NA
  if (is.null(pc.out) || is.null(norm)) {
    # no rescaling needed
    if (is_normalized(spct) || is_scaled(spct)) {
      s.rsp.label <- expression(Spectral~~photon~~response~~R[Q](lambda)~~(relative~~units))
      rsp.label.total <- "atop(R[Q], (relative~~units))"
      rsp.label.avg <- "atop(bar(R[Q](lambda)), (relative~~units))"
      scale.factor <- 1
    } else {
      time.unit <- getTimeUnit(spct)
      if (!length(time.unit)) {
        time.unit <- "unkonwn"
      }
      if (time.unit=="second" || time.unit == lubridate::duration(1, "seconds")) {
        s.rsp.label <- expression(Spectral~~photon~~response~~R[Q](lambda)~~(resp.~~unit~~s^{-1}~nm^{-1}))
        rsp.label.total  <- "atop(R[Q], (resp.~~unit~~s^{-1}))"
        rsp.label.avg  <- "atop(bar(R[Q](lambda)), (resp.~~unit~~s^{-1}~nm^{-1}))"
        scale.factor <- 1
      } else if (time.unit=="day" || time.unit == lubridate::duration(1, "days")) {
        s.rsp.label <- expression(Spectral~~photon~~response~~R[Q](lambda)~~(resp.~~unit~~d^{-1}~nm^{-1}))
        rsp.label.total  <- "atop(R[Q], (resp.~~unit~~d^{-1}))"
        rsp.label.avg  <- "atop(bar(R[Q](lambda)), (resp.~~unit~~d^{-1}~nm^{-1}))"
        scale.factor <- 1
      } else if (time.unit=="hour" || time.unit == lubridate::duration(1, "hours")) {
        s.rsp.label <- expression(Spectral~~photon~~response~~R[Q](lambda)~~(resp.~~unit~~h^{-1}~nm^{-1}))
        rsp.label.total  <- "atop(R[Q], (resp.~~unit~~h^{-1}))"
        rsp.label.avg  <- "atop(bar(R[Q](lambda)), (resp.~~unit~~h^{-1}~nm^{-1}))"
        scale.factor <- 1
      } else if (time.unit=="exposure" || lubridate::is.duration(time.unit)) {
        s.rsp.label <- expression(Spectral~~photon~~response~~R[Q](lambda)~~(resp.~~unit~nm^{-1}))
        rsp.label.total  <- "atop(R[Q], (resp.~~unit))"
        rsp.label.avg  <- "atop(bar(R[Q](lambda)), (resp.~~unit~nm^{-1}))"
        exposure.label <- paste("Length of time:",
                                ifelse(lubridate::is.duration(time.unit),
                                       as.character(time.unit), "unknown"))
        scale.factor <- 1
      } else {
        s.rsp.label <- expression(Spectral~~photon~~response~~R[Q](lambda)~~(arbitrary~~units))
        rsp.label.total <- "atop(R[Q], (arbitrary~~units))"
        rsp.label.avg <- "atop(bar(R[Q](lambda)), (arbitrary~~units))"
        scale.factor <- 1
      }
    }
  } else {
    # rescaling needed
    if (!is.null(norm)) {
      if (is.character(norm)) {
        if (norm %in% c("max", "maximum")) {
          idx <- which.max(spct[["s.q.response"]])
        } else if (norm %in% c("min", "minimum")) {
          idx <- which.min(spct[["s.q.response"]])
        } else if (norm %in% c("mean", "average")) {
          summary.value <- average_spct(spct)
          idx <- "summary"
        } else if (norm %in% c("total", "integral")) {
          summary.value <- integrate_spct(spct)
          idx <- "summary"
        } else {
          warning("Invalid character '", norm, "'value in 'norm'")
          return(ggplot())
        }
        if (idx == "summary") {
          scale.factor <- 1 / summary.value
        } else {
          scale.factor <- 1 / spct[idx, "s.q.response"]
          norm <- spct[idx, "w.length"]
        }
      } else if (is.numeric(norm) && norm >= min(spct) && norm <= max(spct)) {
        scale.factor <- 1 / interpolate_spct(spct, norm)$s.q.response
      } else if (is.numeric(norm)) {
        warning("'norm = ", norm, "' value outside spectral data range of ",
                min(spct), " to ", round(max(spct)), " (nm)")
        return(ggplot())
      } else {
        stop("'norm' should be numeric or character")
      }
    }
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
    s.rsp.label <-
      bquote(Spectral~~photon~~response~~R[Q]( italic(lambda) )/R[Q]( .(norm))~~(.(multiplier.label)))
    rsp.label.total  <- bquote(atop(integral(R[Q](lambda), min, max), (.(multiplier.label))))
    rsp.label.avg  <- bquote(atop(bar(R[Q](lambda)/R[Q](lambda = norm)), (.(multiplier.label))))
  }
  spct[["s.q.response"]] <- spct[["s.q.response"]] * scale.factor
  y.max <- max(spct[["s.q.response"]], na.rm = TRUE)
  y.min <- 0

  if (label.qty == "total") {
    rsp.label <- "integral(R[Q](lambda))"
  } else if (label.qty %in% c("average", "mean")) {
    rsp.label <- "bar(R[Q](lambda))"
  } else if (label.qty == "contribution") {
    rsp.label <- "atop(Contribution~~to~~total, R[Q]~~(fraction))"
  } else if (label.qty == "contribution.pc") {
    rsp.label <- "atop(Contribution~~to~~total, R[Q]~~(percent))"
  } else if (label.qty == "relative") {
    rsp.label <- "atop(Relative~~to~~sum, R[Q]~~(fraction))"
  } else if (label.qty == "relative.pc") {
    rsp.label <- "atop(Relative~~to~~sum, R[Q]~~(percent))"
  } else {
    rsp.label <- ""
  }

  plot <- ggplot(spct, aes_(~w.length, ~s.q.response)) +
    scale_fill_identity() + scale_color_identity()
  plot <- plot + geom_line()
  plot <- plot + labs(x = "Wavelength (nm)", y = s.rsp.label)
  plot <- plot + decoration(w.band = w.band,
                            y.max = y.max,
                            y.min = y.min,
                            x.max = max(spct),
                            x.min = min(spct),
                            annotations = annotations,
                            label.qty = label.qty,
                            summary.label = rsp.label)

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



#' Plot a response spectrum, especialization of generic plot function.
#'
#' This function returns a ggplot object with an annotated plot of a
#' response_spct object.
#'
#' @note Note that scales are expanded so as to make space for the annotations.
#'   The object returned is a ggplot objects, and can be further manipulated.
#'
#' @param x a response_spct object
#' @param ... other arguments passed along, such as \code{label.qty}
#' @param w.band a single waveband object or a list of waveband objects
#' @param range an R object on which range() returns a vector of length 2, with
#'   min annd max wavelengths (nm)
#' @param unit.out character string indicating type of radiation units to use
#'   for plotting: "photon" or its synomin "quantum", or "energy"
#' @param pc.out logical, if TRUE use percents instead of fraction of one
#' @param label.qty character string giving the type of summary quantity to use
#'   for labels
#' @param annotations a character vector
#' @param norm numeric normalization wavelength (nm) or character string "max"
#'   for normalization at the wavelength of highest peak.
#'
#' @return a \code{ggplot} object.
#'
#' @export
#'
#' @keywords hplot
#'
#' @examples
#' library(photobiology)
#' plot(photodiode.spct)
#' plot(photodiode.spct, unit.out = "photon")
#'
#'
#' @family plot functions
#'
plot.response_spct <-
  function(x, ...,
           w.band=getOption("photobiology.plot.bands", default=list(UVC(), UVB(), UVA(), PAR())),
           range=NULL,
           unit.out=getOption("photobiology.radiation.unit", default="energy"),
           pc.out=FALSE,
           label.qty="total",
           annotations=getOption("photobiology.plot.annotations",
                                 default = c("boxes", "labels", "summaries", "colour.guide", "peaks")),
           norm = "max" ) {
    if ("color.guide" %in% annotations) {
      annotations <- c(setdiff(annotations, "color.guide"), "colour.guide")
    }
    if (unit.out=="photon" || unit.out=="quantum") {
      out.ggplot <- q_rsp_plot(spct=x, w.band=w.band, range=range,
                               pc.out=pc.out, label.qty=label.qty, annotations=annotations, norm = norm, ...)
    } else if (unit.out=="energy") {
      out.ggplot <- e_rsp_plot(spct=x, w.band=w.band, range=range,
                               pc.out=pc.out, label.qty=label.qty, annotations=annotations, norm = norm, ...)
    } else {
      stop("Invalid 'unit.out' argument value: '", unit.out, "'")
    }
    if ("title" %in% annotations) {
      out.ggplot <- out.ggplot + labs(title = deparse(substitute(x)))
    }
    out.ggplot
  }



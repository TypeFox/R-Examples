#' Side-by-side LV boxplots with ggplot2.
#'
#' An extension of standard boxplots which draws k letter statistics.
#' Conventional boxplots (Tukey 1977) are useful displays for conveying rough
#' information about the central 50\% of the data and the extent of the data.
#' For moderate-sized data sets (\eqn{n < 1000}), detailed estimates of tail
#' behavior beyond the quartiles may not be trustworthy, so the information
#' provided by boxplots is appropriately somewhat vague beyond the quartiles,
#' and the expected number of ``outliers'' and ``far-out'' values for a
#' Gaussian sample of size \eqn{n} is often less than 10 (Hoaglin, Iglewicz,
#' and Tukey 1986). Large data sets (\eqn{n \approx 10,000-100,000}) afford
#' more precise estimates of quantiles in the tails beyond the quartiles and
#' also can be expected to present a large number of ``outliers'' (about
#' \eqn{0.4 + 0.007 n}).
#' The letter-value box plot addresses both these shortcomings: it conveys
#' more detailed information in the tails using letter values, only out to the
#' depths where the letter values are reliable estimates of their
#' corresponding quantiles (corresponding to tail areas of roughly
#' \eqn{2^{-i}}); ``outliers'' are defined as a function of the most extreme
#' letter value shown. All aspects shown on the letter-value boxplot are
#' actual observations, thus remaining faithful to the principles that
#' governed Tukey's original boxplot.
#'
#'
#' @seealso \code{\link{stat_quantile}} to view quantiles conditioned on a
#'   continuous variable.
#' @inheritParams ggplot2::geom_point
#' @param geom,stat Use to override the default connection between
#'   \code{geom_lv} and \code{stat_lv}.
#' @param outlier.colour Override aesthetics used for the outliers. Defaults
#'   come from \code{geom_point()}.
#' @param outlier.shape Override aesthetics used for the outliers. Defaults
#'   come from \code{geom_point()}.
#' @param outlier.size Override aesthetics used for the outliers. Defaults
#'   come from \code{geom_point()}.
#' @param outlier.stroke Override aesthetics used for the outliers. Defaults
#'   come from \code{geom_point()}.
#' @param varwidth if \code{FALSE} (default) draw boxes that are the same size for each group. If
#'   \code{TRUE}, boxes are drawn with widths proportional to the
#'   square-roots of the number of observations in the groups (possibly
#'   weighted, using the \code{weight} aesthetic).
#' @param width.method character, one of 'linear' (default), 'area', or 'height'. This parameter
#' determines whether the width of the box for letter value LV(i) should be proportional to i (linear), proportional to $2^{-i}$ (height), or  whether
#' the area of the box should be proportional to $2^{-i}$ (area).
#' @export
#' @references McGill, R., Tukey, J. W. and Larsen, W. A. (1978) Variations of
#'     box plots. The American Statistician 32, 12-16.
#' @examples
#' library(ggplot2)
#' p <- ggplot(mpg, aes(class, hwy))
#' p + geom_lv(aes(fill=..LV..)) + scale_fill_brewer()
#' p + geom_lv() + geom_jitter(width = 0.2)
#' p + geom_lv(alpha=1, aes(fill=..LV..)) + scale_fill_lv()
#'
#' # Outliers
#' p + geom_lv(varwidth = TRUE, aes(fill=..LV..)) + scale_fill_lv()
#' p + geom_lv(fill = "grey80", colour = "black")
#' p + geom_lv(outlier.colour = "red", outlier.shape = 1)
#'
#' # Plots are automatically dodged when any aesthetic is a factor
#' p + geom_lv(aes(fill = drv))
#'
#' # varwidth adjusts the width of the boxes according to the number of observations
#' ggplot(ontime, aes(UniqueCarrier, TaxiIn + TaxiOut)) +
#'   geom_lv(aes(fill = ..LV..), varwidth=TRUE) +
#'   scale_fill_lv() +
#'   scale_y_sqrt() +
#'   theme_bw()
#'
#' ontime$DayOfWeek <- as.POSIXlt(ontime$FlightDate)$wday
#' ggplot(ontime, aes(factor(DayOfWeek), TaxiIn + TaxiOut)) +
#'   geom_lv(aes(fill = ..LV..)) +
#'   scale_fill_lv() +
#'   scale_y_sqrt() +
#'   theme_bw()
geom_lv <- function(mapping = NULL, data = NULL, stat = "lv",
  position = "dodge", outlier.colour = "black", outlier.shape = 19,
  outlier.size = 1.5, outlier.stroke = 0.5, na.rm = TRUE,
  varwidth = FALSE, width.method = "linear", show.legend = NA, inherit.aes = TRUE, ...)
{
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomLv,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      outlier.colour = outlier.colour,
      outlier.shape = outlier.shape,
      outlier.size = outlier.size,
      outlier.stroke = outlier.stroke,
      varwidth = varwidth,
      width.method = width.method,
      ...
    )
  )
}

#' @rdname geom_lv
#' @importFrom grid grobTree
#' @export
GeomLv <- ggplot2::ggproto("GeomLv", ggplot2::Geom,
  setup_data = function(data, params) {
  #  browser()
    data$width <- data$width %||%
      params$width %||% (resolution(data$x, FALSE) * 0.9)

    if (!is.null(data$outliers)) {
      suppressWarnings({
        out_min <- vapply(data$outliers, min, numeric(1))
        out_max <- vapply(data$outliers, max, numeric(1))
      })

      data$ymin_final <- pmin(out_min, data$ymin)
      data$ymax_final <- pmax(out_max, data$ymax)
    }

    # if `varwidth` not requested or not available, don't use it
    if (is.null(params) || is.null(params$varwidth) || !params$varwidth || is.null(data$relvarwidth)) {
      data$xmin <- data$x - data$width / 2
      data$xmax <- data$x + data$width / 2
    } else {
      # make `relvarwidth` relative to the size of the largest group
      data$relvarwidth <- data$relvarwidth / max(data$relvarwidth)
      data$xmin <- data$x - data$relvarwidth * data$width / 2
      data$xmax <- data$x + data$relvarwidth * data$width / 2
    }
    data$width <- NULL
# don't delete the relative width for LV plots
#    if (!is.null(data$relvarwidth)) data$relvarwidth <- NULL

    data
  },

  draw_group = function(data, panel_scales, coord,
                        outlier.colour = "black", outlier.shape = 19,
                        outlier.size = 1.5, outlier.stroke = 0.5,
                        width.method="linear",
                        varwidth = FALSE) {
    common <- data.frame(
      colour = data$colour,
      size = data$size,
      linetype = data$linetype,
      fill = alpha(data$fill, data$alpha),
      group = data$group,
      stringsAsFactors = FALSE
    )

    i <- seq_len(data$k[1]-1)-1
    data$width <- data$xmax - data$xmin

    lower <- rev(seq_len(data$k[1]-1)) +1
    upper <- seq_len(data$k[1]-1)+1

    if (width.method=="linear") {
      offset <- c(0, (i / data$k[1]))  * data$width/2
    } else {
      if (width.method=="height") {
        height <- 2^(-c(1, as.numeric(data$LV)[-nrow(data)])+1)
        offset <- (1-height)*data$width/2
      } else {
        if (width.method=="area") {
#          browser()
          areas <- 2^(-as.numeric(data$LV))
          lheight <- c(0, -diff(data$lower))
          uheight <- c(0, diff(data$upper))

          offset <- (1- areas/lheight)*data$width/2
          offset[is.infinite(offset)] <- data$width[1]/2
        }
      }
    }
    # boxes for lower letter values
    # bottom rectangles:
    lowbox <- data.frame(
      xmin = data$xmin[lower] + offset[lower],
      xmax = data$xmax[lower] - offset[lower],
      ymin = data$lower[lower-1],
      ymax = data$lower[lower],
      alpha = data$alpha[lower],
      LV = data$LV[lower],
      common[lower,],
      stringsAsFactors = FALSE
    )

    if (width.method == "area") {
      offset <- (1-areas/uheight)*data$width/2
      offset[is.infinite(offset)] <- data$width[1]/2
    }

    # top rectangles:
    hibox <- data.frame(
      xmin = data$xmin[upper] + offset[upper],
      xmax = data$xmax[upper] - offset[upper],
      ymin = data$upper[upper-1],
      ymax = data$upper[upper],
      alpha = data$alpha[upper],
      LV = data$LV[upper],
      common[upper,],
      stringsAsFactors = FALSE
    )
    # medians, not rectangles:
    medians <- data.frame(
      xmin = data$xmin[1],
      xmax = data$xmax[1],
      ymin = data$upper[1],
      ymax = data$upper[1],
      alpha = data$alpha[1],
      LV = data$LV[1],
      common[1,],
      stringsAsFactors = FALSE
    )
    box <- rbind(medians, lowbox, hibox)
#browser()
    medians <- subset(box, LV=="M")
    medians <-   transform(medians,
        colour = fill,
        x = xmin,
        xend = xmax,
        y = ymin,
        yend = ymax
      )

    if (!is.null(data$outliers) && length(data$outliers[[1]] >= 1)) {
      outliers <- data.frame(
        y = data$outliers[[1]],
        x = data$x[1],
        colour = outlier.colour %||% data$colour[1],
        shape = outlier.shape %||% data$shape[1],
        size = outlier.size %||% data$size[1],
        stroke = outlier.stroke %||% data$stroke[1],
        fill = NA,
        alpha = NA,
        stringsAsFactors = FALSE
      )
      outliers_grob <- GeomPoint$draw_panel(outliers, panel_scales, coord)
    } else {
      outliers_grob <- NULL
    }

    ggplot2:::ggname("geom_lv", grobTree(
      outliers_grob,
      GeomRect$draw_panel(box, panel_scales, coord),
      GeomSegment$draw_panel(medians, panel_scales, coord)
    ))
  },

  draw_key = ggplot2::draw_key_rect,

  default_aes = ggplot2::aes(weight = 1, colour = "grey70", fill = "grey70", size = 0.5,
    alpha = 1, shape = 19, linetype = "solid", outlier.colour = "black",
    outlier.shape = 19, outlier.size = 1.5, outlier.stroke = 0.5),

  required_aes = c("x", "k", "LV")
)

#' @export
#' @rdname geom_lv
scale_fill_lv <- function(...) {
  greys <- rev(RColorBrewer::brewer.pal(9, "Greys"))
  oranges <- RColorBrewer::brewer.pal(3, "Oranges")
  colors <- c("white", greys[2:4], oranges[3], greys[4:6], oranges[2],
      greys[6:8], oranges[1], greys[8:9])

  ggplot2::scale_fill_manual(..., values = colors)
}

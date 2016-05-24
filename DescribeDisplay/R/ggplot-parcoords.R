#' Compact pcp data
#' A parallel coordinates is written out as a series of 1D dotplots.  This function
#' compacts it back into one dataset.
#'
#' @param data data to pull points from
#' @author Hadley Wickham \email{h.wickham@@gmail.com}
#' @keywords internal
#' @importFrom plyr ldply
compact_pcp <- function(data) {
  ldply(data$plots, function(p) {
    data.frame(
      id = 1:nrow(p$points),
      variable = p$params$label,
      p$points[c("col", "pch", "cex")],
      x = p$points$y %||% 1,
      y = p$points$x
    )
  })
}

#' Scale the values by range
#' Divide the values of the objects by the range of values
#'
#' @param x values to be worked on
#' @author Hadley Wickham \email{h.wickham@@gmail.com}
#' @keywords internal
#' @export
range01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) {
    rep(0, length(x))
  } else {
    (x - rng[1]) / diff(rng)
  }
}

#' Create a nice plot for parallel coordinates plot
#' Create a nice looking plot complete with axes using ggplot.
#'
#' @param data plot to display
#' @param absoluteX make the sections proportional horizontally to eachother
#' @param absoluteY make the sections proportional vertically to eachother
#' @param edges boolean value to print edges.  Defaults to TRUE.
#' @param ... arguments passed to the grob function
#' @author Barret Schloerke \email{schloerke@@gmail.com}
#' @keywords hplot
#' @export
#' @importFrom plyr ddply
#' @examples
#' library(ggplot2)
#' print(ggplot(dd_example("parcoord")))
ggplot.parcoords <- function(
  data,
  absoluteX = FALSE,
  absoluteY = FALSE,
  edges = TRUE,
  ...
) {
  variable <- NULL
  x <- NULL
  y <- NULL
  id <- NULL
  df <- compact_pcp(data)

  if (absoluteX) {
    std <- transform(df, x = as.numeric(variable) + range01(x) / 2)
  } else {
    # Scale variables individually
    std <- ddply(df, "variable", transform,
      x = as.numeric(variable) + range01(x) / 2)
  }

  if (!absoluteY) {
    std <- ddply(std, "variable", transform, y = range01(y))
    # yscale <- scale_y_continuous(NULL)
  } else {
    # ybreaks <- seq(min(std$y), max(std$y), length = 5)
    # ylabels <- seq(min(df$y), max(df$y), length = 5)

    # yscale <- scale_y_continuous(breaks = ybreaks, labels = ylabels)
  }

  vars <- levels(df$variable)

  ### Make a pretty picture
  p <- ggplot(std, aes(x, y, group = id, colour = col, order = col)) +
    scale_colour_identity() +
    scale_size_identity() +
    scale_shape_identity() +
    scale_linetype_identity() +
    theme(title = element_text(data$title)) +
    # scale_y_continuous("", breaks = ybreaks,
    #   labels = rep("", length(ybreaks))) +
    scale_x_continuous("", breaks = seq_along(vars),
      labels = vars, minor_breaks = FALSE)

  if (edges) {
    p <- p + geom_line(aes_string(size = "cex * 2"))
  }

  # Plot points on top
  if (data$showPoints) {
    p <- p + geom_point(aes_string(shape = "pch", size = "cex * 4.5"))
  }

  p
}

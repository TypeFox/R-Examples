#' Create a nice plot for Histograms
#' Create a nice looking plot complete with axes using ggplot.
#'
#' @param data plot to display, object created by \code{dd_load()}
#' @param spine (not implemented currently) whether to display the barchart as a spine plot
#' @param ... arguments passed through to the ggplot function
#' @author Barret Schloerke \email{schloerke@@gmail.com}
#' @keywords hplot
#' @export
#' @examples
#' library(ggplot2)
#' print(ggplot(dd_example("barchart")))
#' print(ggplot(dd_example("histogram")))
ggplot.histogram <- function(data, spine = FALSE, ...) {
#  cat("\nggplot.histogram\n")

  p <- ggplot(
      data$points,
      aes_string(x = "x", fill = "col", ...)
    ) +
    scale_x_continuous(data$params$label) +
    coord_flip() +
    scale_size_identity() +
    scale_shape_identity() +
    scale_linetype_identity() +
    scale_fill_identity()


#  if (spine){
#    ## not correct yet
  # p <- p + geom_bar(
  #   position = "fill",
  #   binwidth = diff(data$params$breaks[1:2]),
  #   ...
  # )
#    cat("\nspine\n")
#  }else{
    allBreaks <- c(
      data$params$breaks,
      data$params$breaks[length(data$params$breaks)] +
        diff(data$params$breaks[1:2])
    )

    p <- p + geom_histogram(breaks = allBreaks, ...)
#  }

  p
}




#' Create a nice plot for Bar Plots
#' Create a nice looking plot complete with axes using ggplot.
#'
#' @param data plot to display, object created by \code{dd_load()}
#' @param spine (not implemented currently) whether to display the barchart as a spine plot
#' @param ... arguments passed through to the ggplot function
#' @author Barret Schloerke \email{schloerke@@gmail.com}
#' @keywords hplot
#' @export
#' @examples
#' library(ggplot2)
#' print(ggplot(dd_example("barchart")))
ggplot.barplot <- function(data, spine = FALSE, ...){
#  cat("\nggplot.barplot\n")
  levelnames <- data$params$levelnames
  levelNameOrder <- data$params$levelvalues + 1
  xVals <- data$points$x
  for (i in 1:length(levelnames)){
    xVals[xVals == i] <- levelnames[levelNameOrder[i]]
  }

  data$points$splitBy <- xVals
#print(unique(data$points$splitBy))
  temp <- unique(data$points$splitBy)

  p <- ggplot(data$points, aes_string(x = "splitBy", fill = "col", ...)) +
    geom_bar() +
    coord_flip() +
    scale_x_discrete(data$params$label, limits = c(temp)) +
    scale_size_identity() +
    scale_shape_identity() +
    scale_linetype_identity() +
    scale_fill_identity() #+
#    xlim(temp)


# scale_x_continuous(
#   data$params$label,
#   limits = c(unique(data$points$splitBy))
# )+


  p
}




#' Visualization of dataset characteristics
#'
#' @param x A \code{\link{DatasetCharacterization}} object
#' @param y Ignored
#' @param lines Draw observation dependency lines
#' @param points Draw observation points
#' @param null.line Draw null line
#' @param null.line.col Null line color
#' @param basis Draw basis characterization of the dataset
#' @param basis.col Color of basis characterization
#' @param ... Ignored
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @method plot DatasetCharacterization
#'
#' @family dataset-characterization
#'
#' @rdname datachar-visualization
#'
#' @importFrom graphics plot
#' @importFrom plyr ddply
#' @S3method plot DatasetCharacterization
plot.DatasetCharacterization <- function(x, y = NULL, lines = TRUE, points = TRUE,
                                         null.line = TRUE, null.line.col = gray(0.7),
                                         basis = TRUE, basis.col = NULL, ...) {

  stopifnot(nlevels(x$datasets[, drop = TRUE]) == 1)

  x <- ddply(x, "characteristics", dcscale)
  #data <- subset(x, subset = samples != "basis")
  data <- x[x$samples != "basis", ]
  #data.basis <- subset(x, subset = samples == "basis")
  data.basis <- x[x$samples == "basis", ]
                                             
  p <- ggplot(data, aes(characteristics, value, group = samples))

  if ( null.line )
    p <- p + geom_hline(aes(yintercept = 0), colour = null.line.col)

  if ( lines )
    p <- p + geom_line()

  if ( points )
    p <- p + geom_point()

  if ( (nrow(data.basis) > 0) && basis ) {
    if ( is.null(basis.col) )
      basis.col <- default_colors(n = 1)

    p <- p + geom_line(data = data.basis,
                       aes(characteristics, value, group = samples),
                       colour = basis.col)

    p <- p + geom_point(data = data.basis,
                        aes(characteristics, value, group = samples),
                        colour = basis.col)
  }

  p <- p + scale_y_continuous('', breaks = c(-0.2, seq(0, 1, by = 0.2)),
                              labels = c("NA", seq(0, 1, by = 0.2))) +
           scale_x_discrete("Characteristics") +
           theme_update(axis.text.x = theme_text(angle = 90, hjust = 1))

  p
}



### Internal functions: ##############################################

dcscale <- function(x) {
  x$value <- dcscale0(x$value)
  x
}



dcscale0 <- function(x) {
  rx <- range(x, na.rm = TRUE)

  if ( rx[1] == rx[2] )
    return(rep(1, length = length(x)))

  sx <- (x - min(x, na.rm = TRUE)) /
      (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))

  sx[is.na(sx)] <- -0.2

  sx
}


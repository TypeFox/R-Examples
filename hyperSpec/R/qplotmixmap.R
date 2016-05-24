##' map plot with colour overlay.
##'
##' 
##' @title qplotmap with colour mixing for multivariate overlay
##' @param object hyperSpec object
##' @param ... handed over to \code{\link[hyperSpec]{qmixlegend}} and \code{\link[hyperSpec]{qmixtile}}
##' @return invisible list with ggplot2 objects map and legend
##' @seealso \code{\link[hyperSpec]{qmixtile}}
##' @author Claudia Beleites
##' @importFrom grid pushViewport viewport popViewport grid.layout unit
##' @import ggplot2
##' @export
##' @examples
##' chondro <- chondro - spc.fit.poly.below (chondro)
##' chondro <- sweep (chondro, 1, apply (chondro, 1, mean), "/")
##' chondro <- sweep (chondro, 2, apply (chondro, 2, quantile, 0.05), "-")
##' 
##' qplotmixmap (chondro [,,c (940, 1002, 1440)],
##'              purecol = c (colg = "red", Phe = "green", Lipid = "blue"))
##' 
##' qplotmixmap (chondro [,,c (940, 1002, 1440)],
##'              purecol = c (colg = "red", Phe = "green", Lipid = "blue"))
##' 
qplotmixmap <- function (object, ...){

	p <- qmixtile (object@data, ...) +
       coord_equal ()
  
  p <- p +
       xlab (labels (object)[[as.character (p$mapping$x)]]) + 
       ylab (labels (object)[[as.character (p$mapping$y)]]) 

  l <- qmixlegend (object@data$spc, ...)

  legendright (p, l)

  invisible (list (map = p, legend = l))
}

##' Plot multivariate data into colour channels
##'
##' plot graph with legend right of it
##'
##' @rdname qplotmix
##' @param p plot object
##' @param l legend object
##' @param legend.width,legend.unit size of legend part
##' @return invisible \code{NULL}
##' @author Claudia Beleites
##' @rdname qplotmix
##' @export 
legendright <- function (p, l, legend.width = 8, legend.unit = "lines") {
  plot.new ()
  pushViewport (viewport (layout = grid.layout (1, 2, 
                            widths = unit (c (1, legend.width), c("null", legend.unit))
                            )))
  print (p, viewport (layout.pos.col = 1), newpage = FALSE)
  print (l, viewport (layout.pos.col = 2), newpage = FALSE)
  popViewport ()
}
  
##' plot multivariate data into colour channels using \code{\link[ggplot2]{geom_tile}}
##' @rdname qplotmix
##' @param object matrix to be plotted with mixed colour channels
##' @param purecol pure component colours, names determine legend labels
##' @param mapping see \code{\link[ggplot2]{geom_tile}}
##' @param ... \code{qmixtile}: handed to \link[hyperSpec]{colmix.rgb}
##'
##' \code{qmixlegend} and \code{colmix.rgb} hand further arguments to the \code{normalize} function
##' @param map.tileonly if \code{TRUE}, \code{mapping} will be handed to
##' \code{\link[ggplot2]{geom_tile}} instead of \code{\link[ggplot2]{ggplot}}.
##' 
qmixtile <- function (object,
                      purecol = stop ("pure component colors needed."),
                      mapping = aes_string (x = "x", y = "y", fill = "spc"),
                      ...,
                      map.tileonly = FALSE) {

  ## calculate fill colours
  fill <- colmix.rgb (object  [[as.character (mapping$fill)]], purecol, ...)
  object [[as.character (mapping$fill)]] <- fill

  if (map.tileonly)
      p <- ggplot (object) + geom_tile (mapping = mapping, data = object)
  else
      p <- ggplot (object, mapping = mapping) + geom_tile ()

  p + scale_fill_identity () + theme (legend.position = "none")  
}

##' \code{normalize.colrange} normalizes the range of each column to [0, 1]
##' @rdname qplotmix
##' @export
##' 
##' @param na.rm see \code{link[base]{min}}
##' @param legend should a legend be produced instead of normalized values?
##' @param n of colours to produce in legend
##' @return list with components ymin, max and fill to specify value and fill colour value (still
##' numeric!) for the legend, otherwise the normalized values
normalize.colrange <- function (x, na.rm = TRUE, legend = FALSE, n = 100, ...){
  ## legend
  if (legend){
    y <- apply (x, 2, function (x) seq (min (x), max (x), length.out = n))
    dy2 <- abs (y [2,] - y [1,]) / 2

    list (ymin = sweep (y, 2, dy2, `-`),
          ymax = sweep (y, 2, dy2, `+`),
          fill = apply (x, 2, function (x) seq (      0,       1, length.out = n)))
  } else {
    ## normalized values
    x <- sweep (x, 2, apply (x, 2, min, na.rm = na.rm), `-`)
    sweep (x, 2, apply (x, 2, max, na.rm = na.rm), `/`)
  }
}

##' \code{normalize.range} normalizes the range of all columns to [0, 1]
##' @rdname qplotmix
##' @export
##' 
normalize.range <- function (x, na.rm = TRUE, legend = FALSE, n = 100, ...){
  if (legend){
    y <- matrix (seq (min (x), max (x), length.out = n), nrow = n, ncol = ncol (x))
     dy2 <- abs (y [2,] - y [1,]) / 2

    list (ymin = sweep (y, 2, dy2, `-`),
          ymax = sweep (y, 2, dy2, `+`),
          fill = apply (x, 2, function (x) seq (      0,       1, length.out = n)))
  } else {
  x <- x - min (x, na.rm = na.rm)
       x / max (x, na.rm = na.rm)
  }
}

##' \code{normalize.null} does not touch the values 
##' @rdname qplotmix
##' @export
##' 
normalize.null <- function (x, na.rm = TRUE, legend = FALSE, n = 100, ...){
  if (legend){
    y <- apply (x, 2, function (x) seq (min (x), max (x), length.out = n))

    list (ymin = sweep (y, 2, min),
          ymax = sweep (y, 2, max),
          fill = apply (x, 2, function (x) seq (      0,       1, length.out = n)))
  } else {
    x 
  }
}
##' \code{normalize.minmax} normalizes the range of each column j to [min_j, max_j]
##' @rdname qplotmix
##' @export
##' @param min numeric with value corresponding to "lowest" colour for each column
##' @param max numeric with value corresponding to "hightest" colour for each column
normalize.minmax <- function (x, min = 0, max = 1, legend = FALSE, n = 100, ...){
  if (legend){
    y <- matrix (seq (0, 1, length.out = n), nrow = n, ncol = ncol (x))
    y <- sweep (y, 2, max - min, `*`)
    y <- sweep (y, 2, min, `+`)

    dy2 <- abs (y [2,] - y [1,]) / 2

    l <- list (ymin = sweep (y, 2, dy2, `-`),
               ymax = sweep (y, 2, dy2, `+`),
               ymax = y + dy2,
               fill = matrix (seq (0, 1, length.out = n), nrow = n, ncol = ncol (x)))
  
    l$ymin [1, ] <- pmin (l$ymin [1,], apply (x, 2, min, na.rm = TRUE))
    l$ymax [n, ] <- pmax (l$ymax [n,], apply (x, 2, max, na.rm = TRUE))

    l
  } else {
  x <- sweep (x, 2, min, `-`)
       sweep (x, 2, max, `/`)
  }
}

##' legends for mixed colour plots
##' @rdname qplotmix
##' @param dx width of label bar
##' @param ny number of colours in legend
##' @param labels component names
##' @return ggplot object with legend
##' @author Claudia Beleites
##' @export 
qmixlegend <- function (x, purecol, dx = 0.33, ny = 100, labels = names (purecol),
                        normalize = normalize.colrange, ...) {
  if (! is.matrix (x))
    x <- matrix (x, ncol = 1)

  if (is.null (labels))
    labels <- colnames (x)
  if (is.null (labels))
    labels <- seq_len (ncol (x))

  if (! is.null (normalize))
      l <- normalize (x, ..., legend = TRUE)
  else
      l <- x
  
  df <- data.frame ()
  for (column in seq_along (purecol)){
    tmp <- colmix.rgb (l$fill [, column, drop = FALSE], purecol [column], normalize = NULL, ...)
    df <- rbind (df, data.frame (column = labels [column],
                                 col = tmp,
                                 ymin = l$ymin [, column],
                                 ymax = l$ymax [, column])
                 )
  }
  df$column <- as.factor (df$column)
  df$xmin <- as.numeric (df$column) - dx
  df$xmax <- as.numeric (df$column) + dx

  l <- ggplot (df, aes (x = column), col = col) +
    geom_point (aes (x=column, y = 1), col = NA) + ylab ("") + xlab ("")
  l <- l + geom_rect (aes_string (xmin = "xmin", xmax = "xmax", ymin = "ymin", ymax = "ymax",
                           fill = "col",  colour = "col"))

  l <- l + theme (plot.margin = unit(c(0.5, 0, 0 ,0), "lines"), legend.position = "none") +
       scale_fill_identity () + scale_colour_identity () 

  l
}

##' @rdname qplotmix
##' @title multi channel colour mixing
##' @param x matrix with component intensities in columns
##' @param against value to mix against (for \code{sub = TRUE} only, 1 = white, 0 = black)
##' @param sub subtractive color mixing?
##' @param normalize function to normalize the values.
##' @return character with colours
##' @author Claudia Beleites
##' @export
colmix.rgb <- function (x, purecol, against = 1, sub = TRUE,
                        normalize = normalize.colrange, ...){
  if (! is.null (normalize))
    x <- normalize (x, ...)

  if (is.character (purecol))
    purecol <- t (col2rgb (purecol)) / 255

  if (sub)
    x <- against - x %*% (against - purecol)
  else
    x <- x %*% purecol
  
  x [x < 0] <- 0
  x [x > 1] <- 1

  cols <- rep (NA, nrow (x))
  cols [! is.na (x [,1])] <-   rgb (x [!is.na (x [, 1]),])

  cols
}

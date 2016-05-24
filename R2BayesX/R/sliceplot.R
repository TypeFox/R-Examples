sliceplot <- function(x, y = NULL, z = NULL, view = 1, c.select = NULL,
  values = NULL, probs = c(0.1, 0.5, 0.9), grid = 100,
  legend = TRUE, pos = "topright", digits = 2, data = NULL,
  rawdata = FALSE, type = "akima", linear = FALSE, extrap = FALSE,
  k = 40, rug = TRUE, rug.col = NULL, jitter = TRUE, ...)
{
  if(is.vector(x) & is.vector(y) & is.vector(z)) {
    nx <- c(
      deparse(substitute(x), backtick = TRUE, width.cutoff = 500),
      deparse(substitute(y), backtick = TRUE, width.cutoff = 500),
      deparse(substitute(z), backtick = TRUE, width.cutoff = 500)
    )
    x <- cbind(x, y, z)
    colnames(x) <- nx
  } else {
    if(inherits(x,"formula")) {
      if(is.null(data))
        data <- environment(x)
      else
        if(is.matrix(data))
          data <- as.data.frame(data)
      x <- model.frame(x, data = data)
      if(ncol(x) < 3L)
        stop("formula is specified wrong!")
      if(ncol(x) > 3L)
        x <- x[, c(2L, 3L, 1L, 4L:ncol(x))]
      else
        x <- x[, c(2L, 3L, 1L)]
    }
  }
  stopifnot(is.matrix(x) || is.data.frame(x))
  nx <- colnames(x)
  if(is.null(c.select))
    c.select <- 3
  if(c.select < 3)
    c.select <- if(c.select < 2) 3 else 4 
  if(c.select > ncol(x))
    stop("column number selected in c.select is larger than the number of existing columns in x!")
  if(is.character(view))
    view <- grep(view, nx, ignore.case = TRUE)
  x <- x[order(x[, view]), ]
  noview <- if(view < 2) 2 else 1
  values <- if(is.null(values)) {
    quantile(x[, noview], probs = probs, type = 1)
  } else values
  if(!rawdata) {
    xo <- seq(min(x[, view]), max(x[, view]), length = grid)
    yo <- seq(min(x[, noview]), max(x[, noview]), length = grid)
    zi <- interp2(x[, view], x[, noview], x[, c.select],
      xo = xo,
      yo = yo,
      type = type, linear = linear, extrap = extrap, k = k)
    yg <- rep(yo, each = grid)
    zg <- as.vector(zi)
    slices <- xo
  } else {
    yg <- x[, noview]
    zg <- x[, c.select]
    slices <- unique(x[, view])
  }
  for(j in values) {
    val <- unique(yg[which.min(abs(yg - j))])
    slices <- cbind(slices, zg[yg == val])
  }
  k <- ncol(slices)
  args <- l.args <- list(...)
  args$lty <- if(is.null(args$lty)) 1:k else rep(args$lty, length.out = k)
  args$col <- if(is.null(args$col)) "black" else rep(args$col, length.out = k)
  args$lwd <- if(is.null(args$lwd)) 1 else rep(args$lwd, length.out = k)
  if(is.null(args$xlab))
    args$xlab <- nx[view]
  if(is.null(args$ylab))
    args$ylab <- paste("Effect of", nx[view])
  args$x <- slices[, 1]
  args$y <- slices[, 2:ncol(slices)]
  args$type = "l"
  do.call("matplot", delete.args("matplot", args, c("axes", "main", "xlab", "ylab")))
  if(legend) {
    l.args$x <- pos
    l.args$legend <- paste(nx[noview], "=", round(values, digits))
    l.args$lty <- args$lty
    l.args$col <- args$col
    l.args$lwd <- args$lwd
    if(is.null(l.args$bg))
      l.args$bg <- NA
    if(is.null(l.args$box.col))
      l.args$box.col <- NA
    do.call("legend", delete.args("legend", l.args))
  }
  if(rug) {
    args$x <- if(jitter) jitter(x[, view]) else x[, view]
    args$col <- rug.col
    do.call(graphics::rug, delete.args(graphics::rug, args))
  }
  invisible(args)
}


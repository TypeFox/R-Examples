plotmap <- function(map, x = NULL, id = NULL, c.select = NULL, legend = TRUE,
  missing = TRUE, swap = FALSE, range = NULL, names = FALSE, values = FALSE, col = NULL,
  ncol = 100, breaks = NULL, cex.legend = 1, cex.names = 1, cex.values = cex.names,
  digits = 2L, mar.min = 2, add = FALSE, interp = FALSE, grid = 200, land.only = FALSE,
  extrap = FALSE, outside = FALSE, type = "akima", linear = FALSE, k = 40,
  p.pch = 15, p.cex = 1, shift = NULL, trans = NULL, ...)
{
  if(missing(map))
    stop("map object is missing!")
  if(inherits(map, "SpatialPolygons"))
    map <- sp2bnd(map)
  if(!is.list(map))
    stop("argument map must be a list() of matrix polygons!")
  args <- list(...)
  map.limits <- find.limits(map, mar.min, ...)
  if(is.null(args$asp)) {
    args$asp <- attr(map, "asp")
    if(is.null(args$asp))
      args$asp <- map.limits$asp
  }
  n <- length(map)
  if(is.null(x))
    legend <- FALSE
  poly.names.orig <- names(map)
  if(!any(is.na(poly.names <- x2int(names(map))))) {
    op <- order(poly.names)
    poly.names <- poly.names[op]
    poly.names <- as.character(poly.names)
  } else {
    poly.names <- names(map)
    op <- order(poly.names)
    poly.names <- poly.names[op]
  }
  poly.names.orig <- poly.names.orig[op]
  map <- map[op]
  if(length(upn <- unique(poly.names)) < length(poly.names)) {
    nn <- NULL
    for(i in upn) {
      j <- poly.names == i
      poly.names[j] <- paste(poly.names[j],
        if(sum(j) > 1) paste(".", 1:sum(j), sep = "") else NULL,
        sep = "")
    }
    names(map) <- poly.names
  }
  map <- map[poly.names]
  poly.names <- names(map)
  surrounding <- attr(map, "surrounding")
  inner.id <- which(sapply(surrounding, length) > 0L)
  if(length(inner.id)) {
    poly.names <- c(poly.names[- inner.id], poly.names[inner.id])
    map <- c(map[- inner.id], map[inner.id])
  }
  if(!is.null(args$ylim))
    map.limits$ylim <- args$ylim
  if(!is.null(args$xlim))
    map.limits$xlim <- args$xlim
  if(is.null(args$symmetric))
    symmetric <- TRUE
  else
    symmetric <- args$symmetric
  if(!is.null(x)) {
    if(is.null(col)) {
      col <- colorspace::diverge_hcl
      # col <- colorspace::diverge_hcl(ncol, h = c(130, 10), c = 250,
      #  l = c(30, 90), power = 1.5, gamma = 2.4, fixup = TRUE)
    }
    x <- compute.x.id(x, id, c.select, range, symmetric)
    if(!is.null(shift)) {
      shift <- as.numeric(shift[1])
      x$x <- x$x + shift
    }
    if(!is.null(trans)) {
      if(!is.function(trans)) stop("argument trans must be a function!")
      x$x <- trans(x$x)
    }
    map_fun <- make_pal(col = col, ncol = ncol, data = as.numeric(x$x), 
      range = range, breaks = breaks, swap = swap, 
      symmetric = symmetric)$map
    colors <- map_fun(as.numeric(x$x))
  } else {
    if(is.null(col))
      colors <- rep(NA, length.out = n)
    else {
      if(is.function(col))
        colors <- col(ncol)
      else colors <- col
      colors <- rep(colors, length.out = n)
    }
  }
  if(is.null(args$pos))
    args$pos <- "right"
  if(legend && !is.null(args$pos) && args$pos[1L] == "right") {
    par.orig <- par(c("mar", "las", "mfrow"))
    mar.orig <- mar <- par.orig$mar
    mar[4L] <- 0
    mar[c(1, 3)] <- 1
    on.exit(par(par.orig))
    par(mar = mar)
    w <- (3 + mar[2L]) * par("csi") * 2
    w <- max(c(2.84, w))
    layout(matrix(c(1, 2), nrow = 1), widths = c(1, lcm(w)))
  }
  if(!is.null(map.limits$mar) && is.null(args$asp) && !add)
    par(mar = map.limits$mar)
  args$x <- map.limits$x
  args$y <- map.limits$y
  if(is.null(args$type))
    args$type <- "n"
  if(is.null(args$axes))
    args$axes <- FALSE
  if(is.null(args$xlab))
    args$xlab <- ""
  if(is.null(args$ylab))
    args$ylab <- ""
  if(!add)
    do.call(graphics::plot.default, delete.args(graphics::plot.default, args))
  if(interp & !is.null(x)) {
    cdata <- data.frame(centroids(map), "id" = names(map))
    cdata <- merge(cdata, data.frame("z" = x$x, "id" = x$id), by = "id")
    cdata <- unique(cdata)

    xo <- seq(map.limits$x[1], map.limits$x[2], length = grid)
    yo <- seq(map.limits$y[1], map.limits$y[2], length = grid)

    ico <- with(cdata, interp2(x = x, y = y, z = z,
      xo = xo,
      yo = yo,
      type = type, linear = linear, extrap = extrap,
      k = if(is.null(k)) ceiling(length(map) * 0.8) else as.integer(k)))
    
    yco <- rep(yo, each = length(xo))
    xco <- rep(xo, length(yo))

    d4x <- abs(diff(xco))
    d4x <- min(d4x[d4x != 0], na.rm = TRUE)
    d4y <- abs(diff(yco))
    d4y <- min(d4y[d4y != 0], na.rm = TRUE)
    res <- c(d4x, d4y)
    pp <- NULL
    if(length(res))
      pp <- cbind(xco, yco)

    cvals <- as.numeric(ico)
    cvals[cvals < min(cdata$z)] <- min(cdata$z)
    cvals[cvals > max(cdata$z)] <- max(cdata$z)
    icolors <- map_fun(cvals)

    if(!outside) {
      maptools::gpclibPermit()
      class(map) <- "bnd"
      mapsp <- bnd2sp(map)
      ob <- maptools::unionSpatialPolygons(mapsp, rep(1L, length = length(mapsp)), avoidGEOS  = TRUE)

      nob <- length(slot(slot(ob, "polygons")[[1]], "Polygons"))
      pip <- NULL
      for(j in 1:nob) {
        oco <- slot(slot(slot(ob, "polygons")[[1]], "Polygons")[[j]], "coords")
        pip <- cbind(pip, sp::point.in.polygon(xco, yco, oco[, 1L], oco[, 2L], mode.checked = FALSE) < 1L)
      }
      pip <- apply(pip, 1, function(x) all(x))
    
      icolors[pip] <- NA
    }

    if(land.only) {
      icolors[is.na(maps::map.where("world", xco, yco))] <- NA
    }

    if(length(res)) {
     rect(pp[, 1] - res[1] / 2, pp[, 2] - res[2] / 2, pp[, 1] + res[1] / 2, pp[, 2] + res[2] / 2,
       col = icolors, border = NA, lwd = 0)
    } else {
      points(sp::SpatialPoints(cbind(xco, yco)), col = icolors, pch = p.pch, cex = p.cex)
    }
    colors <- rep(NA, length = length(colors))
  }
  args$ylab <- args$xlab <- args$main <- ""
  args$type <- NULL
  args$axes <- NULL
  lwd.p <- if(!is.null(args$lwd)) rep(args$lwd, length.out = n) else NULL
  if(is.null(lwd.p))
    lwd.p <- rep(1, length.out = n)
  lty.p <- if(!is.null(args$lty)) rep(args$lty, length.out = n) else NULL
  if(is.null(lty.p))
    lty.p <- rep(1, length.out = n)
  border.p <- if(!is.null(args$border)) rep(args$border, length.out = n) else NULL
  if(is.null(border.p))
    border.p <- rep("black", length.out = n)
  density.p <- if(!is.null(args$density)) rep(args$density, length.out = n) else NULL
  angle.p <- if(!is.null(args$angle)) rep(args$angle, length.out = n) else NULL
  if(is.null(angle.p))
    angle.p <- rep(90, length.out = n)

  for(poly in unique(poly.names.orig)) {
    for(i in which(poly.names.orig == poly)) {
      args$x <- map[[i]][, 1L]
      args$y <- map[[i]][, 2L]
      args$border <- border.p[i]
      args$angle <- angle.p[i]
      args$lwd <- lwd.p[i]
      args$lty <- lty.p[i]
      if(!is.null(density.p))
        args$density <- density.p[i]
      if(!is.null(x)){ 
        if(!is.na(k <- pmatch(poly, x$id))) {
          args$col <- colors[k]
          args$density <- NULL
        } else {
          if(!missing) next
          args$col <- NULL
          if(is.null(args$density))
            args$density <- 20L
        }
      } else args$col <- colors[i]
      do.call(graphics::polygon, 
        delete.args(graphics::polygon, args, 
        c("lwd", "cex", "lty")))
      if(names && !values) {
        args$polygon <- map[[i]]
        args$poly.name <- poly.names.orig[i]
        args$counter <- i
        args$cex <- cex.names
        do.call(centroidtext, delete.args(centroidtext, args, "font"))
      }
      if(values && !names) {
        args$polygon <- map[[i]]
        args$poly.name <- as.character(round(x$x[k], digits = digits))
        args$counter <- k
        args$cex <- cex.values
        do.call(centroidtext, delete.args(centroidtext, args, "font"))
      }
    }
  }

  if(legend) {
    if(is.null(args$pos))
      args$pos <- "topleft"
    if(args$pos[1L] == "right") {
      args$full <- TRUE
      args$side.legend <- 2L
      args$side.ticks <- 2L
      mar <- mar.orig
      mar[2L] <- 0.5
      mar[4L] <- 3.1
      par(mar = mar, xaxs = "i", yaxs = "i")
      args$plot <- TRUE
      args$add <- FALSE
    } else {
      args$plot <- FALSE
      if(is.null(args$xpd))
        args$xpd <- TRUE
      args$add <- TRUE
    }
    args$shift <- args$legend.shift
    args$xlim <- map.limits$xlim
    args$ylim <- map.limits$ylim
    args$color <- col
    args$ncol <- ncol
    args$x <- x$x
    args$breaks <- breaks
    args$swap <- swap
    args$digits <- digits
    args$cex.labels <- cex.legend
    args$symmetric <- symmetric
    if(is.null(range)) {
      range <- range(args$x)
      if(diff(range) == 0)
        range <- unique(range) + c(-1, 1)
    }
    args$range <- range
    if(is.null(args$lrange))
      args$lrange <- args$range
    do.call(colorlegend, delete.args(colorlegend, args, c("font")))
  }
  if(!is.null(args$xlab))
    mtext(args$xlab, side = 1L)
  if(!is.null(args$ylab))
    mtext(args$ylab, side = 2L)

  return(invisible(NULL))
}


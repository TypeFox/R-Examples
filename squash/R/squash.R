# squash.R
#
# Aron Charles Eklund
#
# A part of the "squash" R package


######### color palettes ###########

rainbow2 <- function(n) rainbow(n, end = 0.8)

grayscale <- function(n, start = 0.9, end = 0) gray(seq(start, end, length = n))

greyscale <- grayscale

blueorange <- colorRampPalette(c('blue', 'lightgrey', '#FF7A00'), interpolate = 'spline')

bluered <- colorRampPalette(c('blue', 'lightgrey', 'red'), interpolate = 'spline')

darkbluered <- colorRampPalette(c(rgb(0.1, 0.1, 0.4), grey(0.9), rgb(0.4, 0.1, 0.1)), interpolate = 'spline')

jet <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

heat <- colorRampPalette(c('black', 'darkred', 'orange', 'yellow'))

coolheat <- colorRampPalette(c('cyan', '#00A5FF', '#00008B', 'black', 'black', 'darkred', 'orange', 'yellow'))


######### essential color mapping functions ###########

.extendrange2 <- function (x, r = range(x, na.rm = TRUE), f = 0.05) 
{
    if (!missing(r) && length(r) != 2) 
        stop("'r' must be a \"range\", hence of length 2")
    if(all(r == 0))
      return(c(-f, f))
    if(r[1] == r[2]) {
     return(r + c(-f, f) * r[1])
    } else {
     return(r + c(-f, f) * diff(r))
    }
}


makecmap <- function(x, n = 10, 
    breaks = pretty, symm = FALSE, base = NA, 
    colFn = jet, col.na = NA,
    right = FALSE, include.lowest = FALSE, ...) {
  if(missing(x) && (missing(breaks) || !is.numeric(breaks))) 
    stop("either 'x' or numeric 'breaks' must be supplied")
  if((missing(breaks) || !is.numeric(breaks))) {  ## Need to calculate breakpoints
    if (!is.numeric(x)) 
       stop("'x' must be numeric")
    if(!is.na(base)) {
      x <- log(x, base = base)
    }
    ## define "lim" to extend the range if necessary:
    lim <- range(x, finite = TRUE)
    if(!include.lowest) {
      if(right) 
        lim[1] <- .extendrange2(r = lim, f = 0.001)[1]
      else
        lim[2] <- .extendrange2(r = lim, f = 0.001)[2]
    }  
    if(symm) {
      maxabs <- max(abs(lim))
      lim <- c(-maxabs, maxabs)
    }
    breakpoints <- breaks(c(x, lim), n = n, ...)
    if(!is.na(base)) {
      breakpoints <- base ^ breakpoints
    }  
  } else {                         ## Breakpoints were supplied
    breakpoints <- breaks
  }
  if (length(breakpoints) < 2)
    stop("there must be at least two breakpoints")
  if (any(duplicated(breakpoints))) 
    stop("breakpoints are not unique")
  myColors <- colFn(length(breakpoints) - 1)
  list(breaks = breakpoints, colors = myColors, 
    base = base, col.na = col.na, right = right, 
    include.lowest = include.lowest)
}


cmap <- function(x, map, outlier = NULL, ...) {
  if(missing(map)) map <- makecmap(x, ...)
  out.i <- cut(x, breaks = map$breaks, right = map$right, 
      include.lowest = map$include.lowest)
  out <- map$colors[as.integer(out.i)]
  dim(out) <- dim(x)
  dimnames(out) <- dimnames(x)
  names(out) <- names(x)
  n.outliers <- sum(is.na(x) != is.na(out))
  if(n.outliers > 0) {
    if (is.null(outlier)) 
      stop('Found ', n.outliers, ' values outside map range.')
    else 
      warning(n.outliers, ' values outside map range will be colored ', outlier, '.')
    out[is.na(x) != is.na(out)] <- outlier
  }
  out[is.na(x)] <- map$col.na
  out
}

######### functions similar to "pretty" ###########


prettyInt <- function(x, n = 5, ...) {
  x <- x[is.finite(x <- as.numeric(x))]
  if (length(x) == 0L) return(x)
  r <- range(x)
  if (r[2L] == r[1L]) n <- 5  ## preempts certain problems
  r.adj <- c(floor(r[1L]), ceiling(r[2L]))
  p <- pretty(r.adj, n = n, ...)
  unique(as.integer(round(p)))
}



prettyLog <- function(x, n = 5, small = NA, 
     logrange = c(-100, 100)) {
  x <- x[is.finite(x <- as.numeric(x))]
  if (length(x) == 0L) return(x)
  r <- range(x)
  if(all(r == 0)) return(pretty(0))
  rmin <- range(abs(x[x != 0]), finite = TRUE)[1L]
  if(!is.na(small))
    rmin <- max(rmin, small) 
  expVals <- seq(logrange[1L], logrange[2L])
  ## narrow down candidate sequences one step at a time
  cand <- list(base2 = c(1, 2, 5) * 10 ^ rep(expVals, each = 3L),
               base3 = c(1, 3) * 10 ^ rep(expVals, each = 2L),
               base10 = 10 ^ rep(expVals))
  cand2 <- lapply(cand, function(y) y[(which(y > rmin)[1L] - 1L):length(y)])
  cand3 <- lapply(cand2, function(y) c(sort(-y), 0, y))
  cand4 <- lapply(cand3, function(y) y[(which(y > r[1L])[1L] - 1L):(which(y >= r[2L])[1L])])
  out.len <- sapply(cand4, length)
  wh <- which.min(abs(out.len - n))
  cand5 <- cand4[[wh]]
  if(length(cand5) == 1L) {
    if(cand5 > 0) return(c(0, cand5))
    if(cand5 < 0) return (c(cand5, 0))
    if(cand5 == 0) return (c(-1, 0))
  } else 
    return(cand5)
}


######### data manipulation functions ###########

matapply <- function(x, y = NULL, z = NULL, FUN, 
    nx = 50, ny = nx, 
    xlim = NULL, ylim = NULL, 
    xbreaks = NULL, ybreaks = NULL, 
    right = FALSE, include.lowest = TRUE, ...) {
  if(is.matrix(z) || (is.null(z) && is.matrix(x))) {
    xyz <- xyzmat2xyz(x, y, z)
  } else {
    xyz <- xyz.coords(x, y, z)
  }
  if(is.null(xlim)) 
    xlim <- range(xyz$x, finite = TRUE)
  if(is.null(ylim)) 
    ylim <- range(xyz$y, finite = TRUE)
  if(is.null(xbreaks))
    xbreaks <- pretty(xlim, nx)
  if(is.null(ybreaks))
    ybreaks <- pretty(ylim, ny)
  xcut <- cut(xyz$x, breaks = xbreaks, 
    right = right, include.lowest = include.lowest)
  ycut <- cut(xyz$y, breaks = ybreaks, 
    right = right, include.lowest = include.lowest)
  z.out <- tapply(xyz$z, list(xcut, ycut), FUN, ...)
  if(is.null(xyz$xlab)) xyz$xlab <- deparse(substitute(x))
  if(is.null(xyz$ylab)) xyz$ylab <- deparse(substitute(y))
  if(is.null(xyz$zlab)) xyz$zlab <- deparse(substitute(z))  
  list(x = xbreaks, y = ybreaks, z = z.out, 
    xlab = xyz$xlab, ylab = xyz$ylab, 
    zlab = paste(deparse(substitute(FUN)), '(', xyz$zlab, ')', sep = ''))
}



# returns custom "xyzmat.coords" from various types of input
xyzmat.coords <- function(x = NULL, y = NULL, z = NULL, 
    xlab = NULL, ylab = NULL, zlab = NULL, 
    xds = NULL, yds = NULL, zds = NULL) {
  if(is.null(z)) {
    if(is.matrix(x)) {
      z <- x
      if(is.null(zlab)) zlab <- xds
      x <- NULL
    } else if (is.list(x)) {
      if(is.null(xlab)) 
        if(!is.null(x$xlab)) xlab <- x$xlab 
        else xlab <- paste(xds, 'x', sep = '$')
      if(is.null(ylab)) 
        if(!is.null(x$ylab)) ylab <- x$ylab 
        else ylab <- paste(xds, 'y', sep = '$')
      if(is.null(zlab)) 
        if(!is.null(x$zlab)) zlab <- x$zlab 
        else zlab <- paste(xds, 'z', sep = '$')
      z <- x$z
      y <- x$y
      x <- x$x
    }
  } else if (is.null(y) && !is.null(x) && is.list(x)) {
    if(is.null(xlab)) 
      if(!is.null(x$xlab)) xlab <- x$xlab 
      else xlab <- paste(xds, 'x', sep = '$')
    if(is.null(ylab)) 
      if(!is.null(x$ylab)) ylab <- x$ylab 
      else ylab <- paste(xds, 'y', sep = '$')
    y <- x$y
    x <- x$x
  } else {  # z is specified, x and y maybe
    if(is.null(xlab) && !is.null(x)) xlab <- xds
    if(is.null(ylab) && !is.null(y)) ylab <- yds
    if(is.null(zlab)) zlab <- zds
  }
  if(!is.matrix(z)) stop ("a matrix must be specified")
  nr <- nrow(z)
  nc <- ncol(z)
  if(is.null(x)) {
    x <- 1:nr
    if(is.null(xlab)) xlab <- 'Row index'
  }
  if(is.null(y)) {
    y <- 1:nc
    if(is.null(ylab)) ylab <- 'Column index'
  }
  if(!length(x) %in% nr:(nr + 1)) 
    stop ("length of 'x' must be ", nr, " or ", nr + 1, " to match nrow(z)")
  if(!length(y) %in% nc:(nc + 1)) 
    stop ("length of 'y' must be ", nc, " or ", nc + 1, " to match ncol(z)")
  list(x = x, y = y, z = z, xlab = xlab, ylab = ylab, zlab = zlab)
}

# returns standard "xyz.coords" from a matrix input
xyzmat2xyz <- function(...) {
  m <- xyzmat.coords(...)
  if(length(m$x) != nrow(m$z)) {
    m$x <- 1:nrow(m$z)
    warning("Cannot get x coordinates; using indices instead.")
  }
  if(length(m$y) != ncol(m$z)) {
    m$y <- 1:ncol(m$z)
    warning("Cannot get y coordinates; using indices instead.")
  }
  xyz.coords(x = m$x[row(m$z)], y = m$y[col(m$z)], z = as.vector(m$z), 
    xlab = m$xlab, ylab = m$ylab, zlab = m$zlab)
}



######### plotting functions, sorted from low-level to high-level ###########
## (but not key-drawing functions, which are in the following section)


cimage <- function(x = NULL, y = NULL, zcol = NULL, zsize = 1, 
    xlab= NULL, ylab = NULL, xlabels = NULL, ylabels = NULL, 
    border = NA, add = FALSE, axes = TRUE, useRaster = FALSE, ...) {
  xyzmat <- xyzmat.coords(x = x, y = y, z = zcol, 
    xlab = xlab, ylab = ylab,
    xds = deparse(substitute(x)), yds = deparse(substitute(y)))
  x <- xyzmat$x
  y <- xyzmat$y
  zcol <- xyzmat$z
  xlab <- xyzmat$xlab
  ylab <- xyzmat$ylab
  if(is.null(xlab)) xlab <- NA
  if(is.null(ylab)) ylab <- NA
  nr <- nrow(zcol)
  nc <- ncol(zcol)
  if(min(zsize, na.rm = TRUE) < 0 || max(zsize, na.rm = TRUE) > 1)
    stop ("expected 0 <= 'zsize' <= 1") 
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
  	stop("increasing 'x' and 'y' values expected")
  if (length(x) == nr) {   # "x" supplies midpoints
    xmid <- x
    if(nr > 1) {
      dx <- 0.5 * diff(xmid)
      xbreaks <- c(xmid[1] - dx[1], xmid[-length(xmid)] + dx, 
                   xmid[length(xmid)] + dx[length(xmid) - 1])
    } else {
      xbreaks <- c(xmid - 0.5, xmid + 0.5)
    }
  } else {   # "x" supplies breakpoints
    xbreaks <- x
    xmid <- x[-1] - diff(x)/2
  }
  if (length(y) == nc) {   # "y" supplies midpoints
    ymid <- y
    if(nc > 1) {
      dy <- 0.5 * diff(ymid)
      ybreaks <- c(ymid[1] - dy[1], ymid[-length(ymid)] + dy, 
                   ymid[length(ymid)] + dy[length(ymid) - 1])
    } else {
      ybreaks <- c(ymid - 0.5, ymid + 0.5)
    }
  } else {   # "y" supplies breakpoints
    ybreaks <- y
    ymid <- y[-1] - diff(y)/2
  }
  xadj <- diff(xbreaks)[row(zcol)] * (1 - zsize) / 2
  yadj <- diff(ybreaks)[col(zcol)] * (1 - zsize) / 2
  xleft <- xbreaks[row(zcol)] + xadj
  xright <- xbreaks[row(zcol) + 1] - xadj
  ybottom <- ybreaks[col(zcol)] + yadj
  ytop <- ybreaks[col(zcol) + 1] - yadj
  if( !add ) {  # draw frame, axes, labels
    plot(range(xbreaks), range(ybreaks), 
      type = 'n', xaxs = 'i', yaxs = 'i', 
      xlab = xlab, ylab = ylab, axes = FALSE, ...)
    if(axes) {
      if(is.null(xlabels))
        axis(1, ...)    
      else if(length(xlabels) == 1 && is.logical(xlabels) && xlabels)
        axis(1, at = xmid, labels = rownames(xyzmat$z), ...)
      else  
        axis(1, at = xmid, labels = xlabels, ...)
      if(is.null(ylabels))
        axis(2, ...)    
      else if(length(ylabels) == 1 && is.logical(ylabels) && ylabels)
        axis(2, at = ymid, labels = colnames(xyzmat$z), ...)
      else  
        axis(2, at = ymid, labels = ylabels, ...)
      box(...)
    }
  }
  canUseRaster <- (length(zsize) == 1) && (zsize == 1) &&
      diff(range(diff(xbreaks))) / diff(range(xbreaks)) < 0.01 &&
      diff(range(diff(ybreaks))) / diff(range(ybreaks)) < 0.01 
  if(is.na(useRaster)) useRaster <- canUseRaster
  if(useRaster) {               # use rasterImage
    if(!canUseRaster) stop ("cannot use useRaster = TRUE with non-uniform breakpoints")
    flip <- function(x) t(x)[ncol(x):1, ] #rasterImage changes orientation
    rasterImage(flip(zcol), xleft = min(xleft), ybottom = min(ybottom), 
      xright = max(xright), ytop = max(ytop), 
      interpolate = FALSE)
  } else {                       # draw individual rectangles
    rect(xleft = xleft, ybottom = ybottom, 
       xright = xright, ytop = ytop, 
       col = zcol, border = border )
  }
}


colorgram <- function(x = NULL, y = NULL, z = NULL, zsize = 1, 
    map, nz = 10, breaks = pretty, symm = FALSE, base = NA, colFn = jet,
    key = hkey, key.args = list(), 
    xlab = NULL, ylab = NULL, zlab = NULL, 
    outlier = NULL, ...) {
  xyzmat <- xyzmat.coords(x = x, y = y, z = z, 
    xlab = xlab, ylab = ylab, zlab = zlab, 
    xds = deparse(substitute(x)), yds = deparse(substitute(y)), 
    zds = deparse(substitute(z)))
  if(missing(map)) {
    map <- makecmap(xyzmat$z, n = nz, breaks = breaks, 
      symm = symm, base = base, colFn = colFn)
  }
  zcol <- cmap(xyzmat$z, map = map, outlier = outlier)
  cimage(x = xyzmat$x, y = xyzmat$y, zcol = zcol, zsize = zsize, 
    xlab = xyzmat$xlab, ylab = xyzmat$ylab, ... )
  if(is.function(key) || is.character(key)) {
    if('title' %in% names(key.args))
      stop('the color key title cannot be specified in "key.args"; please set "zlab" instead')
    if('map' %in% names(key.args))
      stop('the color map cannot be specified in "key.args"; please set "map" directly')
    do.call(key, args = c(list(map = map, title = xyzmat$zlab), key.args))
  }
  invisible(map)
}  


squashgram <- function(x, y = NULL, z = NULL, FUN, 
    nx = 50, ny = nx, xlim = NULL, ylim = NULL,
    xbreaks = NULL, ybreaks = NULL,
    xlab = NULL, ylab = NULL, zlab = NULL, 
    shrink = 0, ...) {
  xyz <- xyz.coords(x = x, y = y, z = z, 
    xlab = xlab, ylab = ylab, zlab = zlab)
  if(is.null(xlab)) xlab <- xyz$xlab 
  if(is.null(xlab)) xlab <- deparse(substitute(x)) 
  if(is.null(ylab)) ylab <- xyz$ylab
  if(is.null(ylab)) ylab <- deparse(substitute(y)) 
  if(is.null(zlab)) zlab <- xyz$zlab
  if(is.null(zlab)) zlab <- deparse(substitute(z)) 
  zlab <- paste(deparse(substitute(FUN)), '(', zlab, ')', sep = '')
  sq <- matapply(x = xyz$x, y = xyz$y, z = xyz$z, 
    FUN = FUN, nx = nx, ny = ny, 
    xlim = xlim, ylim = ylim, xbreaks = xbreaks, ybreaks = ybreaks)
  if (shrink > 0) {
    h <- hist2(x = xyz$x, y = xyz$y, 
      xbreaks = sq$x, ybreaks = sq$y, plot = FALSE)
    zsize <- sqrt(pmin(h$z, shrink) / shrink)
  } else {
    zsize <- 1
  }
  colorgram(sq, zsize = zsize, 
    xlab = xlab, ylab = ylab, zlab = zlab, ...)
}


hist2 <- function(x, y = NULL, 
       nx = 50, ny = nx,  
       xlim = NULL, ylim = NULL, 
       xbreaks = NULL, ybreaks = NULL, 
       plot = TRUE, 
       xlab = NULL, ylab = NULL, zlab = 'Counts', 
       colFn = heat, breaks = prettyInt, ...) {
  xy <- xy.coords(x, y)
  firstNonNull <- function(...) unlist(list(...))[1]
  xlab <- firstNonNull(xlab, xy$xlab, deparse(substitute(x)), 'x')
  ylab <- firstNonNull(ylab, xy$ylab, deparse(substitute(y)), 'y')
  z <- rep(1, length(xy$x)) # dummy variable for tapply
  h <- matapply(xy$x, xy$y, z, FUN = length, 
    nx = nx, ny = ny, xlim = xlim, ylim = ylim,
    xbreaks = xbreaks, ybreaks = ybreaks)
  h$xlab <- xlab
  h$ylab <- ylab
  h$zlab <- zlab
  if(plot) {
    colorgram(h, colFn = colFn, breaks = breaks, ...)
  }
  invisible(h)  
}


dendromat <- function(x, mat, 
    labRow = rownames(mat), labCol = colnames(mat),
    height = NA, gap = 0, matlabside = 2, border = NA, 
    cex.lab = par('cex.axis'), ...) {
  stopifnot(matlabside %in% c(2, 4))
  if(is(x, 'hclust')) x <- as.dendrogram(x)
  if(is.null(labRow)) labRow <- labels(x)
  n <- attr(x, 'members')
  if(nrow(mat) != n) stop("'nrow(mat)' must equal the number of leaves in 'x'")
  ord <- order.dendrogram(x)
  mat <- as.matrix(mat)
  par(usr = c(0,1,0,1))  ## to ensure a consistent state after function finishes
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  on.exit(clip(0,1,0,1), add = TRUE)  ## removes clipping region
  mar <- par('mar')
  if(is.na(height)) {
    h.mat = par('csi') * ncol(mat)
    h.lab = max(strwidth(labRow, 'inches')) * cex.lab
    h.mar = 2 * par('csi') # set lower margin to 2 lines
    h.tot = par('din')[2]
    height = (h.mat + h.lab + h.mar) / h.tot
    mar[1] <- 2 + (h.lab / par('csi'))
  }
  stopifnot(height < 1 && height > 0)
  layout(matrix(1:2, ncol = 1), heights = c(1 - height, height))
  par(mar = c(gap, mar[2:4]))
  plot(x, leaflab = 'none', ...)
  usr <- par('usr')
  par(mar = c(mar[1:2], 0, mar[4]), las = 2)
  plot(0, 0, type = 'n', axes = FALSE, xlab = '', ylab = '')
  par(usr = c(usr[1:2], 0.5, ncol(mat) + 0.5))
  cimage(mat[ord, , drop = FALSE], border = border, add = TRUE) 
  axis(matlabside, at = 1:ncol(mat), labels = labCol, 
    tick = FALSE, line = -0.5)
  axis(1, at = 1:nrow(mat), labels = labRow[ord], 
    tick = FALSE, line = -0.5, cex.axis = cex.lab)
}


######### functions to draw color keys ###########


vkey <- function(map, title = NA, side = 2, stretch = 1.4, x, y, skip, wh) {
  if(!missing(skip) && !missing(wh)) stop ("cannot specify both 'skip' and 'wh'")
  opar <- par(xpd = NA)
  on.exit(par(opar))
  n <- length(map$breaks)
  dy <- strheight("A")
  aspect <- diff(grconvertX(1:2, from ='inches')) / diff(grconvertY(1:2, from ='inches'))
  dx <- dy * aspect
  if(missing(wh)) {
    if(missing(skip)) { # put as many labels as fit nicely
      for (i in 1:min(n, 20)) {
        if((n - 1) %% i == 0) {
          step <- (n - 1) / i
          wh.tmp <- seq(1, n, by = step)
          if(strheight("A") * 1.2 < dy * step * stretch) wh <- wh.tmp
        }
      }
    } else {
      wh <- seq(1, n, by = skip)
    }
  }
  labs <- format(map$breaks[wh])
  maxlabwidth <- max(strwidth(labs))
  if(missing(x)) {
    x <- grconvertX(1, from = 'nfc') - (2 * dx)
    if(side == 4) x <- x - maxlabwidth - dx
  } else {
    if(is.list(x)) {
      y <- x$y
      x <- x$x
    }
  }
  if(missing(y)) y <- par('usr')[3] + dy
  ybord <- y + ((0:(n-1)) * dy * stretch)
  rect(x, ybord[-n], x + dx, ybord[-1], col = map$colors, border = NA)
  if(side == 4) {
    xtext <- x + dx
    text(x = x, y = ybord[n] + (1.5 * dy), title, adj = c(0, 0))
  }
  if(side == 2) {
    xtext <- x
    text(x = x + dx, y = ybord[n] + (1.5 * dy), title, adj = c(1, 0))
  }
  text(x = xtext, y = ybord[wh], 
    labels = labs, pos = side)
}


hkey <- function(map, title = NA, side = 1, stretch = 1.4, x, y, skip, wh) {
  if(!missing(skip) && !missing(wh)) stop ("cannot specify both 'skip' and 'wh'")
  opar <- par(xpd = NA)
  on.exit(par(opar))
  n <- length(map$breaks)
  dy <- strheight("A")
  aspect <- diff(grconvertX(1:2, from ='inches')) / diff(grconvertY(1:2, from ='inches'))
  dx <- dy * aspect
  labs <- format(map$breaks)
  labwidth <- strwidth(labs)
  if(missing(x)) {
    x <- grconvertX(0, from = 'nfc') + dx + (0.5 * strwidth(format(map$breaks[1])))
  } else {
    if(is.list(x)) {
      y <- x$y
      x <- x$x
    }
  }
  if(missing(y)) y <- grconvertY(0, from = 'nfc') + (2 * dy)
  xbord <- x + ((0:(n-1)) * dx * stretch)
  if(missing(wh)) {
    if(missing(skip)) { # put as many labels as fit nicely
      for (i in 1:min(n, 20)) {
        if((n - 1) %% i == 0) {
          step <- (n - 1) / i
          wh.tmp <- seq(1, n, by = step)
          maxlabwidth <- max(strwidth(format(map$breaks[wh.tmp])))
          if(maxlabwidth + dx < dx * step * stretch) wh <- wh.tmp
        }
      }
    } else {
      wh <- seq(1, n, by = skip)
    }
  }
  rect(xbord[-n], y, xbord[-1], y + dy, col = map$colors, border = NA)
  if(side == 1) {
    ytext <- y
    text(x = x, y = y + (1.5 * dy), title, adj = c(0, 0))
  }
  if(side == 3) {
    ytext <- y + dy
    text(x = x, y = y - (0.5 * dy), title, adj = c(0, 1))
  }
  text(x = xbord[wh], y = ytext, 
    labels = format(map$breaks[wh]), pos = side)
}


######### file-writing function ###########


savemat <- function(x, filename, map = NULL, outlier = NULL, 
    dev = c('png', 'pdf', 'bmp', 'tiff', 'jpeg'), do.dev.off = TRUE, ...) {
  dev <- match.arg(dev)
  if(is.list(x)) x <- x$z
  if(dev != 'pdf') {
    do.call(dev, args = list(filename = filename, 
      width = nrow(x), height = ncol(x),
      antialias = 'none', ...))
  } else if(dev == 'pdf') {
    do.call(dev, args = list(file = filename, ...))
  }
  par(mar = c(0, 0, 0, 0), ann = FALSE, 
    xaxt = 'n', yaxt = 'n', bty = 'n')
  if(!is.null(map)) {
    x <- cmap(x, map = map, outlier = outlier)
  }
  cimage(x, useRaster = TRUE)
  if(do.dev.off)
    dev.off()
}



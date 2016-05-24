## TODO: I still do not deal with the case where there is only a
## single unique y value very nicely.  I should probably draw a delta
## function?
profiles.plot <- function(y, col.line, col.fill,
                          xlim=NULL, ymax=NULL,
                          n.br=50, opacity=.5,
                          xlab="Parameter estimate",
                          ylab="Probability density",
                          legend.pos=NULL,
                          with.bar=TRUE,
                          col.bg=NA, lwd=1, lines.on.top=TRUE,
                          ...) {
  if ( !is.data.frame(y) && is.matrix(y) )
    y <- as.data.frame(y)
  if ( missing(col.fill) )
    col.fill <- add.alpha(col.line, opacity)
  if ( missing(xlim) )
    xlim <- range(unlist(y))
  dx <- diff(xlim) / (n.br - 1)

  ## OK, so I need to do the histogram breaking better.  n.br is the
  ## number of breaks that would fit into xlim.  There are guaranteed
  ## to be breaks at seq(xlim[1], xlim[2], length.out=n.br).  However,
  ## it should never span more than the range of the observed data.
  f <- function(yi) {
    n <- (range(yi) - xlim)/dx
    n <- c(floor(n[1]), ceiling(n[2]))
    ri <- xlim + n * dx
    hist(yi, seq(ri[1], ri[2], by=dx), plot=FALSE)
  }

  hh <- lapply(y, f)
  ci <- lapply(y, hdr)

  if ( is.null(ymax) )
    ymax <- max(sapply(hh, function(x) max(x$density)))
  ## Change this once bottom arrow is optional...
  if ( with.bar )
    ylim <- c(-0.075, 1.05) * ymax
  else
    ylim <- c(0, 1.05 * ymax) # increase a little?
  
  plot(NA, xlim=xlim, ylim=ylim, type="n", yaxs="i", xlab=xlab,
       ylab=ylab, ...)
  if (!is.na(col.bg))
    for (i in seq_along(y))
      add.profile.shading(hh[[i]], ci[[i]], col.bg)
  if (lines.on.top) {
    for (i in seq_along(y))
      add.profile.shading(hh[[i]], ci[[i]], col.fill[i])
    for (i in seq_along(y))
      add.profile.outline(hh[[i]], col.line[i], lwd=lwd)
  } else {
    for (i in seq_along(y)) {
      add.profile.shading(hh[[i]], ci[[i]], col.fill[i])
      add.profile.outline(hh[[i]], col.line[i], lwd=lwd)
    }
  }
  if ( with.bar ) {
    z <- seq(0, 1, length.out=length(y) + 2)[-1] * par("usr")[3]
    for ( i in seq_along(y) )
      arrows(ci[[i]][1], z[i], ci[[i]][2], z[i], code=3, angle=90,
             length=0.02, col=col.line[i])
  }

  if ( !is.null(legend.pos) )
    legend(legend.pos, names(y), fill=col.fill, border=col.line, bty="n")
}

add.profile.shading <- function(h, ci, col) {
  dx <- diff(h$mids[1:2])
  xx <- c(h$mids - dx / 2, h$mids[length(h$mids)] + dx / 2)
  i <- which(xx > ci[1] & xx < ci[2])
  xs <- rep(c(ci[1], xx[i], ci[2]), each=2)
  j <- if ( length(i) > 1 ) max(1, min(i) - 1) else 1
  ys <- c(0, rep(h$density[c(j, i)], each=2), 0)
  graphics::polygon(xs, ys, col=col, border=NA)
}
add.profile.outline <- function(h, col, vertical=FALSE, ...) {
  dx <- diff(h$mids[1:2])
  if (vertical)
    graphics::lines(h, freq=FALSE, col=col)
  else {
    xx <- rep(with(h, c(mids-dx/2, mids[length(mids)]+dx/2)), each=2)
    yy <- c(0, rep(h$density, each=2), 0)
    graphics::lines(xx, yy, col=col, ...)
  }
}

hdr.uniroot <- function(z, p=0.95, lim=NULL) {
  xx <- sort(c(lim, seq(min(z), max(z), length.out=1024)))
  ez <- ecdf(z)
  f <- suppressWarnings(approxfun(ez(xx), xx))
  fit <- suppressWarnings(optimize(function(x)
                                   f(x + p) - f(x), c(0, 1-p)))
  if ( inherits(fit, "try-error") || is.na(fit$objective) )
    stop("HDR interval failure")
  ci <- fit$min
  f(c(ci, ci+p))
}

add.alpha <- function(col, alpha=.5) {
  if ( length(alpha) > 1 && any(is.na(alpha)) ) {
    n <- max(length(col), length(alpha))
    alpha <- rep(alpha, length.out=n)
    col <- rep(col, length.out=n)
    ok <- !is.na(alpha)
    ret <- rep(NA, length(col))
    ret[ok] <- add.alpha(col[ok], alpha[ok])
    ret
  } else {
    tmp <- grDevices::col2rgb(col)/255
    grDevices::rgb(tmp[1,], tmp[2,], tmp[3,], alpha=alpha)
  }
}

hdr <- function(z, p=0.95, lim=NULL) {
  ## ci <- try(hdr.uniroot(z, p, lim), silent=TRUE)
  ci <- try(hdr.new(z, 1-p, lim), silent=TRUE)
  if ( inherits(ci, "try-error") ) {
    warning("HDR falling back on quantile-based intervals")
    ci <- as.numeric(quantile(z, c((1-p)/2, 1/2 + p/2)))
  }
  ci
}

hdr.new <- function(z, alpha=0.05, lim=NULL) {
  spline.prep <- function(x, y) {
    o <- sort.list(x, method="quick", na.last=NA)
    x <- x[o]
    y <- y[o]
    nx <- length(x)
    if (length(ux <- unique(x)) < nx) {
      y <- as.vector(tapply(y, match(x, x), mean))
      x <- ux
    }
    list(x=x, y=y)
  }

  ez <- ecdf(z)
  sz <- environment(ez)$x
  tmp <- spline.prep(ez(sz), sz)
  .Call("hdr", tmp$x, tmp$y, alpha)
}

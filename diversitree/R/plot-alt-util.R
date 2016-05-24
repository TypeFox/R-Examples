## Things that are not really related to plotting trees, but might be
## more generally useful.

######################################################################
## Utility plotting functions
## TODO: These could probably do with checking, particularlyu around
## repetition of arguments, and the number of arguments that should be
## passed through.  In addition, getting the arc() function to return
## coordinates only when given plot=FALSE would simplify things like
## sectors(), as I would not have to redo the np and sequence
## calculations.
## np is the number of points expected around a full circle at the
## largest radius.
arcs <- function(theta0, theta1, r, col=par("fg"), lty=1,
                 lwd=par("lwd"), np=1000) {
  if ( length(r) != length(theta0) || length(r) != length(theta1) )
    stop("theta0, theta1 and r must be of same length")
  if ( any(lty[!is.na(lty)] != 1) )
    warning("lwd != 1 will probably be ugly for arcs")

  dx <- max(r) * 2 * pi / np
  nn <- pmax(2, ceiling((theta1 - theta0) * r / 2 / dx))

  tmp <- lapply(seq_along(nn), function(i)
                cbind(seq(theta0[i], theta1[i], length.out=nn[i]),
                      rep(r[i], nn[i])))
  tmp0 <- do.call(rbind, lapply(seq_along(nn), function(i)
                                rbind(tmp[[i]][-nn[i],])))
  tmp1 <- do.call(rbind, lapply(seq_along(nn), function(i)
                                rbind(tmp[[i]][-1,])))

  if ( length(lwd) > 1 )
    lwd <- rep(rep(lwd, length.out=length(theta0)), nn - 1)
  if ( length(col) > 1 )
    col <- rep(rep(col, length.out=length(theta0)), nn - 1)

  segments(tmp0[,2] * cos(tmp0[,1]), tmp0[,2] * sin(tmp0[,1]),
           tmp1[,2] * cos(tmp1[,1]), tmp1[,2] * sin(tmp1[,1]),
           lwd=lwd, col=col)
}

sectors <- function(theta0, theta1, r0, r1, ..., np=1000) {
  n <- length(theta0)
  if ( length(theta1) != n || length(r0) != n || length(r1) != n )
    stop("theta0, theta1, r0 and r1 must be same length")

  tm <- (theta0 + theta1) / 2
  dt <- abs(theta1 - theta0)

  dx <- max(r1) * 2 * pi / np
  nn <- pmax(2, ceiling(dt * r1 / 2 / dx))

  tt <- lapply(seq_along(nn), function(i)
               c(tm[i],
                 seq(theta0[i], theta1[i], length.out=nn[i]),
                 tm[i], NA))
  rr <- lapply(seq_along(nn), function(i)
               c(r0[i], rep(r1[i], nn[i]), r0[i], NA))

  xx <- unlist(rr) * cos(unlist(tt))
  yy <- unlist(rr) * sin(unlist(tt))

  graphics::polygon(xx, yy, ...)
}

radial.text <- function(r, theta, labels, cex=1, col="black",
                        font=1, ...) {
  n <- length(labels)
  col <- rep(col, length.out=n)
  x <- r * cos(theta)
  y <- r * sin(theta)
  srt <- theta / (2 * pi) * 360
  adj <- rep(0, n)
  i <- srt > 90 & srt < 270
  adj[i] <- 1
  srt[i] <- (srt[i] + 180) %% 360
  ## TODO: might be worth changing?
  if ( length(font) > 1 ) font <- font[1]
  if ( length(cex) > 1 ) cex <- cex[1]

  for ( i in seq_len(n) )
    text(x[i], y[i], labels[i], cex=cex, col=col[i],
         font=font, srt=srt[i], adj=adj[i], ...)
}

filled.arcs <- function(theta0, theta1, r, w, col=par("fg"), np=1000) {
  if ( length(r) == 1 )
    r <- rep(r, length(theta0))
  if ( length(r) != length(theta0) || length(r) != length(theta1) )
    stop("theta0, theta1 and r must be of same length")

  dx <- max(r) * 2 * pi / np
  nn <- pmax(2, ceiling((theta1 - theta0) * r / dx))

  f <- function(i) {
    t <- seq(theta0[i], theta1[i], length.out=nn[i])
    cbind(c(r[i] * cos(t), (r[i] + w) * cos(rev(t)), NA),
          c(r[i] * sin(t), (r[i] + w) * sin(rev(t)), NA))
  }
  tmp <- do.call(rbind, lapply(which(!is.na(col)), f))
  graphics::polygon(tmp, col=col[!is.na(col)], border=NA)
}

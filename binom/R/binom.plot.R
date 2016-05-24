panel.binom.plot.levelplot <- function(x, y, z, subscripts, breaks = NULL, ...) {
  panel.levelplot(x, y, z, subscripts, ...)
  if(!is.null(breaks)) {
    breaks <- breaks[subscripts, ]
    for(i in seq(nrow(breaks))) {
      x.i <- with(breaks, c(lower[i], upper[i]))
      y.i <- breaks$x[i]
      p.x <- x.i[c(1, 2, 2, 1)]
      p.y <- y.i + c(-0.5, -0.5, 0.5, 0.5)
      lpolygon(p.x, p.y, border = "#cccccc", lwd = 3)
    }
  }
}

panel.binom.plot.xyplot <-  function(x, y, subscripts, conf.level, n, breaks, actual, ...) {
  panel.abline(h = actual, lty = 2, lwd = 2, col = "#880000")
  n <- unique(n[subscripts])
  breaks <- unique(sort(unlist(breaks[breaks$n == n, 2:3])))
  nb <- length(breaks)
  if(any(m <- breaks %in% x))
    breaks[m] <- ifelse(breaks[m] > 0.5, breaks[m] - 1e-8, breaks[m] + 1e-8)
  x <- c(x, breaks)
  y <- c(y, rep(NA, nb))
  x <- x[ord <- order(x)]
  y <- y[ord]
  panel.xyplot(x, y, type = "l", ...)
  xx <- rep(breaks, each = 3)
  xx[seq(3, nb * 3, 3)] <- NA
  na <- which(is.na(y))
  wh.y <- rep(na, each = 3) + rep(c(-1, 1, 0), times = length(na))
  ny <- length(y)
  zero <- wh.y == 0
  ny.plus.1 <- wh.y == ny + 1
  end <- zero | ny.plus.1
  yy <- y[wh.y[!end]]
  if(any(end)) {
    if(any(zero)) yy <- c(NA, yy)
    if(any(ny.plus.1)) yy <- c(yy, NA)
  }
  yy[wh.y %in% na] <- NA
  panel.xyplot(xx, yy, type = "l", lty = 4, lwd = 2, col = "#888888")
}

binom.plot <- function(n,
                       method = binom.lrt,
                       np = 500,
                       conf.level = 0.95,
                       actual = conf.level,
                       type = c("xyplot", "levelplot"),
                       tol = .Machine$double.eps^0.5, ...) {
  stopifnot(require(lattice))
  type <- match.arg(type)
  if(length(n) != 1) {
    if(length(n) > 1 && type == "levelplot") {
      warning(sprintf("n must be of length 1, not %d", length(n)))
      n <- n[1]
    }
  }
  E.pn <- function(x, n, p, lower, upper) (p >= lower & p <= upper) * dbinom(x, n, p)
  p <- seq(tol, 1 - tol, length = np)
  args <- list(...)
  if(type == "levelplot") {
    x <- 0:n
    ci <- method(x, n, conf.level, ...)[c("x", "n", "lower", "upper")]
    z <- merge(ci, data.frame(p = p))
    z$coverage <- with(z, E.pn(x, n, p, lower, upper))
    z$n <- factor(sprintf("n = %d", n))
    args$x <- coverage ~ p * x | n
    if(is.null(args$col.regions))
      args$col.regions <- heat.colors(100)[100:1]
    if(is.null(args$panel))
      args$panel <- panel.binom.plot.levelplot
    args$breaks <- ci
    if(is.null(args$scales)) args$scales <- list(y = list(at = x, labels = x))
  } else {
    x <- unlist(lapply(lapply(n, ":", 0), rev))
    n <- rep(n, n + 1)
    ci <- method(x, n, conf.level, ...)[c("x", "n", "lower", "upper")]
    ci$lower <- ifelse(ci$lower < 0, 0, ci$lower)
    ci$upper <- ifelse(ci$upper > 1, 1, ci$upper)
    z <- merge(ci, data.frame(p = p))
    z$coverage <- with(z, E.pn(x, n, p, lower, upper))
    z <- aggregate(z["coverage"], z[c("p", "n")], sum)
    args$x <- coverage ~ p | n
    args$n <- z$n
    z$n <- factor(z$n, labels = sprintf("n = %d", sort(unique(z$n))))
    args$breaks <- ci[c("n", "lower", "upper")]
    if(is.null(args$panel))
      args$panel <- panel.binom.plot.xyplot
    if(is.null(args$ylab)) {
      args$ylab <- expression(E(paste(p,"|",n)))
    }
    args$conf.level <- conf.level
    args$actual <- actual
  }
  args$data <- z
  if(is.null(args$as.table))
    args$as.table <- TRUE
  do.call(type, args)
}

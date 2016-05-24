
### MLPLOT ###

.mlplot3 <- function (x, ..., main=NULL, add=FALSE) 
{
  x.name <- deparse(substitute(x))
  dx <- dim(x)
  n <- dx[3]
  n.col <- ceiling(sqrt(n))
  n.row <- n %/% n.col
  if (n.col*n.row<n) n.row <- n.row+1
  grid <- c(n.row, n.col)
  #
  d.n <- dimnames(x)[[3]]
  if (length(d.n)!=n) {
    d.n[n] <- NA
    def.names <- paste(x.name, "[,,", 1:n, "]", sep="")
    d.n[is.na(d.n)] <- def.names[is.na(d.n)]
  }
  par(mfrow=grid, oma=c(0,0,5,0), mar=c(2,1,1,0.5)+0.1)
  for (k in 1:n) {
    main.title <- d.n[k]
    mlplot(x[,,k], ..., main=main.title, add=FALSE)
  }
  if (is.null(main)) {
    main <- x.name
  }
  mtext(main, line=1, outer=TRUE)
  invisible()
}

mlplot <- function (X, ...)
{
  UseMethod("mlplot")
}

mlplot.default <- function (X, y.center = TRUE, y.shift = 0, y.map = NULL, mar = par("mar"), left.margin = 3, vline=NULL, top.axis = TRUE, exp.labels=FALSE, x.ticks = NULL, axes = NULL, xlim = NULL, ylim = NULL, xlab=deparse(substitute(X)), ylab=NULL, las = NULL, add = FALSE, ...) 
{
    if (missing(xlab)) {
      xlab <- deparse(substitute(X)) #?
    }
    dx <- dim(X)
    lx <- length(dx)
    if (lx > 3) {
        stop("Too many dimensions: ", lx)
    }
    else if (lx == 3) {
        return(.mlplot3(X, y.center=y.center, y.shift=y.shift, y.map = y.map,
             mar = mar, left.margin = left.margin, top.axis = top.axis, 
             exp.labels=exp.labels, x.ticks=x.ticks, axes = axes, 
             xlim=xlim, ylim = ylim,
             xlab=xlab, ylab = ylab,
             las = las,  
             add = add, ...))
    }
    else if (lx == 2) {
        labels <- dimnames(X)[[1]]
        if (any(rv.all.na(X))) {
            f <- function(x) {
                if (!any(is.na <- rv.all.na(x))) 
                  return(x)
                w <- which(is.na)
                c(x[-w], x[w])
            }
            X <- t(apply.rv(X, 1, f))
        }
    }
    else {
        labels <- names(X)
        dim(X) <- c(length(X), 1)
    }
    y.row.coords <- 1:nrow(X)
    ylim <- rev(range(y.row.coords) + c(-1, 1))
    if (is.null(y.map)) {
        y.map <- function(x) {
            if (y.center) {
                f <- function(i, j, n, nc) {
                  i + (j - 1)/n - 0.5 * (nc - 1)/n
                }
            }
            else {
                f <- function(i, j, n, nc) {
                  i + (j - 1)/n
                }
            }
            f(i = row(x), j = col(x), n = max(ncol(x),10)*1.5, nc = ncol(x))
        }
    }
    if (is.function(y.map)) {
        y <- y.map(X)
    }
    else if (is.numeric(y.map)) {
        y <- y.map
    }
    else {
        stop("Unknown type of y.map")
    }
    if (is.numeric(y.shift)) {
        y <- (y + y.shift)
    }
    else {
        stop("y.shift must be numeric")
    }
    if (length(y) != length(X)) {
        stop("y coordinates are not valid (check y.map!)")
    }
    if (is.null(xlim)) {
        x.sims <- sims(as.rvobj(X))
        rng <- range(x.sims[is.finite(x.sims)])
        if (length(rng) < 2) {
            rng <- c(-1, 1)
        }
        xlim <- rng
    }
    if (is.null(x.ticks)) {
      x.row.coords <- pretty(xlim)
    } else {
      x.row.coords <- x.ticks
    }
    if (! is.null(names(x.ticks))) {
      x.labels <- names(x.ticks)
    } else if (exp.labels) {
      x.labels <- paste(signif(exp(x.row.coords),2))
    } else {
      x.labels <- paste(x.row.coords)
    }
    mar <- (mar + c(0, left.margin, 2, 0))
    oldpar <- par(mar = mar)
    on.exit(par(oldpar))
    las <- if (!is.null(las)) las else 1
    if (add) {
      points(X, y, xlim = xlim, ylim = ylim, ...)
    } else {
      plot(X, y, ..., las = las, xlim = xlim, ylim = ylim, 
           axes = FALSE, xlab=xlab, ylab = "", type="n")
      if (is.null(axes) || axes) {
        axis(1, at = x.row.coords, labels=x.labels)
        if (top.axis) 
          axis(3, at = x.row.coords, labels=x.labels)
        if (is.null(labels)) 
          labels <- paste(y.row.coords)
        axis(2, at = y.row.coords, labels = labels, tick = FALSE, 
             line = FALSE, pos = NA, outer = FALSE, font = NA, 
             las = 1)
        if (! is.null(vline)) {
          if (! is.numeric(vline)) {
            stop("'vline' must be a numeric vector (or NULL if not used)")
          }
          if (is.null(names(vline))) {
            names(vline) <- rep("dotted", length(vline))
          }
          for (i in seq_along(vline)) {
            lty <- names(vline)[i]
            abline(v=vline[i], lty=lty, col="gray")
          }
        }
      }
      points(X, y, xlim = xlim, ylim = ylim, ...)
    }
    invisible(NULL)
}


mlplot_OLD_rvsummary <- function (X, y.center = TRUE, y.shift = 0, y.map = NULL, mar = par("mar"), left.margin = 3, top.axis = TRUE, exp.labels=FALSE, x.ticks = NULL, axes = NULL, xlim = NULL, ylim = NULL, xlab=deparse(substitute(X)), ylab=NULL, las = NULL, add = FALSE, ...) # NOEXPORT
{
    if (missing(xlab)) {
      xlab <- deparse(substitute(X)) #?
    }
    dx <- dim(X)
    lx <- length(dx)
    if (lx > 3) {
        stop("Too many dimensions: ", lx)
    }
    else if (lx == 3) {
        return(.mlplot3(X, y.center=y.center, y.shift=y.shift, y.map = y.map,
             mar = mar, left.margin = left.margin, top.axis = top.axis, 
             exp.labels=exp.labels, x.ticks=x.ticks, axes = axes, 
             xlim=xlim, ylim = ylim,
             xlab=xlab, ylab = ylab,
             las = las,  
             add = add, ...))
    }
    else if (lx == 2) {
        labels <- dimnames(X)[[1]]
        if (any(rv.all.na(X))) {
            f <- function(x) {
                if (!any(is.na <- rv.all.na(x))) 
                  return(x)
                w <- which(is.na)
                c(x[-w], x[w])
            }
            X <- t(apply.rv(X, 1, f))
        }
    }
    else {
        labels <- names(X)
        dim(X) <- c(length(X), 1)
    }
    y.row.coords <- 1:nrow(X)
    ylim <- rev(range(y.row.coords) + c(-1, 1))
    if (is.null(y.map)) {
        myrow <- function (x) row(array(NA, dim(x)))
        mycol <- function (x) col(array(NA, dim(x)))
        y.map <- function(x) {
            if (y.center) {
                f <- function(i, j, n, nc) {
                  i + (j - 1)/n - 0.5 * (nc - 1)/n
                }
            }
            else {
                f <- function(i, j, n, nc) {
                  i + (j - 1)/n
                }
            }
            f(i = myrow(x), j = mycol(x), n = max(ncol(x),10)*1.5, nc = ncol(x))
        }
    }
    if (is.function(y.map)) {
        y <- y.map(X)
    }
    else if (is.numeric(y.map)) {
        y <- y.map
    }
    else {
        stop("Unknown type of y.map")
    }
    if (is.numeric(y.shift)) {
        y <- (y + y.shift)
    }
    else {
        stop("y.shift must be numeric")
    }
    if (length(y) != length(X)) {
        stop("y coordinates are not valid (check y.map!)")
    }
    if (is.null(xlim)) {
        x.sims <- sims(X)
        rng <- range(x.sims[is.finite(x.sims)])
        if (length(rng) < 2) {
            rng <- c(-1, 1)
        }
        xlim <- rng
    }
    if (is.null(x.ticks)) {
      x.row.coords <- pretty(xlim)
    } else {
      x.row.coords <- x.ticks
    }
    if (exp.labels) {
      x.labels <- paste(signif(exp(x.row.coords),2))
    } else {
      x.labels <- paste(x.row.coords)
    }
    mar <- (mar + c(0, left.margin, 2, 0))
    oldpar <- par(mar = mar)
    on.exit(par(oldpar))
    las <- if (!is.null(las)) las else 1
    if (add) {
        points(X, y, xlim = xlim, ylim = ylim, ...)
    }
    else {
        plot(X, y, ..., las = las, xlim = xlim, ylim = ylim, 
            axes = FALSE, xlab=xlab, ylab = "")
        if (is.null(axes) || axes) {
            axis(1, at = x.row.coords, labels=x.labels)
            if (top.axis) 
                axis(3, at = x.row.coords, labels=x.labels)
            if (is.null(labels)) 
                labels <- paste(y.row.coords)
            axis(2, at = y.row.coords, labels = labels, tick = FALSE, 
                line = FALSE, pos = NA, outer = FALSE, font = NA, 
                las = 1)
        }
    }
    invisible(NULL)
}



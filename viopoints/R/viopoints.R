viopoints <- function(x, ...) UseMethod("viopoints")

viopoints.formula <- function(formula, data=NULL, ..., subset, na.action=NULL)
{
  if (missing(formula) || (length(formula) != 3L))
    stop("'formula' missing or incorrect")
  m <- match.call(expand.dots=FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m$... <- NULL
  m$na.action <- na.action
  require(stats, quietly=TRUE)
  m[[1L]] <- as.name("model.frame")
  mf <- eval(m, parent.frame())
  response <- attr(attr(mf, "terms"),  "response")

  groups <- mf[-response]
  if (ncol(groups) == 0L) {
    if (is.matrix(mf[[response]]))
      viopoints(mf[[response]], ...)
    else
      viopoints(mf[response], ...)
  } else {
    if (2L <= ncol(groups))
      groups <- interaction(groups)
    if (is.matrix(mf[[response]]))
      viopoints(mf[[response]], groups=groups, ...)
    else
      viopoints(mf[response], groups=groups, ...)
  }
}

viopoints.default <- function(x, ..., groups, na.group=FALSE,
  method="violin", side="both", amount=jitter, jitter=0.2, offset=1, 
  density.par=list(na.rm=TRUE), horizontal=!vertical, vertical=TRUE, at,
  points=TRUE, pch=par("pch"), cex=par("cex"), col="red", bg="pink",
  lines=FALSE, line.lty=par("lty"), line.lwd=0.5, line.col="lightgray", 
  add=FALSE, axes=TRUE, frame.plot=axes, axis.par, group.names, 
  main="", sub="", xlab, ylab, dlab="", glab="", xlim, ylim, log="") 
{
  method <- pmatch(method, c("overplot", "stack", "jitter", "violin"))[1L]
  if (is.na(method) || method == 0L)
    stop("invalid plotting method")
  side <- pmatch(side, c("negative", "both", "positive"), duplicates.ok=TRUE)
  if (is.na(side) || side == 0L)
    stop("invalid plotting side")

  unbind <- function(x, margin=2) {
    if (is.list(x)) return(lapply(x, unlist))
    if (is.matrix(x)) return(unlist(apply(x, margin, list), recursive=FALSE))
    if (is.vector(x)) return(list(x))
  }
  rebind <- function(L) unlist(lapply(L, unbind), recursive=FALSE)
  seqs <- function(x, side) {
    l <- length(x)
    if (side == 1L) return(0.5 - seq_len(l))
    if (side == 2L) return(seq_len(l) - (l + 1) / 2) 
    if (side == 3L) return(seq_len(l) - 0.5)
  }

  x <- rebind(list(x,...))
  n <- length(x)
  if (missing(groups)) {
    m <- 1L
    mn <- n
    reps <- function(x) rep(x, length.out=n)
  } else {
    if (!is.factor(groups)) {
      if (is.list(groups))
        groups <- unlist(groups)
      groups <- factor(groups, exclude=NULL)
    }
    f <- as.integer(groups)
    m <- nlevels(groups)
    if (is.na(levels(groups)[m]) && !na.group) 
      m <- m - 1
    mn <- m * n
    reps <- function(x) {
      if (length(x) <= m) 
        rep(rep(x, length.out=m), each=n) 
      else
        rep(x, length.out=mn)
    }
  }
  if (points)
    pp <- lapply(list(pch=pch, cex=cex, col=col, bg=bg), reps)
  if (lines) 
    lp <- lapply(list(lty=line.lty, lwd=line.lwd, col=line.col), reps)
  side <- rep(side, length.out=mn)
  at <- if (missing(at)) seq_len(mn) else rep(at, length.out=mn)

  if (!add) {
    uat <- unique(at)
    if (missing(group.names)) {
      if ((n == 1L) && (1L < length(uat)) && !missing(groups))
        group.names <-  levels(groups)
      else
        group.names <- if (is.null(names(x))) uat else names(x)
    }
    group.names <- rep(group.names, length=length(uat))
    group.names[is.na(group.names)] <- "<NA>"
    plot.new()
    if (horizontal) {
      if (missing(xlim)) 
        xlim <- range(unlist(x, use.names=FALSE), finite=TRUE)
      if (missing(ylim)) 
        ylim <- range(uat) + 
          ifelse(length(uat) < 2L, 1, min(diff(sort(uat)))) * c(-0.5,0.5)
      plot.window(xlim=xlim, ylim=ylim, log=log)
      if (axes) {
        axis(1)
        if (missing(axis.par))
          axis.par <- list(at=uat, labels=group.names)
        do.call("axis", c(list(2), axis.par))
      }
    } else {
      if (missing(xlim))    
        xlim <- range(uat) + 
          ifelse(length(uat) < 2L, 1, min(diff(sort(uat)))) * c(-0.5,0.5)
      if (missing(ylim))
        ylim <- range(unlist(x, use.names=FALSE), finite=TRUE)
      plot.window(xlim=xlim, ylim=ylim, log=log)
      if (axes) {
        if (missing(axis.par))
          axis.par <- list(at=uat, labels=group.names)
        do.call("axis", c(list(1), axis.par))
        axis(2)
      }
    }
  }

  if (lines) {
    for (i in seq_len(m)) {
      k <- (i - 1) * n
      for (j in seq_len(n - 1)) {
        k <- k + 1
        z1 <- x[[j]]
        z2 <- x[[j+1]]
        if (!missing(groups)) {
          f1 <- rep(f, length=length(z1))
          z1 <- z1[f1==i]
          f2 <- rep(f, length=length(z2))
          z2 <- z2[f2==i]
        }
        l <- min(length(z1), length(z2))
        if (l <= 0L) next
        if (horizontal) {
          arrows(z1[1:l], at[k], z2[1:l], at[k+1],
            length=0, lty=lp$lty[k], lwd=lp$lwd[k], col=lp$col[k])
        } else {
          arrows(at[k], z1[1:l], at[k+1], z2[1:l],
            length=0, lty=lp$lty[k], lwd=lp$lwd[k], col=lp$col[k])
        }
      }
    }
  }

  if (points) {
    cxy <- par("cxy")
    offset <- offset * cex * if (horizontal) cxy[2L]/3 else cxy[1L]/2

    for (i in seq_len(m)) {
      k <- (i - 1) * n
      for (j in seq_len(n)) {
        k <- k + 1
        z1 <- x[[j]]
        if (!missing(groups)) {
          f1 <- rep(f, length=length(z1))
          z1 <- z1[f1==i]
        }
        l <- sum(!is.na(z1))
        if (l <= 0L) next
        if (l <= 1L || method == 1L) { # overplot
          d <- rep(0, length(z1))
        } else if (method == 2L) { # stack
          z2 <- split(z1, factor(z1))
          z1 <- unlist(z2, use.names = FALSE)
          d <- unlist(lapply(z2, seqs, side[k]), use.names = FALSE) * offset
          if (!is.na(amount)) {
            r <- max(abs(d), na.rm=TRUE)
            if (0 < r)
              d <- d / r * amount
          }
        } else {
          if (method == 3L) { # jitter
            d <- amount      
          } else { # violin
            d <- do.call(getFromNamespace("density", "stats"), 
              c(list(z1), density.par))
            d <- stats::approx(d$x, d$y, z1)$y
            d <- amount * d / max(d, na.rm=TRUE)
          }
          r <- if (side[k]==1L) c(-1,0) else if (side[k]==2L) c(-1,1) else c(0,1)
          d <- d * stats::runif(length(z1), r[1], r[2])
        }
        if (horizontal) {
          points(z1, at[k] + d,
            pch=pp$pch[k], cex=pp$cex[k], col=pp$col[k], bg=pp$bg[k])
        } else {
          points(at[k] + d, z1,
            pch=pp$pch[k], cex=pp$cex[k], col=pp$col[k], bg=pp$bg[k])
        }
      }
    }
  }

  if (!add) {
    if (missing(xlab)) xlab <- if (horizontal) dlab else glab
    if (missing(ylab)) ylab <- if (horizontal) glab else dlab
    title(main=main, sub=sub, xlab=xlab, ylab=ylab)  
    if (frame.plot) box()
  }

  invisible(at)
}


sp.plot <- function(...) UseMethod("sp.plot") 

sp.plot.formula <- function(formula, data=NULL, n.cut=5, sort.vars=c("original","wa", "alphabetical"), subset=NULL, sp.scale=c("fixed", "free", "exaggerated"), xlab=NULL, ylab=NULL, sp.max=NULL, cex.max=10, as.table=TRUE, ...) {
  warning("Function sp.plot is deprecated.  It will be removed from the next version of rioja.")
  Terms <- terms(formula)
  form <- formula(paste("~", deparse(Terms[[3]])))
  mf <- model.frame(form, data=data) 
  sp.plot.default(y=eval.parent(Terms[[2]]), x=mf, n.cut=n.cut, sort.vars=sort.vars, subset=subset, sp.scale=sp.scale, xlab=xlab, ylab=ylab, sp.max=sp.max, cex.max=cex.max, as.table=as.table, ...)
}

sp.plot.default <- function(y, x, n.cut=5, sort.vars=c("original","wa", "alphabetical"), subset=NULL, sp.scale=c("fixed", "free", "exaggerated"), xlab=NULL, ylab=NULL, sp.max=NULL, cex.max=NULL, as.table=TRUE, ...) {
#  require(lattice, quietly=TRUE)
  nm <- substitute(x)
  warning("Function sp.plot is deprecated.  It will be removed from the next version of rioja.")
  if (is.null(subset))
     subset <- 1:ncol(y)
  .sp.plot.panel1 <- function(x, y, ...) {
     lim <- current.panel.limits()
     exag <- ""
     if (max(y, na.rm=TRUE) < lim$ylim[2]/10) {
        y <- y * 10
        exag <- "x10"
     }
     panel.xyplot(x, y, ...)
     if (nchar(exag) > 2) {
       rx <- diff(range(lim$xlim))
       ry <- diff(range(lim$ylim))
       xp <- lim$xlim[1] + rx / 20
       yp <- lim$ylim[2] - ry / 20
       panel.text(xp, yp, "x10", col="red", adj=0, cex=0.8)
     }
  }
  scale0mx <- function(x, mx, sp.max=max(x, na.rm=TRUE)) {
   	 p0 <- min(x, na.rm=TRUE)
   	 ra <- max(x, na.rm=TRUE) - p0
	   x1 <- (x) / ra * mx
	   x1 * x / sp.max
  }
  .sp.plot.panel2 <- function(x, y, cex, subscripts, sp.scale, cex.max, sp.max, ...) {
    exag <- ""
    if (sp.scale=="fixed") {
       sz <- scale0mx(cex[subscripts], cex.max, sp.max)
    } else {
      if (sp.scale=="exaggerated") {
        sz <- scale0mx(cex[subscripts], cex.max, sp.max)
        if (max(cex[subscripts], na.rm=TRUE) < sp.max/10) {
          sz <- sz * 10
          exag <- "x10"
        }
      } else {
         sz <- scale0mx(cex[subscripts], cex.max)
      }
    }
    panel.xyplot(x, y, cex=sz, ...)
    if (nchar(exag) > 2 | sp.scale=="free") {
      lim <- current.panel.limits()
      rx <- diff(range(lim$xlim))
      ry <- diff(range(lim$ylim))
      xp <- lim$xlim[1] + rx / 20
      yp <- lim$ylim[2] - ry / 20
      if (sp.scale=="free")
        panel.text(xp, yp, round(max(cex[subscripts], na.rm=TRUE), 1), col="red", adj=0, cex=0.8)
      else
        panel.text(xp, yp, "x10", col="red", adj=0, cex=0.8)
    }
  }
  y <- y[, subset, drop=FALSE]
  if (class(x) == "numeric") {
     nm <- deparse(substitute(x))
     x <- as.matrix(x, ncol=1)
     colnames(x) <- nm
  }
  if (ncol(x) > 2)
     stop("Too many x-variables: give one or two")
  x1 <- x[, 1]
  if (is.null(xlab))
    xlab <- colnames(x)[1]
  if (ncol(x) == 2) {
    x2 <- x[, 2]
    if (is.null(ylab))
       ylab <- colnames(x)[2]
  } else {   
    x2 <- NULL
  } 
  if (nrow(y) != length(as.numeric(x1)))
     stop("Number of rows do not match in x & y")
   n <- colSums(y>0)  
   if (any(n < n.cut))
      y <- y[, n >= n.cut]
   sort.vars <- match.arg(sort.vars)
   sp.scale <- match.arg(sp.scale)
   if (sort.vars == "original") {
      or1 <- order(colnames(y))
      ord <- order(or1)
   } else {
      if (sort.vars == "wa") {
          or1 <- order(colnames(y))
          wa.sc <- apply(y[, or1], 2, function(x, env) { sum(x*env, na.rm=TRUE) / sum(x, na.rm=TRUE) }, env=x1)
          ord <- order(wa.sc)
      }
      else {
         ord <- 1:ncol(y)
      }
   }
   res <- stack(y)
   res$x1 <- rep(x1, by=ncol(y))
   if (is.null(sp.max))
     sp.max <- max(y, na.rm=TRUE)
   if (is.null(x2)) {
     ylim = range(res$value, na.rm=TRUE)
     ylim[2] <- sp.max
     if (sp.scale=="free")
        xyplot(values ~ x1 | ind, data=res, ylab=ylab, xlab=xlab, scales=list(x=list(alternating=FALSE), y="free"), as.table=as.table, index.cond=list(ord), ...)
     else if (sp.scale=="exaggerated")
        xyplot(values ~ x1 | ind, data=res, ylab=ylab, xlab=xlab, ylim=ylim, panel=.sp.plot.panel1, as.table=as.table, index.cond=list(ord), ...)
     else
        xyplot(values ~ x1 | ind, data=res, ylab=ylab, xlab=xlab, ylim=ylim, as.table=as.table, index.cond=list(ord), ...)
   } else {
      res$x2 <- rep(x2, by=ncol(y))
      if (is.null(cex.max))
        cex.max <- 10
      xyplot(x2 ~ x1 | ind, data=res, ylab=ylab, xlab=xlab, cex=res$value, panel=.sp.plot.panel2, as.table=as.table, index.cond=list(ord), sp.scale=sp.scale, cex.max=cex.max, sp.max=sp.max, ...)
   }
}



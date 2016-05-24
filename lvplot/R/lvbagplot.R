determine2dDepth <- function(n,k,alpha,perc) {
# by default, subtracts 3 levels from one dimensional solution
# --- no mathematical foundation for doing this so far

	k <- determineDepth(n,k,alpha,perc)-4
    if (k < 1)
        k <- 1
    return(k)
}

LVbagplot <- function(x, ...) UseMethod("LVbagplot",x)

LVbagplot.formula <- function(formula,alpha=0.95, k=NULL, perc=NULL, method="convex", horizontal=TRUE, xlab=NULL, ylab=NULL, col="grey30", bg="grey90", ...) {
    deparen <- function(expr) {
        while (is.language(expr) && !is.name(expr) && deparse(expr[[1]]) ==
            "(") expr <- expr[[2]]
        expr
    }
    bad.formula <- function() stop("invalid formula; use format y ~ x")
    bad.lengths <- function() stop("incompatible variable lengths")

    formula <- deparen(formula)
    if (!inherits(formula, "formula"))
        bad.formula()
    z <- deparen(formula[[2]])
    x <- deparen(formula[[3]])
    rhs <- deparen(formula[[3]])
    if (is.language(rhs) && !is.name(rhs) && (deparse(rhs[[1]]) ==
        "*" || deparse(rhs[[1]]) == "+")) {
        bad.formula()
    }
    z.name <- deparse(z)
    z <- eval(z,  parent.frame())
    x.name <- deparse(x)
    x <- eval(x,  parent.frame())

	LVbagplot.numeric(x,z, alpha, k, perc, col, method, xlab=x.name, ylab=z.name, ...)
}


LVbagplot.numeric <- function(x,y, alpha=0.95, k=NULL, perc=NULL, col="grey30", method="convex", horizontal=TRUE, xlab=NULL, ylab=NULL, bg="grey95", ...) {
    win <- function(dx, dy) {
        atan2(y = dy, x = dx)
    }

  if (missing(y)) {
  	print("don't have y")
  }
  args <- as.list(match.call())[-1]
  x.name <- args$x
  y.name <- args$y

  n <- length(x)
  if (length(y) != n) stop("x and y do not have the same length")

  k <- determine2dDepth(n,k,alpha,perc)
  src.col <- col

  if (length(src.col)==1) col <- color_scale(src.col, k)
  if ((length(src.col)==0)) col <- rep("grey",k)

  xy <- data.frame(x,y, order=1:length(x))

	if (method=="convex") {
	# compute halfspace depth
	  i <- 1

	  m <- nrow(xy)
	  res <- numeric(0)

	  while ((!is.null(m)) && (m > 0)) {
	    pts <- chull(xy)
	    res <- rbind(res, cbind(xy[pts,], rep(i,length(pts))))
	    xy <- xy[-pts,]
	    m <- dim(xy)[1]
	    i <- i+1
	  }
	}
	if (method=="apl") {

	}
	if (method=="depth") {
		dep <- vector(length=nrow(xy))
		for (i in 1:nrow(xy))
			dep[i] <- depth::depth(xy[i,],xy)
		res <- cbind(xy,round(dep*nrow(xy)))
	}
	# compute median as average of points with maximal halfspace depth
		med <- res[which(res[,3]==max(res[,3])),]
		medx <- mean(med[,1])
		medy <- mean(med[,2])

	# draw LV polygons
		if(is.null(xlab)) xlab=x.name
		if(is.null(ylab)) ylab=y.name
		plot(x,y, type="n", xlab=xlab, ylab=ylab)
	    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = bg)
	    box()
	    axis(1)
	    axis(2)
	    grid(lty=1, col="white")

		for (i in k:1) {
			dd <- 2^i
			Q <- res[n/dd,3]
			tp <- res[which(res[,3]==Q),1:2]
			angle <- win(tp[,1]-medx, tp[,2]-medy)
			ord <- order(angle)
			polygon(tp[ord,], col=col[i], density=-1)
		}
		points(medx, medy, pch=20)
		Qmin <- res[n/2^k,3]

		points(res[which(res[,3]>=Qmin),1:2], pch=".")
		points(res[which(res[,3]<Qmin),1:2], ...)

  invisible(list(outliers=res[which(res[,3]<Qmin),1:2]))
}

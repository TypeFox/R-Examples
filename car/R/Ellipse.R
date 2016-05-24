# Ellipses (orignally by J. Fox and G. Monette)

# added grid lines, 25 May 2010 by S. Weisberg
# arguments more consistent with other functions; ... passes args to plot, 5 Sept 2010 by J. Fox
# confidenceEllipse.lm and .glm can add to current plot, applying patch from Rafael Laboissiere, 17 Oct 2010 by J. Fox
# added fill and fill.alpha arguments for translucent fills (suggested by Michael Friendly), 14 Nov 2010 by J. Fox
# modified 2 May 2011 by Michael Friendly
#   - allow pivot=TRUE (with warning)
#   - barf on non-symmetric shape
#   - return coordinates of ellipse invisibly
# dataEllipse() and confidenceEllipse() invisibly return coordinates,  3 May 2011 by J. Fox
# Modified 5 May 2011 by Michael Friendly
#   - dataEllipse now honors add=FALSE, plot.points=FALSE
# Modified 16 May 2011 by Michaell Friendly
#   - corrected bug introduced in dataEllipse via allowing pivot=TRUE 
# Modified 7 Aug 2011 by J. Fox: added draw argument
# Modified 28 Nov 2011  by J. Fox (suggested by Michael Friendly):
#   - corrected bug in xlab, ylab in confidenceEllipse()
#   - added dfn argument to .lm and .glm methods for confidenceEllipse()
# Modified 14&16 Dec 2011 by J. Fox (suggested by Michael Friendly) to add weights argument to dataEllipse().
# Modified 2 Feb 2012 by J. Fox: Improved handling of center.pch argument to ellipse() (suggestion of Rob Kushler).
# 16 July 2012 added showLabels to dataEllipse
# 2014-02-16: prevent dataEllipse() from opening a graphics device when draw=FALSE (fixing bug reported by Rafael Laboissiere).
# 2015-09-04: throw error if there are too few colors for groups (fixing bug reported by Ottorino Pantani). J. Fox

ellipse <- function(center, shape, radius, log="", center.pch=19, center.cex=1.5, segments=51, draw=TRUE, add=draw, 
		xlab="", ylab="", col=palette()[2], lwd=2, fill=FALSE, fill.alpha=0.3,
		grid=TRUE, ...) {
	trans.colors <- function(col, alpha=0.5, names=NULL) {
		# this function by Michael Friendly
		nc <- length(col)
		na <- length(alpha)
		# make lengths conform, filling out to the longest
		if (nc != na) {
			col <- rep(col, length.out=max(nc,na))
			alpha <- rep(alpha, length.out=max(nc,na))
		}
		clr <-rbind(col2rgb(col)/255, alpha=alpha)
		col <- rgb(clr[1,], clr[2,], clr[3,], clr[4,], names=names)
		col
	}
	logged <- function(axis=c("x", "y")){
		axis <- match.arg(axis)
		0 != length(grep(axis, log))
	}
	if (! (is.vector(center) && 2==length(center))) stop("center must be a vector of length 2")
	if (! (is.matrix(shape) && all(2==dim(shape)))) stop("shape must be a 2 by 2 matrix")
	if (max(abs(shape - t(shape)))/max(abs(shape)) > 1e-10) stop("shape must be a symmetric matrix")
	angles <- (0:segments)*2*pi/segments 
	unit.circle <- cbind(cos(angles), sin(angles)) 
#	ellipse <- t(center + radius*t(unit.circle %*% chol(shape,pivot=TRUE))) 
	Q <- chol(shape, pivot=TRUE)
	order <- order(attr(Q, "pivot"))
	ellipse <- t( center + radius*t( unit.circle %*% Q[,order]))
	colnames(ellipse) <- c("x", "y")
	if (logged("x")) ellipse[, "x"] <- exp(ellipse[, "x"])
	if (logged("y")) ellipse[, "y"] <- exp(ellipse[, "y"])
	fill.col <- trans.colors(col, fill.alpha)
	if (draw) {
		if (add) {
			lines(ellipse, col=col, lwd=lwd, ...) 
			if (fill) polygon(ellipse, col=fill.col, border=NA)
		}
		else {
			plot(ellipse, type="n", xlab = xlab, ylab = ylab, ...) 
			if(grid){
				grid(lty=1, equilogs=FALSE)
				box()}
			lines(ellipse, col=col, lwd=lwd, ... )
			if (fill) polygon(ellipse, col=fill.col, border=NA)
		} 	
		if ((center.pch != FALSE) && (!is.null(center.pch))) points(center[1], center[2], pch=center.pch, cex=center.cex, col=col)
	}
	invisible(ellipse)
}


dataEllipse <- function(x, y, groups,
    group.labels=group.levels, ellipse.label,
    weights, log="", levels=c(0.5, 0.95), center.pch=19, 
    center.cex=1.5, draw=TRUE,
    plot.points=draw, add=!plot.points, segments=51, robust=FALSE, 
    xlab=deparse(substitute(x)), ylab=deparse(substitute(y)), 
    col=if (missing(groups)) palette()[1:2] else palette()[1:length(group.levels)],
    pch=if (missing(groups)) 1 else seq(group.levels),
    lwd=2, fill=FALSE, fill.alpha=0.3, grid=TRUE,
    labels, id.method = "mahal", id.n = if(id.method[1]=="identify") Inf else 0,
    id.cex=1, id.col=if (missing(groups)) palette()[1] else palette()(1:length(groups)), 
    ...) {
    label.ellipse <- function(ellipse, label, col, ...){
        # This sub-function from Michael Friendly
        if (cor(ellipse)[1,2] >= 0){         # position label above top right
            index <- which.max(ellipse[,2])
            x <- ellipse[index, 1] + 0.5 * strwidth(label)
            y <- ellipse[index, 2] + 0.5 * strheight("A")
            adj <- c(1, 0) 
        }
        else {                               # position label below bot left
            index <- which.min(ellipse[,2])
            x <- ellipse[index, 1] - 0.5 * strwidth(label)
            y <- ellipse[index, 2] - 0.5 * strheight("A")
            adj <- c(0, 1) 
        }
        text(x, y, label, adj=adj, col=col, ...)
    }
    if(missing(y)){
        if (is.matrix(x) && ncol(x) == 2) {
            if (missing(xlab)) xlab <- colnames(x)[1]
            if (missing(ylab)) ylab <- colnames(x)[2]
            y <- x[,2]
            x <- x[,1]
        }
        else stop("x and y must be vectors, or x must be a 2 column matrix")
    }
    else if(!(is.vector(x) && is.vector(y) && length(x) == length(y)))
        stop("x and y must be vectors of the same length")
    if (missing(weights)) weights <- rep(1, length(x))
    if (length(weights) != length(x)) stop("weights must be of the same length as x and y")
    if (!missing(groups)){
        xlab
        ylab
        if (!is.factor(groups)) stop ("groups must be a factor")
        if (!(length(groups) == length(x))) stop ("groups, x, and y must all be of the same length")
        if(missing(labels)) labels <- seq(length(x))
        valid <- complete.cases(x, y, groups)
        x <- x[valid]
        y <- y[valid]
        weights <- weights[valid]
        groups <- groups[valid]
        labels <- labels[valid]
        group.levels <- levels(groups)
        col <- col[!is.na(col)]
        if (length(col) < length(group.levels)) stop("too few colors for number of groups")
        result <- vector(length(group.levels), mode="list")
        names(result) <- group.levels
        if(draw) {
            if (!add) {
                plot(x, y, type="n", xlab=xlab, ylab=ylab,  ...) 
                if(grid){
                    grid(lty=1, equilogs=FALSE)
                    box()
                }
            }
        }
        for (lev in 1:length(group.levels)){
            level <- group.levels[lev]
            sel <- groups == level
            result[[lev]] <- dataEllipse(x[sel], y[sel],
                weights=weights[sel], log=log, levels=levels, center.pch=center.pch,
                center.cex=center.cex, draw=draw, plot.points=plot.points, add=TRUE, segments=segments,
                robust=robust, col=rep(col[lev], 2), pch=pch[lev], lwd=lwd, fill=fill, fill.alpha=fill.alpha, 
                labels=labels[sel], id.method=id.method, id.n=id.n, id.cex=id.cex, id.col=col[lev], 
                ellipse.label=group.labels[lev], ...)
        }
        return(invisible(result))
    }
    if (length(col) == 1) col <- rep(col, 2)
    if(draw) {
        if (!add) {
            plot(x, y, type="n", xlab=xlab, ylab=ylab,  ...) 
            if(grid){
                grid(lty=1, equilogs=FALSE)
                box()}
        }
        if (plot.points)  points(x, y, col=col[1], pch=pch[1], ...)
    }
    dfn <- 2
    dfd <- length(x) - 1
    if (robust) {
        use <- weights > 0
        v <- cov.trob(cbind(x[use], y[use]), wt=weights[use])
        shape <- v$cov
        center <- v$center
    }
    else {
        v <- cov.wt(cbind(x, y), wt=weights)
        shape <- v$cov
        center <- v$center
    }
    result <- vector("list", length=length(levels))
    names(result) <- levels
    for (i in seq(along=levels)) {
        level <- levels[i]
        radius <- sqrt(dfn * qf(level, dfn, dfd ))
        result[[i]] <- ellipse(center, shape, radius, log=log,
            center.pch=center.pch, center.cex=center.cex, segments=segments, 
            col=col[2], lwd=lwd, fill=fill, fill.alpha=fill.alpha, draw=draw, ...)
        if (!missing(ellipse.label)) {
            lab <- if (length(ellipse.label) < i) ellipse.label[1] else ellipse.label[i]
            label.ellipse(result[[i]], lab, col[2], ...)
        }
    }
    if (missing(labels)) labels <- seq(length(x))
    if (draw) showLabels(x, y, labels=labels,
        id.method=id.method, id.n=id.n, id.cex=id.cex,
        id.col=id.col)
    invisible(if (length(levels) == 1) result[[1]] else result)
}
    
confidenceEllipse <- function (model, ...) {
	UseMethod("confidenceEllipse")
}


confidenceEllipse.lm <- function(model, which.coef, L, levels=0.95, Scheffe=FALSE, dfn,
		center.pch=19, center.cex=1.5, segments=51, xlab, ylab, 
		col=palette()[2], lwd=2, fill=FALSE, fill.alpha=0.3, draw=TRUE, add=!draw, ...){
	if (missing(dfn)) dfn <- if (Scheffe) sum(df.terms(model)) else 2
	dfd <- df.residual(model)
	if (missing(L)){
		which.coef <- if(length(coefficients(model)) == 2) c(1, 2)
				else{
					if (missing(which.coef)){
						if (has.intercept(model)) c(2,3) else c(1, 2)
					} else which.coef
				}
		coef <- coefficients(model)[which.coef]
		if (missing(xlab)) xlab <- paste(names(coef)[1], "coefficient")
		if (missing(ylab)) ylab <-  paste(names(coef)[2], "coefficient")
		shape <- vcov(model)[which.coef, which.coef]
	}
	else {
		res <- makeLinearCombinations(L, coef(model), vcov(model))
		coef <- res$coef
		xlab <- res$xlab
		ylab <- res$ylab
		shape <- res$shape
	}
	levels <- rev(sort(levels))
	result <- vector("list", length=length(levels))
	names(result) <- levels
	for (i in seq(along=levels)){
		level <- levels[i]
		radius <- sqrt(dfn*qf(level, dfn, dfd))
		add.plot <- !level==max(levels) | add
		result[[i]] <- ellipse(coef, shape, radius, add=add.plot, xlab=xlab, ylab=ylab,
				center.pch=center.pch, center.cex=center.cex, segments=segments, 
				col=col, lwd=lwd, fill=fill, fill.alpha=fill.alpha, draw=draw, ...)
	}
	invisible(if (length(levels) == 1) result[[1]] else result)
}


confidenceEllipse.default <- function(model, which.coef, L, levels=0.95, Scheffe=FALSE, dfn,
		center.pch=19, center.cex=1.5, segments=51, xlab, ylab,
		col=palette()[2], lwd=2, fill=FALSE, fill.alpha=0.3, draw=TRUE, add=!draw, ...){
	if (missing(L)){
		which.coef <- if(length(coefficients(model)) == 2) c(1, 2)
				else{
					if (missing(which.coef)){
						if (has.intercept(model)) c(2, 3) else c(1, 2)
					} else which.coef
				}
		coef <- coefficients(model)[which.coef]
		shape <- vcov(model)[which.coef, which.coef]
		xlab <- if (missing(xlab)) paste(names(coef)[1], "coefficient")
		ylab <- if (missing(ylab)) paste(names(coef)[2], "coefficient")
	}
	else {
		res <- makeLinearCombinations(L, coef(model), vcov(model))
		coef <- res$coef
		xlab <- res$xlab
		ylab <- res$ylab
		shape <- res$shape
	}
	df <- if (!missing(dfn)) dfn
			else if (Scheffe) sum(df.terms(model)) else 2
	levels <- rev(sort(levels))
	result <- vector("list", length=length(levels))
	names(result) <- levels
	for (i in seq(along=levels)){
		level <- levels[i]
		radius <- sqrt(qchisq(level, df))
		add.plot <- !level==max(levels) | add
		result[[i]] <- ellipse(coef, shape, radius, add=add.plot, xlab=xlab, ylab=ylab,
				center.pch=center.pch, center.cex=center.cex, segments=segments,
				col=col, lwd=lwd, fill=fill, fill.alpha=fill.alpha, draw=draw, ...)
	}
	invisible(if (length(levels) == 1) result[[1]] else result)
}

confidenceEllipse.glm <- function (model, chisq, ...) {
	sumry <- summary(model)
	if (missing(chisq)) chisq <- is.null(sumry$dispersion)
	if (chisq) confidenceEllipse.default(model, ...)
	else confidenceEllipse.lm(model, ...)
}

makeLinearCombinations <- function(L, coef, V){
	nms <- names(coef)
	if (is.character(L)){
		L <- makeHypothesis(nms, L)
		L <- L[, -ncol(L)]
	}
	if (nrow(L) != 2 || ncol(L) != length(coef))
		stop("the hypothesis matrix is the wrong size")
	coef <- as.vector(L %*% coef)
	shape <- L %*% V %*% t(L)
	L.nms <- printHypothesis(L, c(0, 0), nms)
	names(coef) <- sub(" =.*", "", L.nms)
	xlab <- names(coef)[1]
	ylab <- names(coef)[2]
	list(coef=coef, shape=shape, xlab=xlab, ylab=ylab)
}

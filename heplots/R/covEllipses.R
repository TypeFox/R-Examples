covEllipses <-function(x, ...) {
	UseMethod("covEllipses")
}


covEllipses.boxM <-
		function(x, ...) {

	cov <- c(x$cov, pooled=list(x$pooled))
	mns <- x$means
	df <- x$df
	covEllipses.default(cov, mns, df, ...)
}

covEllipses.data.frame <-
		function(x, group,
		         pooled=TRUE, 
		         method = c("classical", "mve", "mcd"), ...) {

 method <- match.arg(method)
 
 if (missing(group)) {
   group <- factor(rep(1, nrow(x)))
   pooled <- FALSE
 }

 if (!is.factor(group)) {
   warning(deparse(substitute(group)), " coerced to factor.")
   group <- as.factor(group)
 }
 
 p <- ncol(x)
 nlev <- nlevels(group)
 lev <- levels(group)
 dfs <- tapply(group, group, length) - 1
 mats <- list()
 means <- matrix(0, nrow=nlev, ncol=p) 
 for(i in 1:nlev) {
 		rcov <- MASS::cov.rob(x[group == lev[i], ], method=method)
    mats[[i]] <- rcov$cov
    means[i,] <- rcov$center
 }

 names(mats) <- lev
 rownames(means) <- lev
 colnames(means) <- colnames(x)

 if(pooled) {
	 rcov <- MASS::cov.rob(x, method=method)
	 pooled <- rcov$cov
	 mats <- c(mats, pooled=list(pooled))
	 means <- rbind(means, pooled=rcov$center)
	 dfs <- c(dfs, sum(dfs))
 }

   covEllipses.default(mats, means, dfs, ...)

}

covEllipses.matrix <- covEllipses.data.frame

		  
covEllipses.default <-
		function ( 
		    x,                 # a list of covariance matrices
		    means,             # a matrix of means
		    df,                # vector of degrees of freedom
		    labels=NULL,
				variables=1:2,     # x,y variables for the plot [variable names or numbers]
				level=0.68,
				segments=40,       # line segments in each ellipse
				center = FALSE,    # center the ellipses at c(0,0)?
				center.pch="+",  
				center.cex=2,
				col=getOption("heplot.colors", c("red", "blue", "black", "darkgreen", "darkcyan","magenta", "brown","darkgray")),
				# colors for ellipses
				lty=1,
				lwd=2,
				fill=FALSE,         ## whether to draw filled ellipses (vectorized)
				fill.alpha=0.3,     ## alpha transparency for filled ellipses
				label.pos=0,    # label positions: NULL or 0:4
				xlab,
				ylab,
				main="",
				xlim,           # min/max for X (override internal min/max calc) 
				ylim,
				axes=TRUE,      # whether to draw the axes
				offset.axes,    # if specified, the proportion by which to expand the axes on each end (e.g., .05)
				add=FALSE,      # add to existing plot?
				warn.rank=FALSE,  
				...) 
{

	ell <- function(center, shape, radius) {
		angles <- (0:segments)*2*pi/segments
		circle <- radius * cbind( cos(angles), sin(angles))
		if (!warn.rank){
			warn <- options(warn=-1)
			on.exit(options(warn))
		}
		Q <- chol(shape, pivot=TRUE)
		order <- order(attr(Q, "pivot"))
		t( c(center) + t( circle %*% Q[,order]))
	}

	if (!is.list(x)) stop("Argument 'x' must be a list of covariance matrices")
	cov <- x
	response.names <- colnames(cov[[1]])
	p <- ncol(cov[[1]])

	if (!is.numeric(variables)) {
		vars <- variables
		variables <- match(vars, response.names)
		check <- is.na(variables)
		if (any(check)) stop(paste(vars[check], collapse=", "), 
					" not among response variables.") 
	}
	else {
		if (any (variables > p)) stop("There are only ", 	p, " response variables among", variables)
		vars <- response.names[variables]
	}
	n.ell <- length(cov)
	if (n.ell == 0) stop("Nothing to plot.")
	if (n.ell != nrow(means)) 
	    stop( paste0("number of covariance matrices (", n.ell, ") does not conform to rows of means (", nrow(means), ")") )
	if (n.ell != length(df)) 
	  stop( paste0("number of covariance matrices (", n.ell, ") does not conform to df (", length(df), ")") )
	
	if (missing(xlab)) xlab <- vars[1]
	if (missing(ylab)) ylab <- vars[2] 

	# assign colors and line styles
	rep_fun <- rep_len
	col <- rep_fun(col, n.ell) 
	lty <- rep_fun(lty, n.ell)
	lwd <- rep_fun(lwd, n.ell)
	# handle filled ellipses
	fill <- rep_fun(fill, n.ell)
	fill.alpha <- rep_fun(fill.alpha, n.ell)
	fill.col <- trans.colors(col, fill.alpha)
	label.pos <- rep_fun(label.pos, n.ell)
	fill.col <- ifelse(fill, fill.col, NA)

	radius <- c(sqrt(2 * qf(level, 2, df)))
	ellipses <- as.list(rep(0, n.ell))
	for(i in 1:n.ell) {
		S <- as.matrix(cov[[i]])
		S <- S[vars, vars]
		ctr <- if (center)  c(0,0)
		       else as.numeric(means[i, vars])
		ellipses[[i]] <- ell(ctr, S, radius[i])
	}
	
	if (!add){
		max <- apply(sapply(ellipses, function(X) apply(X, 2, max)), 1, max)
		min <- apply(sapply(ellipses, function(X) apply(X, 2, min)), 1, min)

		if (!missing(offset.axes)){
			range <- max - min
			min <- min - offset.axes*range
			max <- max + offset.axes*range
		}
		xlim <- if(missing(xlim)) c(min[1], max[1]) else xlim
		ylim <- if(missing(ylim)) c(min[2], max[2]) else ylim
		plot(xlim, ylim,  type = "n", xlab=xlab, ylab=ylab, main=main, axes=axes, ...)
	}

	labels <- if (!is.null(labels)) labels
			else names(cov)

	for (i in 1:n.ell){
			polygon(ellipses[[i]], col=fill.col[i], border=col[i],  lty=lty[i], lwd=lwd[i])
			label.ellipse(ellipses[[i]], labels[i], col=col[i], label.pos=label.pos[i], ...) 
			if (!center) 
				points(means[i,1], means[i,2], pch=center.pch, cex=center.cex, col=col[i], xpd=TRUE)
		}   

	names(ellipses) <- labels
#	result <- if (!add) list(ellipses, center=means, xlim=xlim, ylim=ylim, radius=radius)
#			else list(H=ellipses,  center=gmean, radius=radius)
	result <- ellipses
	class(result) <- "covEllipses"
	invisible(result)
	
}


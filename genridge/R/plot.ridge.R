# for pcaridge objects, default to last 2 variables
plot.pcaridge <-
		function(x, variables=(p-1):p, labels=NULL, ...) {
	p <- dim(coef(x))[2]
	plot.ridge(x, variables, labels=labels, ...)
}

plot.ridge <-
function(x, variables=1:2, radius=1, which.lambda=1:length(x$lambda), 
		labels=lambda, pos=3, cex=1.2,
		lwd=2, lty=1, xlim, ylim,
		col = c("black", "red", "darkgreen", "blue","darkcyan","magenta", "brown","darkgray"), 
		center.pch = 16, center.cex=1.5,
		fill=FALSE, fill.alpha=0.3, ref=TRUE, ref.col=gray(.70), ...) {

	ell <- function(center, shape, radius, segments=60) {
		angles <- (0:segments)*2*pi/segments
		circle <- radius * cbind( cos(angles), sin(angles))
		warn <- options(warn=-1)
		on.exit(options(warn))
		Q <- chol(shape, pivot=TRUE)
		order <- order(attr(Q, "pivot"))
		t( c(center) + t( circle %*% Q[,order]))
	}

	vnames <- dimnames(x$coef)[[2]]
	if (!is.numeric(variables)) {
		vars <- variables
		variables <- match(vars, vnames)
		check <- is.na(variables)
		if (any(check)) stop(paste(vars[check], collapse=", "), 
					" not among the predictor variables.") 
	}
	else {
		if (any (variables > length(vnames))) stop("There are only ", 
					length(vnames), " predictor variables.")
		vars <- vnames[variables]
	}
	if(length(variables)>2) {
		warning("Only two variables will be plotted. Perhaps you want plot3d.ridge()")
		variables <- variables[1:2]
	}
	lambda <- x$lambda[which.lambda]
	coef <- x$coef[which.lambda,variables]
	cov <- x$cov[which.lambda]
	n.ell <- length(lambda)
	lambda <- signif(lambda, 3)   # avoid many decimals when used as labels

	ells <- as.list(rep(0, n.ell))
# generate the ellipses for each lambda, to get xlim & ylim
	for (i in 1:n.ell) {
		ells[[i]] <- ell(center=coef[i,], shape=(cov[[i]])[variables,variables], radius=radius)
	}
	max <- apply(sapply(ells, function(X) apply(X, 2, max)), 1, max)
	min <- apply(sapply(ells, function(X) apply(X, 2, min)), 1, min)
	xlim <- if(missing(xlim)) c(min[1], max[1]) else xlim
	ylim <- if(missing(ylim)) c(min[2], max[2]) else ylim

	col <- rep(col, n.ell)		
	lwd <- rep(lwd, n.ell)		
	lty <- rep(lty, n.ell)		
	# handle filled ellipses
	fill <- rep(fill, n.ell)
	fill.alpha <- rep(fill.alpha, n.ell)
	fill.col <- trans.colors(col, fill.alpha)
	fill.col <- ifelse(fill, fill.col, NA)

	plot(coef, type='b', pch=center.pch, cex=center.cex, col=col, xlim=xlim, ylim=ylim, ...)
	if (ref) abline(v=0, h=0, col=ref.col)
	for (i in 1:n.ell) {
#		lines(ells[[i]], col=col[i], lwd=lwd[i], lty=lty[i])
		polygon(ells[[i]], col=fill.col[i], border=col[i],  lty=lty[i], lwd=lwd[i])
	}
	if(!is.null(labels)) text(coef, labels=labels, pos=pos, cex=cex)
}


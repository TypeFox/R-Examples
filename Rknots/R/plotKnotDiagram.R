##Script generated in:
# 2011
# 9:25:11 AM
#by: 
# Author: Maurizio Rinaldi @ University of Piemonte Orientale
# Author: Federico Comoglio @ D-BSSE, ETH Zurich
###############################################################################

defineCut <- function(points3D, size) {	
	cut <- min( diff( apply(points3D, 2, range) ) ) / size
	return(cut)
}

splitUndercrossing <- function(points2D, i, cut)
{
	r1 <- norM(points2D[i - 1, ] - points2D[i, ])
	r2 <- norM(points2D[i + 1, ] - points2D[i, ])
	radius <- min(c(r1, r2, cut))
	V1 <- points2D[i - 1, ] - points2D[i, ]
	V2 <- points2D[i + 1, ] - points2D[i, ]
	matrix(c(points2D[i, ] + 0.5 * radius * (V1) / norM(V1),
			 points2D[i,] + 0.5 * radius * (V2) / norM(V2)),
				nrow = 2, byrow = TRUE)
}

plotDiagram <- function(points3D, ends, pca = FALSE, size = 1, colors = c(), return.vars = FALSE, ...) {
	#if PCA
	if(pca) points3D <- PCAProjection( points3D )
	vp <- vertexPresentation(points3D, ends)
	if( is.null(vp) ) ##no intersections
		return(plot(size * points3D[, 1:2], xlab = '', ylab = '', type = 'l', axes = FALSE, ...))
	
	p3 <- t( sapply(vp$points3Dout, "[[", 1) )
	points2D <- p3[, 1:2]
	n <- nrow( points2D )
	endsout <- vp$endsout 
	extends <- c(0, endsout, n)
	
	eAT <- auxiliaryAlexander(get2D(vp$points3Dout), vp$ends)
	under <- eAT[[7]] 
	over <- eAT[[5]]
	n.under <- length(under)
	comp <- sapply(1 : n, function(x) -1 + localize(x, extends))
	
	if(is.null(colors)) {
		colors <- 1 : max(comp) #not supplied by the user
		colors.set <- comp
	}
	else 
		colors.set <- rep( colors, times = diff(extends) )
	
	cut.k <- defineCut(p3, size)
	
	plot(size * points2D, xlab = '', ylab = '', type = 'n', axes = FALSE, ...)
	
	split.under <- lapply(under, function(x) splitUndercrossing(points2D, x, cut.k))
	dd <- setdiff(setdiff( seq(-1 + n), c(under - 1, under) ), endsout)
	n.max <- 0
	for (i in dd) 
		lines(points2D[i : (i + 1), 1], points2D[i : (i + 1), 2], col = colors.set[i], ...)
	
	tmp <- lapply( c(unique(comp)), function(x) intersect(dd, which(comp == x)))
	normV <- lapply(tmp, function(j) 
				sapply(j, function(x) norM( points2D[x + 1, ] - points2D[x, ])) )
	longest.idx <- lapply(normV, which.max)
	longest.edge <- sapply( 1 : max(comp), function(x) tmp[[x]][longest.idx[[x]]])	
	for (i in seq(length(under))) {
		lines(c(split.under[[i]][2, 1], points2D[under[i] + 1, 1]),
				c(split.under[[i]][2, 2], points2D[under[i] + 1, 2]), col = colors.set[under[i]], lwd = 1)
		lines(c(split.under[[i]][1, 1], points2D[under[i] - 1, 1]),
				c(split.under[[i]][1, 2], points2D[under[i] - 1, 2]), col = colors.set[under[i]], lwd = 1)
	}
	for(i in 1 : max(comp)) {
		s <- longest.edge[i]
		arrows(points2D[s, 1], points2D[s, 2],
				.5 * (points2D[s, 1] + points2D[s + 1, 1]),
				.5 * (points2D[s, 2] + points2D[s + 1, 2]), length = 0.1, col = colors[i], ...)
	}
	
	if(return.vars)
		return(list(points3D = p3, ends = endsout, undercross = under, overcross = over, component = comp) )
}












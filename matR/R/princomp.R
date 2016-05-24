
#-----------------------------------------------------------------------------------------
#  principal coordinates analysis of biom object.
#-----------------------------------------------------------------------------------------

princomp.biom <- function(
		x,
		method="euclidean",
		dim=1:3,
		..., 
		map=NULL, 
		rows=TRUE, 
		columns=TRUE, 
		rerender=NULL) {

	x <- x [rows, columns]
#	method <- match.arg (method)						# need to resolve this
	arg <- list (...)

	if (inherits (rerender, "pco")) {
		P <- rerender
	} else {
		if (inherits (rerender, "dist")) {
			D <- rerender
		} else if (is.null (rerender)) {
			D <- distx (x, method = method)
		} else
			stop ("\'rerender\' is of an unsupported class")

		P <- ecodist::pco(D)
		scaled <- P$values / sum(P$values)
		names (scaled) <- paste0 ("PCO", 1:length(scaled))
		rownames (P$vectors) <- colnames(x)
		P <- list (values = scaled, vectors = P$vectors, dist = D)
		}

####  apply metadata substitution to "labels"
####  apply metadata mapping
####  save arguments "labels.*" for later

	arg$labels <- subColumn (arg$labels, x)
	names (arg) [names (arg) == "labels"] <- "label.labels"

	arg [names (map)] <- parMap (x, map, arg)

	arg.plot <- arg [substr(names(arg),1,6) != 'label.']
	arg.labels <- arg [substr(names(arg),1,6) == 'label.']
	names (arg.labels) <- substring (names (arg.labels), 7)

#-----------------------------------------------------------------------------------------
#  one-dimensional plot
#-----------------------------------------------------------------------------------------
	if (length (dim) == 1) {
		par <- resolve (arg.plot, list(
			x = P$vectors [ ,dim],
			y = rep_len (0, ncol(x)),
			xlab = paste0 ("PC", dim, ", R^2 = ", format (P$values [dim], dig=3)),
			ylab = "",
			ylim = c(-0.5,1),
			yaxt = "n",
			bty = "n",
			pch = 19,
			cex = 0.7,
			col = "black"))
		do.call (plot, par)
		do.call (points, par)
		abline (h=0, lty="dotted")

#-----------------------------------------------------------------------------------------
#  two-dimensional plot
#-----------------------------------------------------------------------------------------
	} else if (length (dim) == 2) {
		par <- resolve (arg.plot, list(
			x = P$vectors [ ,dim [1]],
			y = P$vectors [ ,dim [2]],
			xlab = paste0 ("PC", dim[1], ", R^2 = ", format (P$values [dim[1]], dig = 3)),
			ylab = paste0 ("PC", dim[2], ", R^2 = ", format (P$values [dim[2]], dig = 3)),
			pch = 19,
			cex = 0.7,
			col = "black"))
		do.call (plot, par)
		do.call (points, par)
		grid ()

#-----------------------------------------------------------------------------------------
#  three-dimensional plot.
#-----------------------------------------------------------------------------------------
	} else if (length (dim) == 3) {

####  tweak to enable consistent use of "col" and "cex".
####  and a hack:  update x and y after the call for label placement (below).

		names (arg.plot) [names(arg.plot) == "col"] <- "color"
		names (arg.plot) [names(arg.plot) == "cex"] <- "cex.symbols"
		par <- resolve (arg.plot, list(
			x = P$vectors [ ,dim [1]],
			y = P$vectors [ ,dim [2]],
			z = P$vectors [ ,dim [3]],
			xlab = paste0 ("PC", dim[1], ", R^2 = ", format (P$values [dim[1]], dig=3)),
			ylab = paste0 ("PC", dim[2], ", R^2 = ", format (P$values [dim[2]], dig=3)),
			zlab = paste0 ("PC", dim[3], ", R^2 = ", format (P$values [dim[3]], dig=3)),
			pch = 19,
			cex.symbols = 0.7,
			color = "black",
			type = "h",
			lty.hplot = "dotted",
			axis = TRUE,
			box = FALSE))
		xy <- do.call (scatterplot3d::scatterplot3d, par) $ xyz.convert (par$x, par$y, par$z)
 		par [c('x','y')] <- xy [c('x','y')]
		}

####  add labels, now using only arguments beginning with "label."
####  this includes "labels" itself, per earlier tweak

	do.call (text, resolve (arg.labels, list(
		x = par$x,
		y = par$y,
		labels = colnames(x),
		cex = 0.7,
		pos = if (length (dim) == 1) ifelse (P$vectors [,dim] > 0, 2, 4) else 4,
		srt = if (length (dim) == 1) 60)))

####  structure return value

	P$call <- match.call()
	class (P) <- c("pco", "list")
	invisible (P)
	}

#' Plotting function for vector decomposition and remainder fields
#'
#' This function plots various vector fields
#' @param x.field a two-dimensional array containing the x-values for the vector field, generated from \code{\link{VecDecomAll}}, \code{\link{VecDecomVec}}, \code{\link{VecDecomGrad}}, or \code{\link{VecDecomRem}}.
#' @param y.field a two-dimensional array containing the y-values for the vector field, generated from \code{\link{VecDecomAll}}, \code{\link{VecDecomVec}}, \code{\link{VecDecomGrad}}, or \code{\link{VecDecomRem}}.
#' @param dens two-element vector respectively specifying the number of respective arrows in the x and y directions.
#' @param x.bound the x boundaries denoted at c(minimum, maximum) for the quasi-potential simulation.
#' @param y.bound the y boundaries denoted at c(minimum, maximum) for the quasi-potential simulation.
#' @param xlim numeric vectors of length 2, giving the x coordinate range.
#' @param ylim numeric vectors of length 2, giving the y coordinate range.
#' @param arrow.type sets the type of line segments plotted. If set to "proportional" the length of the line segments reflects the magnitude of the derivative. If set to "equal" the line segments take equal lengths, simply reflecting the gradient of the derivative(s). Defaults to "equal".
#' @param tail.length multiplies the current length of the tail (both proportional and equal arrow.types) by the specified factor.  The argument defaults to 1, which is length of the longest vector within the domain boundaries (i.e., the entire field).
#' @param head.length length of the edges of the arrow head (in inches).
#' @param ... passes arguments to both \code{\link{plot}}.
#' @details If \code{arrow.type = "proportional"}, a common warning, passed from \code{\link{arrows}}, will appear: "The direction of a zero-length arrow is indeterminate, and hence so is the direction of the arrowheads. To allow for rounding error, arrowheads are omitted (with a warning) on any arrow of length less than 1/1000 inch."  Either increase \code{tail.length} or increase the plot window to avoid this warning.
#' @keywords vector field plot, detrministic skeleton vector field plot, gradient vector field plot, remainder vector field plot
#'
#' @examples
#' # First, system of equations
#' 	equationx <- "1.54*x*(1.0-(x/10.14)) - (y*x*x)/(1.0+x*x)"
#' 	equationy <- "((0.476*x*x*y)/(1+x*x)) - 0.112590*y*y"
#' 
#' # Second, shared parameters for each quasi-potential run
#' 	xbounds <- c(-0.5, 10.0)
#' 	ybounds <- c(-0.5, 10.0)
#' 	xstepnumber <- 100
#' 	ystepnumber <- 100
#' 
#' # Third, first local quasi-potential run
#' 	xinit1 <- 1.40491
#' 	yinit1 <- 2.80808
#' 	storage.eq1 <- QPotential(x.rhs = equationx, x.start = xinit1, 
#'		x.bound = xbounds, x.num.steps = xstepnumber, y.rhs = equationy, 
#'		y.start = yinit1, y.bound = ybounds, y.num.steps = ystepnumber)
#' 
#' # Fourth, second local quasi-potential run
#' 	xinit2 <- 4.9040
#' 	yinit2 <- 4.06187
#' 	storage.eq2 <- QPotential(x.rhs = equationx, x.start = xinit2, 
#'		x.bound = xbounds, x.num.steps = xstepnumber, y.rhs = equationy, 
#'		y.start = yinit2, y.bound = ybounds, y.num.steps = ystepnumber)
#' 
#' # Fifth, determine global quasi-potential 
#' 	unst.x <- c(0, 4.2008)
#' 	unst.y <- c(0, 4.0039)
#' 	ex1.global <- QPGlobal(local.surfaces = list(storage.eq1, storage.eq2), 
#'		unstable.eq.x = unst.x, unstable.eq.y = unst.y, x.bound = xbounds, 
#'		y.bound = ybounds)
#' 
#' # Sixth, decompose the global quasi-potential into the 
#' # deterministic skeleton, gradient, and remainder vector fields
#' 	VDAll <- VecDecomAll(surface = ex1.global, x.rhs = equationx, y.rhs = equationy, 
#'		x.bound = xbounds, y.bound = ybounds)
#' 
#' # Seventh, plot all three vector fields
#' 	# The deterministic skeleton vector field
#' 	VecDecomPlot(x.field = VDAll[,,1], y.field = VDAll[,,2], dens = c(25,25), 
#'		x.bound = xbounds, y.bound = ybounds, tail.length = 0.25, head.length = 0.05)
#' 	# The gradient vector field
#' 	VecDecomPlot(x.field = VDAll[,,3], y.field = VDAll[,,4], dens = c(25,25), 
#'		x.bound = xbounds, y.bound = ybounds, tail.length = 0.15, head.length = 0.05)
#' 	# The remainder vector field
#' 	VecDecomPlot(x.field = VDAll[,,5], y.field = VDAll[,,6], dens = c(25,25), 
#'		x.bound = xbounds, y.bound = ybounds, tail.length = 0.15, head.length = 0.05)

VecDecomPlot <- function(x.field, y.field, dens, x.bound, y.bound, xlim = 'NULL', ylim = 'NULL', arrow.type = "equal", tail.length = 0.25, head.length = 0.25, ...){

		if(any(dim(x.field) != dim(y.field))){stop("x.field and y.field have at least one unequal dimension length")}
		x.range <- max(x.bound)-min(x.bound)
		y.range <- max(y.bound)-min(y.bound)
		row.range <- nrow(x.field)-1
		col.range <- ncol(x.field)-1

		if(missing(xlim) == F & missing(ylim) == F) {
		row.min <- (min(xlim)-min(x.bound))/x.range*row.range + 1
		row.max <- (max(xlim)-min(x.bound))/x.range*row.range + 1
		col.min <- (min(ylim)-min(y.bound))/y.range*col.range + 1
		col.max <- (max(ylim)-min(y.bound))/y.range*col.range + 1
		} else {
			if(missing(xlim)) {
				row.min <- min(which(x.field != 0 , arr.ind = T)[,1])
				row.max <- max(which(x.field != 0 , arr.ind = T)[,1])
				x.min <- ((row.min-1)/row.range)*x.range + min(x.bound)
				x.max <- ((row.max-1)/row.range)*x.range + min(x.bound)
				xlim <- c(x.min,x.max)
				}
			if(missing(ylim)) {
				col.min <- min(which(x.field != 0 , arr.ind = T)[,2])
				col.max <- max(which(x.field != 0 , arr.ind = T)[,2])
				y.min <- ((col.min-1)/col.range)*y.range + min(y.bound)
				y.max <- ((col.max-1)/col.range)*y.range + min(y.bound)
				ylim <- c(y.min,y.max)
				}
			}

 	sub.x <- seq(row.min, row.max, length.out=dens[1])
	sub.y <- seq(col.min, col.max, length.out=dens[2])

	sub.x.val <- ((sub.x-1)/row.range)*x.range + min(x.bound)
	sub.y.val <- ((sub.y-1)/col.range)*y.range + min(y.bound)

	dx.sub <- x.field[sub.x, sub.y]
	dy.sub <- y.field[sub.x, sub.y]

	if(arrow.type == "proportional"){
	dx.rel <- (dx.sub/max(((dx.sub^2)+(dy.sub^2))^0.5, na.rm = T))
	dy.rel <- (dy.sub/max(((dx.sub^2)+(dy.sub^2))^0.5, na.rm = T))
	dx.plot <- dx.rel*tail.length
	dy.plot <- dy.rel*tail.length
	}

	if(arrow.type == "equal"){
	dx.even <- dx.sub/((dx.sub^2)+(dy.sub^2))^0.5
	dy.even <- dy.sub/(((dx.sub^2)+(dy.sub^2))^0.5)
	dx.plot <- dx.even*tail.length
	dy.plot <- dy.even*tail.length
	}


	if(arrow.type != "proportional" & arrow.type != "equal"){
	dx.plot <- dx.sub*tail.length
	dy.plot <- dy.sub*tail.length
	}

 	qpr <- nrow(dx.plot)
	qpc <- ncol(dy.plot)

	plot(0 , type = "n" , xlim = xlim , ylim = ylim , las = 1, ...)
	for (j in 1:qpr){
			for (i in 1:qpc){
				x0 <- sub.x.val[j] - (dx.plot[j,i]/2)
				x1 <- sub.x.val[j] + (dx.plot[j,i]/2)
				y0 <- sub.y.val[i] - (dy.plot[j,i]/2)
				y1 <- sub.y.val[i] + (dy.plot[j,i]/2)
				arrows(x0, y0, x1, y1, length = head.length)
			}
		}
}

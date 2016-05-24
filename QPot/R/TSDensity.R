#' Density plot from simulation of two-dimensional stochastic differential equations
#'
#' This function creates density plots for the simulation of two-dimensional stochastic differential equations from \code{\link{TSTraj}}
#' @param mat a matrix output from \code{\link{TSTraj}}.
#' @param dim dimensions of the plot; \code{dim = 1} plots simple density histogram or \code{dim = 2} plots the density in state space (i.e., X and Y respectively on the abscissa and ordinate axes).
#' @param xlim numeric vectors of length 2, giving the x coordinate range. Default \code{= 'NULL'} automatically sizes plot window.
#' @param ylim numeric vectors of length 2, giving the y coordinate range. Default \code{= 'NULL'} automatically sizes plot window.
#' @param contour.levels the number of contour levels for the two-dimensional plots (i.e., when \code{dim = 2}).
#' @param col2d vector of colors to be used in the plot.
#' @param contour.lwd line width of contour lines if \code{contour.lines = TRUE}.
#' @param contour.lines if \code{TRUE}, then black countour lines added to the graph.
#' @param kde2d.n number of grid points in each direction. Can be scalar or a length-2 integer vector.  Passes to argument \code{n} in \code{\link[MASS]{kde2d}}. 
#' @param ... passes arguments to \code{\link{plot}}.
#' @keywords plot stochastic simulations
#' 
#' @examples
#'\dontrun{
#' # First, the parameter values, as found in TSTraj
#'	model.state <- c(x = 3, y = 3)
#'	model.sigma <- 0.2
#'	model.deltat <- 0.005
#'	model.time <- 100
#'
#' # Second, write out the deterministic skeleton of the equations to be simulated, 
#' # as found in TSTraj
#'	#Example 1 from article
#'	equationx <- "1.54*x*(1.0-(x/10.14)) - (y*x*x)/(1.0 + x*x)"
#'	equationy <- "((0.476*x*x*y)/(1 + x*x)) - 0.112590*y*y"
#'
#' # Third, run it, as found in TSTraj
#'	ModelOut <- TSTraj(y0 = model.state, time = model.time, deltat = model.deltat, 
#'		x.rhs = equationx, y.rhs = equationy, sigma = model.sigma)
#' # Fourth, plot it:
#'	# in 1D
#'	TSDensity(ModelOut, dim = 1)
#'	# in 2D
#'	TSDensity(ModelOut, dim = 2, kde2d.n = 20, xlab = "")
#' }

	TSDensity <- function(mat, dim = 1, xlim = 'NULL', ylim = 'NULL', contour.levels = 15,  col2d = c("blue", "yellow", "orange", "red") , contour.lwd = 0.5, contour.lines = TRUE, kde2d.n = 100, ...){
		if (dim ==1){
			densA <- density(mat[,2] , na.rm = T)
			densB <- density(mat[,3] , na.rm = T)
			if (missing(xlim)) {xlim = c(min(c(densA$x, densB$x)) , max(c(densA$x, densB$x)))}
			if (missing(ylim)) {ylim = c(min(c(densA$y, densB$y)) , max(c(densA$y, densB$y)))}
			plot(0 , type = "n" , ylab = "Density" , xlab = "State variables" , las = 1 , xlim = xlim , ylim = ylim, ...)
			polygon(densA$x , densA$y, col = rgb(255,0,0,75,maxColorValue=255) , border = rgb(255,0,0,130,maxColorValue=255))
			polygon(densB$x , densB$y, col = rgb(0,0,255,75,maxColorValue=255) , border = rgb(0,0,255,130,maxColorValue=255))
			}
		if (dim == 2) {
#			require("MASS")	#Mass called in DESCRIPTION, Depends
			parms <- as.list(match.call())
			if (any(names(parms) == "xlab") == F) {xlab = "State variable 1"}
			if (any(names(parms) == "ylab") == F) {ylab = "State variable 2"}
			# if (any(names(parms)) == "dim") {print("yay!")}
			if(missing(xlim) | missing(ylim)){
				kern.2d <- MASS::kde2d(mat[,2], mat[,3], n = kde2d.n)
				} else {
				kern.2d <- MASS::kde2d(mat[,2] , mat[,3], lims = c(min(xlim), max(xlim), min(ylim), max(ylim)), n = kde2d.n)
				}
			x.max <- length(kern.2d$x)
			y.max <- length(kern.2d$y)
			x.range <- 1:x.max
			y.range <- 1:y.max
			contour.breaks <- seq(min(kern.2d$z) , max(kern.2d$z), length = contour.levels)
			myRmap <- colorRampPalette(col2d)(contour.levels)
			plot(0 , type = "n" , xlim = c(1 , x.max), ylim = c(1 , y.max),  xaxt = "n" , yaxt = "n", xaxs = "i", yaxs = "i", xlab ="", ylab = "", ...)
			.filled.contour(x.range , y.range , kern.2d$z , levels = contour.breaks , col = myRmap)
			contour(x.range , y.range , kern.2d$z , levels = contour.breaks , col = myRmap , add = T , drawlabels=F)
			if (contour.lines == T) {contour(x.range , y.range , kern.2d$z , levels = contour.breaks, drawlabels = F ,  add = TRUE , col = "black" , lwd = contour.lwd)}

			par(new = TRUE)
			if(missing(xlim) != missing(ylim)) {stop("Both xlim and ylim must be specified, or neither")}
			if(missing(xlim)) {xlim <- c(min(kern.2d$x), max(kern.2d$x))}
			if(missing(ylim)) {ylim <- c(min(kern.2d$y), max(kern.2d$y))}
			plot(0, type = "n" , xlim = xlim , ylim = ylim, xlab = xlab, ylab = ylab, ...)
		}
	}
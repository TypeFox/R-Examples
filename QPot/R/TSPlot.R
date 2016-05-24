#' Plot simulation of two-dimensional stochastic differential equations
#'
#' This function plots the simulation of two-dimensional stochastic differential equations from \code{\link{TSTraj}}
#' @param mat a matrix output from \code{\link{TSTraj}}.
#' @param deltat numeric value indicating the frequency of stochastic perturbation, as \eqn{\Delta t}, used in the function to recaluculate axes if applicable.
#' @param dim dimensions of the plot; \code{dim = 1} to plot a timeseries with X and Y on the ordinate axis or \code{dim = 2} to plot the trjectories in state space (i.e., X and Y respectively on the abscissa and ordinate axes).
#' @param xlim numeric vectors of length 2, giving the x coordinate range. Default \code{= 'NULL'} automatically sizes plot window.
#' @param ylim numeric vectors of length 2, giving the y coordinate range. Default \code{= 'NULL'} automatically sizes plot window.
#' @param x.lab for \code{dim = 1}, allows user to specify the axis as "time" or "steps," with steps being \eqn{time \times \Delta t}
#' @param dens if \code{dens = TRUE}, plots a horizontal one-dimensional density plot adjacent to the timerseries.
#' @param lwd line width.  Defaults to 1.
#' @param line.alpha transparency of lines from 0--255.
#' @param zero.axes if TRUE, then axes plotted at \code{X = 0} and \code{Y = 0}.
#' @param ... passes arguments to \code{\link{plot}}.
#' @keywords plot stochastic simulations
#' 
#' @examples
#' # First, the parameter values, as found in TSTraj
#'	model.state <- c(x = 3, y = 3)
#'	model.sigma <- 0.2
#'	model.deltat <- 0.05
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
#'	TSPlot(ModelOut, deltat = model.deltat, dim = 1)
#'	# in 2D
#'	TSPlot(ModelOut, deltat = model.deltat, dim = 2)

TSPlot <- function(mat, deltat, dim = 1, xlim = 'NULL', ylim = 'NULL', x.lab = "time", dens = TRUE, lwd = 2, line.alpha = 130, zero.axes = TRUE, ...) {
		if (missing(deltat) == TRUE) {stop("deltat is missing and needed to compute steps and 1D hist.  Please specify.")}
		orig.par <- par()
		global.min <- min(mat[,2:3], na.rm = T)
		global.max <- max(mat[,2:3], na.rm = T)
		x.min <- min(mat[,2], na.rm = T)
		x.max <- max(mat[,2], na.rm = T)
		y.min <- min(mat[,3], na.rm = T)
		y.max <- max(mat[,3], na.rm = T)
		if (dim == 1) {
			if(any(xlim != 'NULL')) {warning("'xlim' cannot be adjusted in this function because of calculations used in switching between displaying time and step")}
			if (table(is.infinite(mat))["FALSE"] != nrow(mat)*ncol(mat) ) { # if Inf values in the timeseries
					warning("Simulation -> Inf. Try: (i) set exact y-axis limits using the ylim argument and (ii) dens = FALSE")
				if (missing(ylim)) {
					ylim = c(global.min, global.max)
						if(dens == TRUE) {par(fig=c(0,0.775,0,1), new=FALSE , oma = rep(1,4))}
						ifelse(x.lab == "step" , x.label <- "Step" , x.label <- "Time")
					plot(mat[,1], type = "l" , ylim = ylim , las = 1 , xlab = x.label , ylab = "State variables" , col = rgb(255, 0, 0, line.alpha, NULL, 255) , lwd = lwd , xaxt = "n", ...)
					lines(mat[,1] , mat[,3] , col = rgb(0,0,255,line.alpha,NULL,255) , lwd = lwd)
					if (x.lab == "step") {axis(1)
						} else {
							par(new = TRUE)
							plot(0, type = "n" , xlim = c(0,(nrow(mat)*deltat)) , ylim = ylim , ylab = "" , xlab ="" , xaxt = "n" , yaxt = "n", ...)
							axis(1)
						}
						if(dens == TRUE) {par(fig=c(0.65,1,0,1),new=TRUE)
						densA <- density(mat[,2] , na.rm = T)
						densB <- density(mat[,3] , na.rm = T)
						y.lims <- ylim
						x.lims <- c(min(c(densA$y, densB$y)) , max(c(densA$y, densB$y)))
						plot(0 , type = "n" , xlab = "Density" , ylab = "" , las = 1 , xlim = x.lims , ylim = y.lims , yaxt = "n")
						polygon(densA$y , densA$x, col = rgb(255,0,0,75,NULL,255) , border = rgb(255,0,0,130,NULL,255))
						polygon(densB$y , densB$x, col = rgb(0,0,255,75,NULL,255) , border = rgb(0,0,255,130,NULL,255))
						axis(4 , las = 1)}
				} else {
						if(dens == TRUE) {par(fig=c(0,0.775,0,1), new=FALSE , oma = rep(1,4))}
						ifelse(x.lab == "step" , x.label <- "Step" , x.label <- "Time")
					plot(mat[,1], type = "l", las = 1 , ylim = ylim, xlab = x.label , ylab = "State variables" , col = rgb(255,0,0,line.alpha,NULL,255) , lwd = lwd , xaxt = "n", ...)
					lines(mat[,1] , mat[,3] , col = rgb(0,0,255,line.alpha,NULL,255) , lwd = lwd)
						if (x.lab == "step") {axis(1)
						} else {
							par(new = TRUE)
							plot(0, type = "n" , xlim = c(0,(nrow(mat)*deltat)), ylim = ylim , ylab = "" , xlab ="" , xaxt = "n" , yaxt = "n", ...)
							axis(1)
						}
						if(dens == T) {par(fig=c(0.65,1,0,1),new=TRUE)
						densA <- density(mat[,2] , na.rm = T)
						densB <- density(mat[,3] , na.rm = T)
						x.lims <- c(min(c(densA$y, densB$y)) , max(c(densA$y, densB$y)))
						plot(0 , type = "n" , xlab = "Density" , ylab = "" , las = 1 , xlim = x.lims , yaxt = "n", ylim = ylim)
						polygon(densA$y , densA$x, col = rgb(255,0,0,75,NULL,255) , border = rgb(255,0,0,130,NULL,255))
						polygon(densB$y , densB$x, col = rgb(0,0,255,75,NULL,255) , border = rgb(0,0,255,130,NULL,255))
						axis(4 , las = 1)}
				}
			} else { # No Inf values in timeseries
					if (missing(ylim)) {ylim <- c(global.min, global.max)}
						if(dens == TRUE) {par(fig=c(0,0.775,0,1), new = FALSE , oma = rep(1,4))}
						ifelse(x.lab == "step" , x.label <- "Step" , x.label <- "Time")
					plot(mat[,1], type = "n", ylim = ylim , las = 1 , ylab = "State variables" , xlab = x.label , lwd = lwd , xaxt = "n", ...)
					lines(mat[,1] , mat[,2] , col = rgb(255,0,0,line.alpha,NULL,255) , lwd = lwd)
					lines(mat[,1] , mat[,3] , col = rgb(0,0,255,line.alpha,NULL,255) , lwd = lwd)
					if (x.lab == "step") {axis(1)
						} else {
							par(new = TRUE)
							plot(0, type = "n" , xlim = c(0,(nrow(mat)*deltat))  , ylab = "" , xlab ="" , xaxt = "n" , yaxt = "n", ...)
							axis(1)
						}
						if(dens == TRUE) {par(fig=c(0.65,1,0,1),new=TRUE)
						densA <- density(mat[,2] , na.rm = T)
						densB <- density(mat[,3] , na.rm = T)
						x.lims <- c(min(c(densA$y, densB$y)) , max(c(densA$y, densB$y)))
						plot(0 , type = "n" , xlab = "Density" , ylab = "" , las = 1 , xlim = x.lims , ylim = ylim , yaxt = "n")
						polygon(densA$y , densA$x, col = rgb(255,0,0,75,NULL,255) , border = rgb(255,0,0,130,NULL,255))
						polygon(densB$y , densB$x, col = rgb(0,0,255,75,NULL,255) , border = rgb(0,0,255,130,NULL,255))
						axis(4 , las = 1)}
			}
		} else { # dim != 1
			if (table(is.infinite(mat))["FALSE"] != nrow(mat)*ncol(mat) ) {
				if(missing(xlim)) {stop("'Inf' detected and the plot cannot set axis limits. Try manually setting axis limits with xlim and ylim.")}
				plot(mat[,2] , mat[,3] , type = "n", las = 1 , ylim = ylim , xlim = xlim, ...)
				if (zero.axes == TRUE) {abline(v = 0 , h = 0 , col = "grey75" , lwd = 0.75 , lty = 1)}
				lines(mat[,2] , mat[,3] , col = rgb(50,50,50,line.alpha,NULL,255) , lwd = lwd)
			} else {
				if(missing(xlim)) {xlim = c(x.min, x.max)}
				if(missing(ylim)) {ylim = c(y.min, y.max)}
				plot(mat[,2] , mat[,3] , type = "n", las = 1 , ylim = ylim , xlim = xlim, ...)
				if (zero.axes == TRUE) {abline(v = 0 , h = 0 , col = "grey75" , lwd = 0.75 , lty = 1)}
				lines(mat[,2] , mat[,3] , col = rgb(50,50,50,line.alpha,NULL,255) , lwd = lwd)				
			}
		}
	par(new = F, fig = orig.par$fig)
	}

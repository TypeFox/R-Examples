#' Contour plot of quasi-potential surfaces
#'
#' This function allows users to create a contour plot of quasi-potential surfaces from \code{\link{QPGlobal}}
#' @param surface the surface to be plotted, from \code{\link{QPGlobal}}.
#' @param dens vector respectively for the number of \code{x} and \code{y} points to be plotted.
#' @param x.bound a two-element vector with the minimum and maximum x values used for computing the quasi-potential.
#' @param y.bound a two-element vector with the minimum and maximum y values used for computing the quasi-potential.
#' @param xlim numeric vectors of length 2, giving the x coordinate range.
#' @param ylim numeric vectors of length 2, giving the y coordinates range.
#' @param n.filled.contour numeric value for the nubmber of breaks in the filled contour.
#' @param n.contour.lines numeric value for the number of breaks in the contour lines.
#' @param c.parm contour line adjustment (see details).
#' @param col.contour colors to interpolate; must be a valid argument to \code{\link{col2rgb}}.
#' @param contour.lines if \code{TRUE}, then contour lines plotted over filled contour; vice versa if \code{FALSE}.
#' @param ... passes arguments to \code{\link{plot}}.
#' @details Because, in general, capturing the topological features of a surface can be subtle, we implemented a feature in \code{\link{QPContour}} to keep the filled contour region while changing the contour lines.  Specifically, \code{\link{filled.contour}} takes the range of the surface values (\eqn{\phi}), divides by the number of the specified contours (i.e., \code{n.filled.contour}), and creates a contour at each break, which happenes to be equal across the range.  But because visualizing some topology may (i) require looking between contour breaks and (ii) adding contour lines would overload the plot with lines, we use an equation to modify the distribution of contour lines.  Namely, adjusting the \code{c} argument in the \code{\link{QPContour}} function adjusts the \eqn{c} paramter in the following equation: \deqn{max_\phi \times \left(\frac{x}{n-1}\right)^c.}  This allows the user to keep the same number of contour lines (i.e., specified with \code{n.contour.lines}), but focus them toward the troughs or peaks of the surfaces. At \eqn{c=1}, the contour lines correspond to the filled.contour breaks.  If \eqn{c > 1}, then the contour lines become more concentrated towards the trough.  Similarly, if \eqn{c < 1}, then the contour lines are more focused at the peaks of the surface.  As an example, we change \eqn{c} : \cr \figure{Example3.png}.
#'
#' @examples
#' # First, System of equations
#' 	equationx <- "1.54*x*(1.0-(x/10.14)) - (y*x*x)/(1.0+x*x)"
#' 	equationy <- "((0.476*x*x*y)/(1+x*x)) - 0.112590*y*y"
#' 
#' # Second, shared parameters for each quasi-potential run
#' 	xbounds <- c(-0.5, 10.0)
#' 	ybounds <- c(-0.5, 10.0)
#' 	xstepnumber <- 150
#' 	ystepnumber <- 150
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
#' # Sixth, contour of the quasi-potential
#' 	QPContour(ex1.global, dens = c(100,100), x.bound = xbounds, 
#'		y.bound = ybounds, c.parm = 5)

QPContour <- function(surface, dens, x.bound, y.bound, xlim = 'NULL', ylim = 'NULL', n.filled.contour=25, n.contour.lines=25, c.parm=1, col.contour, contour.lines = TRUE, ...){
	x.range <- max(x.bound)-min(x.bound)
	y.range <- max(y.bound)-min(y.bound)

	row.range <- nrow(surface)-1
	col.range <- ncol(surface)-1

	row.min <- min(which(surface != 0 , arr.ind = T)[,1])
	row.max <- max(which(surface != 0 , arr.ind = T)[,1])
	col.min <- min(which(surface != 0 , arr.ind = T)[,2])
	col.max <- max(which(surface != 0 , arr.ind = T)[,2])

	x.min <- ((row.min-1)/row.range)*x.range + min(x.bound)
	x.max <- ((row.max-1)/row.range)*x.range + min(x.bound)
	y.min <- ((col.min-1)/col.range)*y.range + min(y.bound)
	y.max <- ((col.max-1)/col.range)*y.range + min(y.bound)
	
	if(missing(xlim)) {xlim <- c(x.min,x.max)}
	if(missing(ylim)) {ylim <- c(y.min,y.max)}

 	sub.x <- seq(row.min, row.max, length.out=dens[1])
	sub.y <- seq(col.min, col.max, length.out=dens[2])

	sub.x.val <- ((sub.x-1)/row.range)*x.range + min(x.bound)
	sub.y.val <- ((sub.y-1)/col.range)*y.range + min(y.bound)

	eq.sub <- surface[sub.x, sub.y]

	plot(0 , type = "n" , xlim = xlim , ylim = ylim , las = 1, ...)
	min.eq.sub <- min(eq.sub , na.rm = T)
	max.eq.sub <- max(eq.sub , na.rm = T)
	contour.breaks <- seq(min.eq.sub , max.eq.sub , length = n.filled.contour)
	eq.max <- max(surface, na.rm = T)
	line.contour.breaks <- (eq.max)*(((0:n.contour.lines)/(n.contour.lines-1)))^c.parm
	myRamp <- if(missing(col.contour)){colorRampPalette(c("#FDE725FF","#E3E418FF","#C7E020FF","#ABDC32FF","#8FD744FF","#75D054FF","#5DC963FF","#47C06FFF","#35B779FF","#28AE80FF","#20A486FF","#1F9A8AFF","#21908CFF","#24868EFF","#287C8EFF","#2C728EFF" ,"#31688EFF","#355D8DFF","#3B528BFF","#404688FF","#443A83FF","#472D7BFF","#481F71FF","#471163FF","#440154FF"))(n.filled.contour)}else{colorRampPalette(col.contour)(n.filled.contour)}
	.filled.contour(sub.x.val , sub.y.val , eq.sub , levels = contour.breaks , col = myRamp)
	if(contour.lines==TRUE){contour(sub.x.val , sub.y.val , eq.sub , levels = line.contour.breaks, drawlabels = F ,  add = TRUE , col = "black" , ...)}
	}
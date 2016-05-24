#' Vector decomposition and remainder fields
#'
#' This function calculates the vector, gradient, and remainder fields.
#' @param surface matrix output from \code{\link{QPGlobal}} or \code{\link{QPotential}}.
#' @param x.rhs a string containing the right hand side of the equation for x.
#' @param y.rhs a string containing the right hand side of the equation for y.
#' @param x.bound the x boundaries denoted at c(minimum, maximum).
#' @param y.bound the y boundaries denoted at c(minimum, maximum).
#'
#' @return returns an array of all three vector fields: the deterministic skeleton, the negative gradient of the quasi-potential, and the remainder.  The array has three dimensions with the respective lengths of xstepnumber, ystepnumber, and 6.  The six are the x- and y-values for each of the three vector fields, as x-deterministic skeleton, y-deterministic skeleton, x-negative gradient of the quasi-potential, y-negative gradient of the quasi-potential, x-remainder, and y-remainder.
#'
#' @examples
#' # First, the system of equations
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
#' 	x.bound = xbounds, y.bound = ybounds)

VecDecomAll <- function(surface,x.rhs,y.rhs,x.bound,y.bound){
	qpr <- nrow(surface)
	qpc <- ncol(surface)
	equations <- list(x.rhs,y.rhs)

	# column derivative
	r.dc <- (surface[,qpc] - surface[,(qpc-1)])
	l.dc <- (surface[,2] - surface[,1])
	int.dc <- matrix(data = NA , nrow = qpr , ncol = qpc)
	for (i in (qpc-1):2){
		int.dc[,i] <- (surface[,(i+1)]-surface[,(i-1)])/2
	}
	int.dc[,1] <- l.dc
	int.dc[,qpc] <- r.dc
	dc <- -int.dc

	# row derivative
	t.dr <- (surface[1,] - surface[2,])
	b.dr <- (surface[(qpr-1),] - surface[qpr,])
	int.dr <- matrix(data = NA , nrow = qpr , ncol = qpc)
	for (i in 2:(qpr-1)){
		int.dr[i,] <- (surface[(i-1),]-surface[(i+1),])/2
	}
	int.dr[1,] <- t.dr
	int.dr[qpr,] <- b.dr
	dr <- int.dr

	# vector and remainder fields
	x.val <- seq(min(x.bound),max(x.bound),length.out=qpc)
	y.val <- seq(min(y.bound),max(y.bound),length.out=qpr)
	n.eq <- length(equations)
	z.list <- vector(mode="list" , length = n.eq)
	for(i in 1:n.eq){
		z.list[[i]] <- 
		outer(x.val,y.val,function(x,y){eval(parse(text=equations[i]))})
		}

		# verctor field
		vx <- z.list[[1]]
		vy <- z.list[[2]]

		#remainder field
		rx <- z.list[[1]]+dr
		ry <- z.list[[2]]+dc

	array(data=c(vx,vy,dr,dc,rx,ry),dim=c(qpr,qpc,6))
}

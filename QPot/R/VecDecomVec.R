#' Vector decomposition and remainder fields
#'
#' This function calculates the vector field.
#' @param x.num.steps the number of steps between the minimum and maximum x value defined in x range.
#' @param y.num.steps the number of steps between the minimum and maximum y value defined in y range.
#' @param x.rhs a string containing the right hand side of the equation for x.
#' @param y.rhs a string containing the right hand side of the equation for y.
#' @param x.bound the x boundaries denoted at c(minimum, maximum).
#' @param y.bound the y boundaries denoted at c(minimum, maximum).
#'
#' @return returns an array of the deterministic skeleton vector field.  The array has three dimensions with the respective lengths of x.num.steps, y.num.steps, and the number of variables (always 2).  The two variables are the x-deterministic skeleton and the y-deterministic skeleton.
#'
#' @examples
#' # First, the system of equations
#' 	equationx <- "1.54*x*(1.0-(x/10.14)) - (y*x*x)/(1.0+x*x)"
#' 	equationy <- "((0.476*x*x*y)/(1+x*x)) - 0.112590*y*y"
#' 
#' # Second, shared parameters for each quasi-potential run
#' 	xbounds <- c(-0.5, 20.0)
#' 	ybounds <- c(-0.5, 20.0)
#' 	xstepnumber <- 1000
#' 	ystepnumber <- 1000
#'
#' # Third, create the deterministic skeleton vector field
#' 	VDV <- VecDecomVec(x.num.steps = xstepnumber, y.num.steps = ystepnumber, x.rhs = equationx, 
#' 	y.rhs = equationy, x.bound = xbounds, y.bound = ybounds)

VecDecomVec <- function(x.num.steps,y.num.steps,x.rhs,y.rhs,x.bound,y.bound){
	equations <- list(x.rhs,y.rhs)
	x.val <- seq(min(x.bound),max(x.bound),length.out=x.num.steps)
	y.val <- seq(min(y.bound),max(y.bound),length.out=y.num.steps)
	n.eq <- length(equations)
	z.list <- vector(mode="list" , length = n.eq)
	for(i in 1:n.eq){
		z.list[[i]] <- 
		outer(x.val,y.val,function(x,y){eval(parse(text=equations[i]))})
		}
	array(data=c(z.list[[1]],z.list[[2]]),dim=c(x.num.steps,y.num.steps,2))
}

###################################################################################################
#' Single objective test functions for SPOT
#' 
#' The following target functions are available:\cr \cr
#' spotSphereFunction:\cr
#' multi-dimensional sphere function, one global optimum: Sum[x^2] \cr \cr
#' spotSphere1Function:\cr
#' multi-dimensional sphere function, one global optimum: Sum[i^2*(x[[i]]-i)^2, {i, 1, ndim}] \cr\cr
#' spotSixHumpFunction:\cr
#' Two dimensional target function, two global optima the 6-hump camel back function\cr \cr
#' spotRosenbrockFunction:\cr
#' Two dimensional Rosenbrocks function with one global optimum, see: \url{http://en.wikipedia.org/wiki/Rosenbrock_function}.
#' Maple: 
#' f := (1-x)^2+100*(y-x^2)^2
#' plot3d(f, x = -1.5 .. 1.5, y = -.5 .. 2) \cr \cr
#' spotRosenbrockGradientFunction:\cr
#' Gradient of Rosenbrocks function.\cr \cr
#' spotRastriginFunction:\cr
#' Multi-dimensional rastrigin function, one global optimum see 
#'   \url{http://en.wikipedia.org/wiki/Rastrigin_function}\cr \cr
#' spotMexicanHatFunction:\cr
#' Two dimensional MexicanHat function, with a circular valley of global optima\cr \cr
#' spotBraninFunction:\cr
#' Two dimensional Branin function implementation, 3 global optima\cr \cr
#' spotWildFunction: \cr
#' Another test function, y=10*sin(0.3*x)*sin(1.3*x^2) + 0.00001*x^4 + 0.2*x+80
#' 
#'
#' @name Testfunctions
#' @usage spotSphereFunction(x)
#' spotSphere1Function(x)
#' spotSixHumpFunction(x)
#' spotRosenbrockFunction(x)
#' spotRosenbrockGradientFunction(x)
#' spotRastriginFunction(x)
#' spotMexicanHatFunction(x)
#' spotBraninFunction(x)
#' spotWildFunction(x)
#'			
#' @param x	vector that will be evaluated by the test-function
#'
#' @return number \code{y} \cr
#' - \code{y} is the response value of the corresponding \code{x} vector
#'
#' @seealso \code{\link{SPOT}} \code{\link{spot}} \code{\link{demo}}
#' @export spotSphereFunction spotSphere1Function  spotSixHumpFunction spotRosenbrockFunction spotRosenbrockGradientFunction spotRastriginFunction spotMexicanHatFunction spotBraninFunction spotWildFunction
#' @aliases spotSphereFunction spotSphere1Function  spotSixHumpFunction 
#' 		spotRosenbrockFunction spotRosenbrockGradientFunction spotRastriginFunction 
#'		spotMexicanHatFunction spotBraninFunction spotWildFunction
###################################################################################################
spotSphereFunction <- function (x) {
	sum(x^2)
}
spotSixHumpFunction <- function (x) {
	x1 <- x[1] 
	x2 <- x[2]	
	(4-2.1*x1^2+x1^4/3)*x1^2+x1*x2+(-4+4*x2^2)*x2^2	
}
spotRosenbrockFunction <- function (x) {
	x1 <- x[1]
	x2 <- x[2]
	((1-x1)^2)+(100*((x2-(x1^2))^2))
}
spotRastriginFunction <- function (x) {  
	sum(((x^2) - (cos(x*pi*2)*10))) + 10*length(x)
}
spotMexicanHatFunction <- function (x) {	
	x1 <- x[1] 
	x2 <- x[2] 
	
	distance <- sqrt(x1^2 + x2^2)
	if (distance == 0) # stetig ergaenzen
		{y<-1}
	else
		{y<- sin(distance) / distance}
		
	y
}
spotBraninFunction <- function (x) {	
	x1 <- x[1] 
	x2 <- x[2] 
	(x2 - 5.1/(4 * pi^2) * (x1^2) + 5/pi * x1 - 6)^2 + 10 * (1 - 1/(8 * pi)) * cos(x1) + 10
}

spotWildFunction <- function (x){
	10*sin(0.3*x)*sin(1.3*x^2) + 0.00001*x^4 + 0.2*x+80  
}

spotRosenbrockGradientFunction <- function(x) {
	x1 <- x[1]
	x2 <- x[2]
	c(-400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1), 200 * (x2 - x1 * x1))
}
spotSphere1Function <- function(x){
	sum <- 0
	for(i in 1:length(x)){
		sum <- sum+(i^2*x[[i]]^2 - i)^2
	}
	sum
}






#' Adds a line to a ternary diagram.
#' 
#' A low-level plot function which adds a line to a high-level ternary diagram.
#' 
#' This is a small utility function which helps to add a line in a ternary plot
#' from two given points in an isometric transformed space.
#' 
#' @param x Two-dimensional data set in isometric log-ratio transformed space.
#' @param \dots Additional graphical parameters passed through.
#' @return no values are returned.
#' @author Matthias Templ
#' @seealso \code{\link{ternaryDiag}}
#' @keywords aplot
#' @export
#' @examples
#' 
#' data(coffee)
#' x <- coffee[,2:4]
#' ternaryDiag(x, grid=FALSE)
#' ternaryDiagAbline(data.frame(z1=c(0.01,0.5), z2=c(0.4,0.8)), col="red")
#' 
ternaryDiagAbline <- function(x, ...){
	if(ncol(x) > 2) stop("it is assumed that the 2-dim points are provided in the transformed space.")
#	plot(x)
	k <- (x[2,2]-x[1,2])/(x[2,1]-x[1,1])
	a <- c(0,c(x[1,2]-k*(x[1,1])))
	fk <- function(x,a,k){
		y <- k*x+a[2]
		y
	}
	SEQ <- seq(-10,10,length=1000)
	x <- fk(SEQ,a,k)
#	print(x)
	x <- cbind(SEQ,x)
#	lines(x)
#	print(x)
	x <- isomLRinv(x)
	s <- rowSums(x)
	if (any(s <= 0)) 
		stop("rowSums of the input data x must be positive.")
	x <- x/s
	top <- sqrt(3)/2
	xp <- x[, 2] + x[, 3]/2
	yp <- x[, 3] * top
	lines(xp, yp, ...)	
}

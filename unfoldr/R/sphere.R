###############################################################################
# Author: M. Baaske
###############################################################################


#' Get sphere system
#'
#' Get the stored sphere system
#'
#' The sphere system is internally stored as a C data structure.
#' This function copies and converts the sphere objects to the R level.
#'
#' @param S result of \code{\link{simSphereSystem}}
#'
#' @return list of spheres of class \code{spheres}
getSphereSystem <- function(S) {
  if(length(S)==0 && class(attr(S,"eptr"))=="externalptr") {
	if(attr(S,"class")!="spheres")
		stop("Expected 'S' as object of type 'spheres'.")
	.Call(C_GetSphereSystem,attr(S,"eptr"))
  } else S
}


#' Setup sphere system
#'
#' Reinitialize sphere system after R workspace reloading
#'
#' The internally stored sphere system has to be reinitialized
#' after R workspace reloading. Calling this function is needed only in case one desires to get profile
#' sections of the sphere system again which has been previously stored as an R (list) object and afterwards
#' reloaded.
#'
#' @param S      sphere system
#' @param pl	 print level
#' @return		 \code{NULL} 
setupSphereSystem <- function(S,pl=0) {
	if(!(class(attr(S,"eptr"))=="externalptr"))
		warning("'S' has no external pointer attribute, thus we set one.")
	
	box <- attr(S,"box")	
	if(length(box)==0 || !is.list(box))
		stop("Expected argument 'box' as list type.")
	if(length(box)==1)
		box <- rep(box[1],3)
	else if(length(box)!=3)
		stop("Simulation box has wrong dimensions.")
	names(box) <- c("xrange","yrange","zrange")
	structure(.Call(C_SetupSphereSystem,as.character(substitute(S)),.GlobalEnv,
					list("lam"=0),list("box"=box, "pl"=pl)), box = box)
}

#' Simulation of sphere system
#'
#' The function simulates a Poisson sphere system of
#' intensity \code{lam} where each sphere center is uniformly
#' distributed in a box. The function returns a list of spheres with elements
#' \code{id}, \code{center} and radius \code{r}.
#'
#' Any random generating function, passed as a name, for the radii distribution is accepted as long as
#' the formal function parameter names match the actual parameter names exactly as defined in
#' the parameter list \code{theta}.
#'
#' The simulation box is of type list. The vector arguments correspond to the lower and upper points in x,y
#' and z direction. If \code{box} has only one element, i.e. \code{list(c(0,1)}, the same extent is used for the other dimensions.
#' The argument \code{pl} denotes the print level of information during simulation. Currently, only
#' \code{pl=0} for no output and \code{pl}>100 is implemented. Argument \code{cond$rdist} is of type string
#' naming the (user defined) radii random generating function.
#' Setting \code{size} equal to 'rlnorm' generates log normally distributed radii for a stationary Poisson 
#' ball system according to a general approach of exact simulation (see reference below). 
#' 
#' @param theta simulation parameters
#' @param lam   mean number of spheres per unit volume
#' @param rdist string, radii random generating function name
#' @param box 	simualtion box
#' @param pl 	print level
#'
#' @return list of class \code{spheres} if \code{pl}>100 or empty list
#'
#' @references
#'	\itemize{		
#'    \item{} {C. Lantu\eqn{\acute{\textrm{e}}}joul. Geostatistical simulation. Models and algorithms. 
#'             Springer, Berlin, 2002. Zbl 0990.86007}
#' 	 }
#' @examples
#'  theta <- list("meanlog"=-2.5,"sdlog"=0.2)
#'  S <- simSphereSystem(theta,lam=1000,rdist="rlnorm",pl=101)
simSphereSystem <- function(theta,lam,rdist,box=list(c(0,1)), pl=0) {
	theta <- list("lam"=lam,"radii"=theta)
	if(!is.numeric(lam) || !(lam>0) ) 
		stop("Expected 'lam' as non-negative numeric argument")
	if(!is.list(theta))
		stop("Expected 'theta' as list of named arguments.")
	if(!is.list(theta$radii))
		stop("Expected 'radii' as list of named arguments.")

	if(length(box)==0 || !is.list(box))
		stop("Expected 'box' as list.")
	if(length(box)==1)
		box <- rep(box[1],3)	
	if(is.null(names(box)) || !(names(box) %in% c("xrange","yrange","zrange")))
		names(box) <- c("xrange","yrange","zrange")
	
	cond <- list("rdist"=rdist,"box"=box, "pl"=pl,
			      "rho"=.GlobalEnv)

	if(cond$rdist=="const") {
		structure(.Call(C_SphereSystem, theta, cond),box = box)
	} else if(exists(cond$rdist, mode="function")) {
		fargs <- names(formals(cond$rdist))
		if(cond$rdist %in% c("rlnorm","rbeta","rgamma","runif"))
		  fargs <- fargs[-1]

		it <- match(names(theta$radii),fargs)
		if(length(it)==0 || anyNA(it))
			stop(paste("Arguments of 'radii' must match formal arguments of function ",cond$rdist,sep=""))

		structure(.Call(C_SphereSystem, theta, cond),box = box)
	} else
	  stop(paste("The ", cond$rdist, " function must be defined"))

}

#' Sphere planar section
#'
#' Intersect sphere system
#'
#' Given a sphere system obtained from \code{\link{simSphereSystem}}
#' the function returns the section radii from intersecting all spheres
#' stored in \code{S}.
#'
#' @param S list of spheres, see \code{\link{simSphereSystem}}
#' @param d distance of the intersecting xy-plane to the origin
#'
#' @return  vector of circle radii
planarSection <- function(S,d) {
	if(!is.list(S))
		stop("Expected spheres as list argument.")
	unlist(lapply(.Call(C_IntersectSphereSystem,attr(S,"eptr"),c(0,0,1),d),
		function(x) 2.0*x$r))
}

#' Binning numeric values
#'
#' Vector of count data
#'
#' This function provides basic binning (grouping) of numeric values into
#' classes defined by the breaks vector \code{bin}. The values are binned
#' according to bin[i[j]]\eqn{<}x[j]\eqn{\leq} bin[i[j]+1] for interval i=1,...,N-1
#' for \code{length(bin)=N} and value x[j].
#' If x[j] > bin[N] or x[j] < bin[1] then x[j] is not counted at all.
#'
#' @param x   	 numeric values to be binned
#' @param bin    non-decreasingly sorted breaks vector
#' @param na.rm  logical, removing missing values (including NaN) in the argument \code{x}?
#'
#' @return Vector of count data
#'
#' @examples
#' 	x <- runif(100,0,1)
#' 	bin <- seq(0,1,by=0.1)
#' 	binning1d(x,bin)
binning1d <- function(x,bin, na.rm = FALSE) {
	if (anyNA(x)) {
		if(na.rm) x <- x[which (!is.na(x))]
		else stop("Vectors contains missing values or NAs.")
	}
	if(anyNA(bin) )
	  stop("'bin' vector contains missing values or NAs.")

	if (is.unsorted(bin))
		stop("'bin' must be sorted non-decreasingly")

	.Call(C_Binning1d,x,bin)
}

#' Calculate coefficients (spheres)
#'
#' Matrix of coefficients for Wicksell's corpuscle problem
#'
#' The function calculates the matrix of coefficients of the
#' discretized integral equation for Wicksell's corpuscle problem.
#'
#' @param bin  non-decreasing vector of class limits
#' @return 	   array of coefficients
#'
#' @references
#'  Ohser, J. and Muecklich, F. Statistical analysis of microstructures in materials science J. Wiley & Sons, 2000
coefficientMatrixSpheres <- function(bin) {
	## count of bin classes
	n <- length(bin)-1
	if (n<=0)
		stop("'bin' has length zero")
	if (anyNA(bin))
		stop("Vectors contain NA values")
	if (is.unsorted(bin))
		stop("'bin' must be sorted non-decreasingly")

	if(any(bin<0))
		stop("Breaks vector must have non-negative values.")

	p <- .C(C_em_saltykov_p, as.integer(n),as.numeric(bin),
		    p=as.numeric(matrix(0,n,n)))$p
	dim(p) <- c(n,n)
	return (p)
}

#' Expectation Maximization algorithm
#'
#' Estimation of empirical sphere diameter distribution
#'
#' The function performs the EM algorithm.
#'
#' @param y   		vector of observed absolute frequencies of circle diameters
#' @param bin 		non-decreasing vector of class limits
#' @param maxIt		maximum number of iterations used
#' @return 			vector of count data of absolute frequenties of sphere diameters
#'
#' @example inst/examples/sphere.R
#'
#' @references
#' Ohser, J. and Muecklich, F. Statistical analysis of microstructures in materials science J. Wiley & Sons, 2000
em.saltykov <- function(y,bin,maxIt=32) {
	if (length(y)==0) stop("input array 'y' has length zero")
	if (anyNA(y) || anyNA(bin) )
		stop("Vectors contain NA values")
	if (is.unsorted(bin))
		stop("'bin' must be sorted non-decreasingly")

	n <- length(y)
	p <- .C(C_em_saltykov_p, as.integer(n),as.numeric(bin),
			p=as.vector(matrix(0,n,n)))$p

	theta <- y+1.0e-6
	.C(C_em_saltykov,as.integer(n), as.integer(maxIt),
			as.numeric(p),as.numeric(y),theta=as.numeric(theta))$theta
}
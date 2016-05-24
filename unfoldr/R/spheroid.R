###############################################################################
# Author: M. Baaske
###############################################################################

.section2d <- function(S) {
	idx <- ifelse(class(S)=="prolate",2,1)
	AB <- apply(sapply(S, "[[", "ab"),1,function(x) x)
	#cat("check: ",any(!(AB[,1]<AB[,2])),"\n")
	structure(list("A"=AB[,idx],
					"S"=apply(AB,1,function(x) x[2]/x[1]),
					"alpha"=unlist(lapply(S,function(x) .getAngle(x$phi)))),
			class=class(S))
}

.getAngle <- function(phi) {
	# return (abs(asin(sin(s)))) alternatively, but much slower as the following
	if(phi<=pi/2 ) phi
	else {
		if( phi <=pi) pi-phi
		else if(phi<=1.5*pi) phi%%(pi/2)
		else 2*pi-phi
	}
}


#' Maximum radius of exact simulated spheroids
#' 
#' Get maximum (random) radius 
#' 
#' In case of exact simulation the maximum radius of all randomly generated radii is returned.
#' 
#' @param S   spheroid system, simulated by perfect simulation option
#' @return    maximum radius 
getMaxRadius <- function(S) {
  .Call(C_GetMaxRadius,attr(S,"eptr"))
}


#' Updated spheroid intersections
#' 
#' Determine intersections of spheroids with bounding (simulation) box
#' 
#' For a given set of spheroid id numbers calculate if these spheroids intersect 
#' the simulation box or not.
#' 
#' @param S 	spheroid system
#' @param id    integer vector of spheroid id numbers
#' @param cond  conditioning object
#' @return 		integer vector indicating intersection=\code{1} or non intersection by \code{0}
updateIntersections <- function(S,id,cond) {
  .Call(C_UpdateIntersections,S,id,cond)	
}

#' Constuctor section profiles
#'
#' Storing structure for section profiles
#'
#' The function aggregates the necessary data for the trivariate unfolding procedure either for \code{type}
#' \code{prolate} or \code{oblate} spheroids. Argument \code{size} is a numeric matrix containing the
#' axes lengths as columns.  The \code{angle} is the orientation angle between the major axis and the vertical
#' axis direction in the section plane. If the angles range within \eqn{[0,2\pi]} these are transformed
#' to \eqn{[0,\pi/2]}. The function returns a list which consists of either longer or shorter axis \code{A}
#' of section profiles corresponding to the type of spheroids which are intended to be reconstructed,
#' the aspect ratio as the shape parameter \code{S} with values in \eqn{(0,1]}, and the orientation angle \code{alpha}.
#'
#' @param size	  matrix of major and minor axes sizes
#' @param angle   orientation angle
#' @param type    \code{prolate} or \code{oblate}, default class is \code{prolate}
#'
#' @return 		  section profiles object, either of class \code{prolate} or \code{oblate}
sectionProfiles <- function(size,angle,type=c("prolate","oblate")) {
	type <- match.arg(type)
	idx <- ifelse(type=="prolate",1,2)

	if(any(size<0))
		stop("'size' must have non-negative values.")
	if(min(angle)<0 || max(angle)>pi/2)
		stop(paste("'angle' must have values between zero and ",quote(pi/2),sep=""))

	if(max(angle)>pi/2)
		angle <- sapply(angle,.getAngle)

	structure(list("A"=size[,idx],
				   "S"=apply(size,1,function(x) x[1]/x[2]),
				   "alpha"=angle),
		   class=type)
}
#' Setup spheroid system
#'
#' Reinitialize spheroid system after R workspace reloading
#'
#' The internally stored sphere system has to be reinitialized
#' after R workspace reloading. Calling this function is needed only in case one desires to get profile
#' sections of the spheroid system again which has been previously stored as an R (list) object and afterwards
#' reloaded.
#'
#' @param S      sphere system
#' @param pl	 print level
#' @param mu	 main orientation vector
#' @return		 \code{NULL} 
setupSpheroidSystem <- function(S,mu=c(0,0,1),pl=0) {
	if(!(class(attr(S,"eptr"))=="externalptr"))
		warning(paste(substitute(S)," has no external pointer attribute, thus we set one.",sep=""))
	stype <- attr(S,"class")
	it <- pmatch(stype,c("prolate","oblate"))
	if(length(it)==0 || is.na(it))
		stop("Spheroid type 'stype' is either 'prolate' or 'oblate'.")

	box <- attr(S,"box")
	if(!is.list(box) || length(box)==0)
		stop("Expected simulation 'box' as list argument.")
	if(length(box)==1)
		box <- rep(box[1],3)
	else if(length(box)!=3)
		stop("Simulation box has wrong dimensions.")
	names(box) <- c("xrange","yrange","zrange")

	invisible(structure(
	  .Call(C_SetupSpheroidSystem,as.character(substitute(S)),.GlobalEnv,list("lam"=0),
		 	list("stype"=stype,"box"=box,"mu"=mu,"pl"=pl)),
			box = box))
}

#' Get spheroid system
#'
#' Get the internally stored spheroid system
#'
#' The spheroid system is stored as a C structure.
#' This function copies and converts the spheroids to the R level.
#'
#' @param S result of \code{\link{simSpheroidSystem}}
#' @return list of spheroids, either of class \code{prolate} or \code{oblate}
getSpheroidSystem <- function(S) {
	if(length(S)==0 && class(attr(S,"eptr"))=="externalptr") {
		if(class(S) %in% c("prolate","oblate"))
			.Call(C_GetEllipsoidSystem,attr(S,"eptr"))
		else stop("Spheroids type does not match.")
	} else S
}

#' Simulation of spheroid system
#'
#' Simulation of Poisson spheroid system
#'
#' The function simulates a Poisson spheroid system according to the supplied
#' simulation parameter \code{theta} in a predefined simulation box.
#' The argument \code{size} is of type string and denotes the major-axis length random generating 
#' function name.
#'
#' Further the function simulates either \code{stype}="prolate" or \code{stype}='oblate' spheroids.
#' For the directional orientation of the spheroid's major-axis one has the choice of a uniform
#' (\code{runifdir}), isotropic random planar (\code{rbetaiso}, see reference) or von Mises-Fisher 
#' (\code{rvMisesFisher}) distribution. The simulation box is a list containing of vector arguments 
#' which correspond to the lower and upper points in each direction. If the argument \code{box} has 
#' only one element, i.e. \code{list(c(0,1)}, the same extent is used for the other dimensions.
#' If \code{rjoint} names a joint random generating function then argument \code{size} is ignored.
#' For the purpose of exact simulation setting \code{size} equal to \code{rbinorm} declares a bivarite 
#' size-shape distribution which leads to a log normal distributed semi-major axis \code{a} and a scaled 
#' semi-minor axis length \code{c}. If \eqn{[X,Y]} follow a bivariate normal distribution with correlation parameter 
#' \eqn{\rho} then \eqn{a=exp(x)} defines the sample semi-major axis length together with the scaled semi-minor
#' axis length \eqn{c=a*s} and shape parameter set to \eqn{s=1/(1+exp(-y))}. The parameter \eqn{\rho} defines 
#' the degree of correlation between the semi-axes lengths which must be provided as part of the list of simulation 
#' parameters \code{theta}. The method of exact simulation is tailored to the above described model. For a general 
#' approach please see the given reference below. 
#' 
#' The argument \code{pl} denotes the print level of output information during simulation.
#' Currently, only \code{pl}=0 for no output and \code{pl}>100 for some additional info is implemented.
#'
#' @param theta simulation parameters 
#' @param lam   mean number of spheroids per unit volume
#' @param size  name of random generating function for size distribution
#' @param shape scalar shape parameter
#' @param orientation name of random generating function for orientation distribution
#' @param stype spheroid type
#' @param rjoint name of joint random generating function
#' @param box simulation box
#' @param mu  main orientation axis
#' @param pl  optional: print level
#' @param label optional: set a label to all generated spheroids
#' @return list of spheroids either of class \code{prolate} or \code{oblate}
#'
#' @example inst/examples/sim.R
#'
#' @references
#'	 \itemize{
#'		\item{} {Ohser, J. and Schladitz, K. 3D images of materials structures Wiley-VCH, 2009}
#'      \item{} {C. Lantu\eqn{\acute{\textrm{e}}}joul. Geostatistical simulation. Models and algorithms. 
#' 					Springer, Berlin, 2002. Zbl 0990.86007}
#' 	  }
simSpheroidSystem <- function(theta, lam, size, shape="const",orientation="rbetaiso",
								stype=c("prolate","oblate"),rjoint=NULL,box=list(c(0,1)),mu=c(0,0,1), pl=0, label="N")
{
	it <- pmatch(stype,c("prolate","oblate"))
	if(length(it)==0 || is.na(it)) stop("Spheroid type 'stype' is either 'prolate' or 'oblate'.")

	if(!is.list(theta))
		stop("Expected 'theta' as list of named  arguments.")
	if(!is.numeric(lam) || !(lam>0) ) 
		stop("Expected 'lam' as non-negative numeric argument")	
	
	if(length(box)==0 || !is.list(box))
		stop("Expected argument 'box' as list type.")
	if(length(box)==1)
	  box <- rep(box[1],3)
  	else if(length(box)!=3)
	  stop("Simulation box has wrong dimensions.")
    names(box) <- c("xrange","yrange","zrange")

	stype <- match.arg(stype)
	if(!is.null(rjoint)) {
		if(!exists(rjoint, mode="function"))
			stop("Unknown multivarirate random generating function.")
		#largs <- theta[-(which(it==1))]
		it <- match(names(theta),names(formals(rjoint)))
		if(length(it)==0 || anyNA(it))
			stop(paste("Arguments must match formal arguments of function ",rjoint,sep=""))

		# check function
		funret <- try(do.call(rjoint,theta))
		if(!is.list(funret))
		  stop("Expected list as return type in user defined function.")
		if(inherits(funret,"try-error"))
		  stop(paste("Error in user defined function ",rjoint,".",sep=""))
		if(any(!(c("a","b","u","shape","theta","phi") %in% names(funret))))
		  stop("Argument names of return value list does not match required arguments.")

		structure(.Call(C_EllipsoidSystem, 
					list("lam"=lam,"rmulti"=theta),
					list("stype"=stype,	"rdist"=rjoint,"box"=box,
							"pl"=pl,"mu"=mu,"rho"=.GlobalEnv,"label"=label)),
				   box = box)

	} else  {	
		theta <- c("lam"=lam,theta)
		it <- match(names(theta), c("lam","size","shape","orientation"))
		if(!is.list(theta) || anyNA(it))
			stop("Expected 'theta' as list of named arguments.")
		if(!is.list(theta$size) || !is.list(theta$shape) || !is.list(theta$orientation) )
			stop("Expected 'size','shape' and 'orientation' as lists of named arguments.")
		it <- pmatch(orientation,c("runifdir","rbetaiso","rvMisesFisher"))
		if(is.na(it) && !exists(orientation, mode="function"))
			stop("Unknown random generating function for orientation distribution.")

		if (missing(size))
		  stop("Argument 'size' has to be given if 'rjoint' is 'NULL'!")
		cond <- list("stype"=stype,
				"rdist"=list("size"=size, "shape"=shape,"orientation"=orientation),
				"box"=box, "pl"=pl,"mu"=mu,"rho"=.GlobalEnv,"label"=label)

		if(cond$rdist$size %in% c("const","rbinorm")) {
			structure(.Call(C_EllipsoidSystem, theta, cond), box = box)
		} else if(exists(cond$rdist$size, mode="function")) {
			fargs <- names(formals(cond$rdist$size))
			if(cond$rdist$size %in% c("rlnorm","rbeta","rgamma","runif"))
				fargs <- fargs[-1]

			it <- match(names(theta$size),fargs)
			if(length(it)==0 || anyNA(it))
				stop(paste("Arguments of 'size' must match formal arguments of function ",cond$rdist$size,sep=""))

			structure(.Call(C_EllipsoidSystem, theta, cond), box = box)		
		} else
			stop(paste("The ", cond$rdist$size, "random generating function must be defined"))

	}
}


#' Simulate spheroid system
#' 
#' Simulate a spheroid system by perfect simulation 
#'  
#' Simulate a spheroid system by perfect simulation with log normal sizes and transformed
#' shape parameter, see \code{\link{simSpheroidSystem}}. This function is intended to be a
#' condensed version of \code{\link{simSpheroidSystem}} just for ease of use of exact 
#' simulation with random planar orientation of spheroids. 
#' 
#' @param param  parameters
#' @param cond   condition object
#' 
#' @return 		 spheroid system
#' @example 	 inst/examples/sim.R
simModel3d <- function(param, cond) {	
	theta <- list("size"=list("mx"=param$mx,
							  "sdx"=param$sdx,
							  "my"=param$my,
							  "sdy"=param$sdy,
							  "rho"=param$rho),
				  "orientation"=list("kappa"=param$kappa),
				  "shape"=list())
	
	simSpheroidSystem(theta,cond$lam, size="rbinorm",
			orientation="rbetaiso",	stype=cond$stype,box=cond$box,pl=cond$pl)	
}

#' Calculate coefficients (spheroids)
#'
#' Calculate coefficients of discretized integral equation
#'
#' In order to apply the EM algorithm to the stereological
#' unfolding procedure for the joint size-shape-orientation distribution
#' one first has to calculate the coefficients of the discretized integral
#' equation. This step is the most time consuming part of unfolding
#' the parameters and therefore has been separated in its own function.
#' Further, the number of classes for size, shape and orientation do not
#' need to be equal, whereas the same class limits are used for binning spatial and planar values.
#' This might be changed in future releases.
#' Using multiple cpu cores is controlled by either setting the option \code{par.unfoldr} in the global 
#' R environment or passing the number of cores \code{nCores} directly.
#'
#'  @param   breaks  list of bin vectors
#'  @param   stype   either \code{prolate} or \code{oblate}
#'  @param   check   logical, whether to run some input checks
#'  @param   nCores  number of cores used to calculate the coefficients
#'  @return  coefficient array
#' 
#'  @example inst/examples/coeffarray.R
coefficientMatrixSpheroids <- function(breaks, stype=c("prolate","oblate"),
								check=TRUE,nCores=getOption("par.unfoldr",1))
{			
	stype <- match.arg(stype)	
	it <- match(names(breaks), c("size","angle","shape"))
	if (length(it)==0 || anyNA(it))
		stop("Expected 'breaks' as named list of: 'size','angle','shape' ")
	if (is.unsorted(breaks$size) || is.unsorted(breaks$angle) || is.unsorted(breaks$shape))
		stop("'breaks' list must contain non-decreasingly sorted classes")

	if(check) {
		if(any(breaks$size<0))
			stop("Breaks vector 'size' must have non-negative values.")
		if(min(breaks$angle)<0 || max(breaks$angle)>pi/2)
			stop(paste("Breaks vector 'angle' must have values between zero and ",pi/2,sep=""))
		if(min(breaks$shape)<0 || max(breaks$shape)>1)
			stop("Breaks vector 'shape' must have values between 0 and 1.")
	}
	
	.Call(C_CoefficientMatrixSpheroids,
			breaks$size,breaks$angle,breaks$shape,
			breaks$size,breaks$angle,breaks$shape, list(stype,nCores))

}

#' Spheroid vertical section
#'
#' Vertical section of spheroid system
#'
#' The function performs a vertical intersection of the spheroid system either
#' parallel to yz-plane defined by the normal vector \code{n=c(1,0,0)} or
#' xz-plane defined by \code{n=c(0,1,0)}.
#'
#' @param S list of spheroids, see \code{\link{simSpheroidSystem}}
#' @param d distance of intersecting plane to the origin
#' @param n normal vector
#' @return list of size, shape and angle of section profiles
verticalSection <- function(S,d,n=c(1,0,0)) {
	if(!(class(S) %in% c("prolate","oblate")))
		stop("Spheroids type does not match.")
	ss <- .Call(C_IntersectSpheroidSystem,attr(S,"eptr"),n, d)
	.section2d(ss)
}
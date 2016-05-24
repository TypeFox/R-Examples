###############################################################################
# Author: M. Baaske
###############################################################################

#cot <- function(x) tan(pi/2 - x)

#' Trivariate stereological unfolding
#'
#' Estimate the joint size-shape-orientation distribution of spheroids
#'
#' Given an array of coefficients \code{P} and the input histogram \code{F} of
#' measured planar section profiles, see \code{\link{binning3d}}, the function
#' estimates the spatial joint size-shape-orientation distribution of spheroids
#' as a trivariate histogram, \code{triHist}, by the EM algorithm.
#' If the option 'par.unfoldr' is set to a user chosen amount of cores then
#' parts of the EM iterations are done in parallel.
#'
#' @param P coefficient array
#' @param F input histogram
#' @param maxIt maximum number of EM iterations
#' @param nCores number of cpu cores used 
#' @return trivariate histogram 
#'
#' @example inst/examples/unfold.R
#'
#' @references
#' 	Bene\eqn{\check{\textrm{s}}}, V. and Rataj, J. Stochastic Geometry: Selected Topics Kluwer Academic Publishers, Boston, 2004
em.spheroids <- function(P,F,maxIt,nCores=getOption("par.unfoldr",1)) {	
	.Call(C_EMS,P,F,list("maxSteps"=maxIt,"nCores"=nCores))
}

#' Stereological unfolding
#'
#' Unfolding the (joint) distribution of planar parameters
#'
#' This is a S3 method for either trivariate stereological unfolding or
#' estimation of 3d diameter distribution of spheres (Wicksell's corpuscle problem). The function
#' aggregates all intermediate calculation steps required for the unfolding procedure given the data 
#' in the prescribed format and returning the parameters as count data in histogram form. 
#' The section profiles object \code{sp}, see \code{\link{sectionProfiles}}, is either of class
#' \code{prolate} or \code{oblate} for the reconstruction of spheroids or in case of spheres
#' the \code{sp} is simply a numeric vector of circle diameters. Here, the class of section profiles 
#' corresponds to the type of objects that will be reconstructed. The number of bin classes is set
#' by the argument \code{nclass} which is either a scalar value in case of Wicksell's corpuscle problem
#' or a vector of length three defined in the order of the number of size, angle and shape class limits.
#' Using multiple cpu cores during the calculations is controlled by either setting the option 'par.unfoldr'
#' to a user chosen amount of cores or by the argument \code{nCores} directly. 
#' The return value of the function is an object of class \code{unfold} whose arguments are as follows
#' \itemize{
#' 	\item{N_A}{ (trivariate) histogram of section profile parameters}
#'  \item{N_V}{ (trivariate) histogram of reconstructed parameters}
#'  \item{P}{ array of coefficients}     
#'  \item{breaks}{ list of class limits for binning the parameter values}
#' }
#'
#' @param sp 	  section profiles
#' @param nclass  number of classes
#' @param maxIt   maximum number of EM iterations
#' @param nCores  number of cpu cores used
#' @param ...	  optional arguments passed to \code{\link{setbreaks}}
#' @return        object of class \code{unfold}
#'
#' @seealso \code{\link{setbreaks}}, \code{\link{binning3d}}
unfold <- function(sp,nclass,maxIt=64,nCores=getOption("par.unfoldr",1),...) UseMethod("unfold", sp)
unfold.oblate <- function(sp,nclass,maxIt=64,nCores=getOption("par.unfoldr",1),...) { 
	unfold.prolate(sp,nclass,maxIt,nCores,...)
}

unfold.prolate <- function(sp,nclass,maxIt=64,nCores=getOption("par.unfoldr",1),...) {
	if(length(nclass)!=3)
	  stop("Lenght of 'dims' not equals 3.")
  	if(any(!(c("A","S","alpha") %in%  names(sp))))
		stop("Missing named arguments in 'X'")
	
	breaks <- setbreaks(nclass=nclass,maxSize=max(sp$A),...)
	N_A <- binning3d(sp$A,sp$alpha,sp$S,breaks)
	P <- coefficientMatrixSpheroids(breaks,class(sp),TRUE,nCores)

	N_V <- em.spheroids(P,N_A,maxIt,nCores)
	structure(list("N_A"=N_A,"N_V"=N_V,"P"=P,"breaks"=breaks),
			class=c("unfold",class(sp)))
}

unfold.numeric <- function(sp,nclass,maxIt=64,nCores=getOption("par.unfoldr",1),...) {
	if(anyNA(sp))
		stop("Vector of radii 'sp' has NAs.")
	if(!is.numeric(nclass) || length(nclass)!=1)
		stop("Expected numeric value 'nclass' as number of classes.")
	breaks <- seq(0,max(sp), length.out=nclass)

	# Input histogram
	#N_V <- em.saltykov(N_A,breaks,maxIt)
	y <- binning1d(sp,breaks)
		
	n <- length(y)
	p <- .C(C_em_saltykov_p, as.integer(n),as.numeric(breaks),
			p=as.vector(matrix(0,n,n)))$p
	
	theta <- y+1.0e-6
	theta <- .C(C_em_saltykov,as.integer(n), as.integer(maxIt),
				as.numeric(p),as.numeric(y),theta=as.numeric(theta))$theta

	structure(list("N_A"=y,"N_V"=theta,"P"=p,"breaks"=breaks),
			class=c("unfold",class(sp)))
}

#' Histogram data
#'
#' Count data of size, shape and orientation values.
#'
#' For each value of planar or spatial measured quantities \code{size,shape,orientation}
#' the function counts the number of observations falling into each class.
#' The \code{breaks} list is either obtained from function \code{setbreaks} or the
#' user supplies his own variant, see \code{\link{setbreaks}} for details.
#' If \code{check} is \code{TRUE} some checks on the \code{breaks} list are done.
#'
#' @param size  vector of sizes
#' @param angle vector of angles
#' @param shape	vector of aspect ratios
#' @param breaks list of bin vectors
#' @param check logical, default is \code{TRUE}
#' @param na.rm logical, if NAs are to be removed in the data vectors, default is \code{TRUE}
#' @return 3d array of count data
binning3d <- function(size,angle,shape,breaks,check=TRUE,na.rm = TRUE) {
	if(na.rm) {
		size <- na.omit(size)
		angle <- na.omit(angle)
		shape <- na.omit(shape)
	} else if(anyNA(c(size,angle,shape)))
		stop("Data vectors contain missing values or NAs.")

	if(!is.list(breaks) || length(breaks)!=3)
		stop("Expected 'breaks' argument as list of break vectors of length 3 for 'size', 'angle' and 'shape' data.")

	it <- match(names(breaks), c("size","angle","shape"))
	if (anyNA(it))
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

	# return object of class 'triHist'
	.Call(C_Binning3d,size,angle,shape,breaks$size,breaks$angle,breaks$shape)
}


#' Break vectors
#'
#' Construct class limits vectors
#'
#' The function constructs the class limits for the size, shape and orientation parameters.
#' One can either define \code{linear} class limits of 'size' as
#' \eqn{a_i=i\delta, \delta=maxSize/M } or using exponentially increasing limits: \eqn{base^i, i=1,\dots,M}.
#' The orientation classes are defined as \eqn{\theta_j=j\omega, \omega=\pi/(2N), j=1,\dots,N} in the range
#' \eqn{[0,\pi/2]}, where \eqn{M,N} are the number of size classes and the number of orientation classes, respectively.
#' Argument \code{base} must not be \code{NULL} if \code{sizeType} equals "exp".
#'
#' @param nclass 	number of classes
#' @param maxSize 	maximum of \code{size} values
#' @param base 		constant for size class construction
#' @param kap  	    constant for shape class construction
#' @param sizeType 	either \code{linear} or \code{exp}, default is \code{linear}
#' @return 			list of class limits vectors
setbreaks <- function(nclass,maxSize,base=NULL,kap=1,sizeType=c("linear","exp")) {
	if(length(nclass)==0 || any(nclass==0))
		stop("Number of bin classes 'nclass' must be greater than zero.")
	type <- match.arg(sizeType)
	
	binA <- switch(type,
		linear=seq(from=0,to=max(maxSize),by=max(maxSize)/nclass[1]),
		exp={
			if(is.null(base))
			 stop("Argument 'base' must be non NULL if 'sizeType' equals 'exp'.")
			ll <- unlist(sapply(1:nclass[1],function(i) base^i))
			if(base<1) ll <- sort(ll)
		  	ll
		})
	structure(list("size"=binA,
		"angle"=seq(from=0,to=pi/2,by=pi/(2*nclass[2])),
		"shape"=unlist(lapply(0:nclass[3],function(i) (i/nclass[3])^kap))))
}

#' Original parameters
#'
#' Extract the original 3d size-shape-orientation parameters
#'
#' The function simply extracts the parameters for the original 3d spheroid system
#' either for oblate or prolate spheroids and returns a list which consists of sizes
#' \code{a}, the orientation angles \code{Theta} and the shape parameters \code{s}.
#'
#' @param S list of spheroids
#' @return  list
parameters3d <- function(S) {
	idx <- if(class(S)=="prolate") c(1,2) else c(2,1)
	list("a"=unlist(lapply(S,function(x) x$ab[1])),
		 "Theta"=unlist(lapply(S,function(x) .getAngle(x$angles[1]))),
		 "s"=unlist(lapply(S,function(x) x$ab[idx[1]]/x$ab[idx[2]])))
}

#' Estimated spatial histogram data
#'
#' Get histogram data from estimated joint distribution
#'
#' Given the estimated joint distribution in histogram form the function
#' returns a list of \code{breaks} replicates which equals the number of
#' estimated count data for each class.
#'
#' @param H trivariate histogram
#' @param breaks breaks as obtained from \code{\link{setbreaks}}
#' @return list of size, angle and shape parameters, see \code{\link{parameters3d}}
parameterEstimates <- function(H,breaks) {
	list("a"=unlist(lapply(1:(length(breaks$size)-1),
							function(i) rep(breaks$size[i],sum(H[i,,])))),
		 "Theta"=unlist(lapply(1:(length(breaks$angle)-1),
							function(i) rep(breaks$angle[i],sum(H[,i,])))),
		 "s"=unlist(lapply(1:(length(breaks$shape)-1),
							function(i) rep(breaks$shape[i],sum(H[,,i])))))
}

#' Trivariate histogram
#'
#' Plot trivariate histogram of joint size-shape-orientation distribution
#'
#' The (estimated spatial) joint size-shape-orientation distribution is plotted
#' in a box with corresponding axes shown. The axes intersect in the first class number.
#' The ball volumes visualize the relative frequencies of count data for each class which 
#' can be scaled by the user for non-overlapping spheres.
#' Balls within the same size class have the same color.
#'
#'  @param A 		3d array of count data
#'  @param main 	main title
#'  @param scale 	factor to scale the spheres
#'  @param col 		vector of color values repeatedly used
#'  @param ... 		optional graphic arguments passed to function \code{\link[rgl]{spheres3d}}
#' 
#'  @return 		\code{NULL}
trivarHist <- function(A, main = paste("Trivariate Histogram"),scale = 0.5,col, ...) {	
	N <- sum(A)
	pos <- do.call(rbind,lapply(seq(1:dim(A)[1]),
					function(i) {
						X <- which(A[i,,]!=0,arr.ind=TRUE)
						cbind(X,rep(i,nrow(X)))
					}))
	## scaling of balls
	sz <- apply(pos, 1,function(x) scale*A[x[3],x[1],x[2]]/max(A))

	# coloring the balls
	xt <- as.vector(table(pos[,3]))
	cols2 <- rep(col,length.out=dim(A)[1])
	ncols <- lapply(1:dim(A)[1], function(i) rep(cols2[i], length.out=xt[i]))
	
	#plot spheres	 
	rgl::spheres3d(pos,radius=sz,col=unlist(ncols),...)
	rgl::axes3d(c('x','y','z'), pos=c(1,1,1), tick=FALSE)
	rgl::title3d(main,'',"orientation","shape","size")
}
#************************************************************************
#           Boundary Kernel Density class
#************************************************************************

setClass(
  Class = "BoundaryKernel",
  representation = representation(
	mu = "numeric",
    "VIRTUAL"),
  contains = "KernelDensity"
  )

setValidity(
  Class = "BoundaryKernel",
  method = function(object) {
    if (object@mu != 0 & object@mu != 1 & object@mu != 2 & object@mu != 3){
      stop("Mu can only take values 0,1,2 or 3")
	}else if (object@b < 0 | object@b > 1){
      stop("The parameter b must be 0<=b<=1")
	}else if (object@b > 0.5){
		cat(c("WARNING: b parameter is too large (b =", object@b, "). If b > 0.5 left and right boundary overlap and \n it may result", 
			"in extrange discontinuities in the density function\n"))
    }else{}	
    return(TRUE)
  }
  )

setGeneric (
  name = "getmu",
  def  = function(.Object){standardGeneric("getmu")}
  )

setMethod(
  f = "getmu",
  signature = "BoundaryKernel", 
  definition = function(.Object) {
    return(.Object@mu)
  }
  )


setGeneric (
  name = "leftBoundaryKernelFunction",
  def  = function(.Object,q,u){standardGeneric("leftBoundaryKernelFunction")}
  )

setGeneric (
  name = "interiorKernelFunction",
  def  = function(.Object,u){standardGeneric("interiorKernelFunction")}
  )

setGeneric (
  name = "rightBoundaryKernelFunction",
  def  = function(.Object,q,u){standardGeneric("rightBoundaryKernelFunction")}
  )

setGeneric (
  name = "singlePoint.density",
  def  = function(.Object,x){standardGeneric("singlePoint.density")}
  )

setMethod(
	f = "singlePoint.density",
	signature = "BoundaryKernel",
	definition = function(.Object,x){
		
		if(length(x) != 1){
			stop("This method can only evaluate the density in a single data point")
		}else{}
		
		indices <- abs(.Object@dataPoints - x)/.Object@b <= 1
  		subsamples <- .Object@dataPoints[indices]
		
		u <- (x-subsamples)/.Object@b # evaluation points where the kernel functions are evaluated
		if (x<.Object@b){ # Use left boundary kernel function
			q <- x/.Object@b
			sum(leftBoundaryKernelFunction(.Object,q,u)) / (.Object@b*length(.Object@dataPoints))
		} else if ((x>=.Object@b & x<=(1-.Object@b))){ # Use interior kernel function
			sum(interiorKernelFunction(.Object,u)) / (.Object@b*length(.Object@dataPoints))		
		} else if (x > (1-.Object@b)){ # Use right boundary kernel function
			q <- (1-x)/.Object@b			
			sum(rightBoundaryKernelFunction(.Object,q,u)) / (.Object@b*length(.Object@dataPoints))
		} else {}
		
})


setMethod(
  f = "density", 
  signature = "BoundaryKernel",
  definition = function(x,values,scaled = FALSE) {    
    
    .Object <- x
    x <- values
	
	isMatrix.x <- is.matrix(x)
    #dims = [nrows,ncols]
    dims <- dim(x) 
    
    if(!scaled){
      # scale the data to the 0-1 interval
      x <- getScaledPoints(.Object,x)    
    }
	
	# if any value in x is lower than 0 or grater than 1 its density is 0
	numDataPoints <- length(x)
	index.nozero <- which(x>=0 & x <=1)
	
	x <- x[index.nozero]
	if(length(x) == 0){ # all elements in x were out of bound
		return(rep(0,numDataPoints - length(index.nozero)))
	}
        
    # x is considered as a vector even if it is a matrix(elements taken by columns) 
	x.indices <- rep(0,times = length(x))
    x.densities <- numeric(0)

    if(length(.Object@densityCache) == length(.Object@dataPointsCache)){      
      # if there are density values calculated in the cache, first we look
      # at the cache to check whether some of the values in x have been already calculated      
      x.indices <- match(x, .Object@dataPointsCache, nomatch=0)
      
      if(any(x.indices > 0)){
        # the density of some of the points are already calculated in the cache
        x.densities[x.indices != 0] <- .Object@densityCache[x.indices[x.indices!=0]]
      }else{}
      
    }else{}

    # the data poins whose densities are not calculated in the cache
    x.new <- x[x.indices == 0]
    x.new.length <- length(x.new)
    if(x.new.length > 0){
		# There are densities to be calculated
		densities <- sapply(x.new, FUN = function(x,.Object){
				singlePoint.density(.Object,x)}, .Object= .Object) 
		x.densities[x.indices == 0] <- densities
	}else{}

	 # include the density (density=0) of the out-of-bound x points in the final result
	aux.density <- numeric(numDataPoints)
	aux.density[index.nozero] <- x.densities
	x.densities <- aux.density
    
	#if x is a matrix, we storte the densities as a matrix object
	if(isMatrix.x){
		dim(x.densities) <- dims
	}else{}

	#if data are in another scale (not in the [0,1] interval) we should
	# normalize the density by dividing it by the length
	# of the interval so that the density integrates to 1
    domain.length <- .Object@upper.limit - .Object@lower.limit
    if(!scaled){
      x.densities <- x.densities/domain.length
    }
    return(x.densities)
  })

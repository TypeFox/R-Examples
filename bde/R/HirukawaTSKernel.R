#************************************************************************
# TS Kernel estimator as described in Hirukawas's 2010 paper:
#
#	@article{hirukawa2010,
#	title = {Nonparametric multiplicative bias correction for kernel-type density estimation on the unit interval},
#	author = {Hirukawa, M.},
#	journal = {Computational Statistics & Data Analysis},
#	year = {2010},
#	volume = {54},
#	number = {2}
#	pages = {473--495}
#	}
#
#************************************************************************

setClass(
  Class = "HirukawaTSKernel",
  representation = representation(
	c = "numeric"),
  prototype = prototype(modified = FALSE),
  contains = "Chen99Kernel"
  )

setValidity(
  Class = "HirukawaTSKernel",
  method = function(object) {
    if (object@c <= 0 & object@c >= 1){
		stop("HirukawaTSKernel, parameter c must be 0<c<1")
	}else if (object@b/object@c > 1){
      stop("In HirukawaTSKernel density calculations, a Chen99Kernel with parameter b=b/c is used. The provided values for b and c make\n",
		"b/c >1 and it must be within the (0,1) interval\n")
    }else if (object@modified & object@b/object@c > 0.25){
		cat(c("In HirukawaTSKernel density calculations, a Chen99Kernel with parameter b=b/c is used. The provided values for b and c make\n",
			"b/c > 0.25 and therefore left and right boundaries overlap which may result in extrange discontinuities in the density function\n"))
    }else{}
    return(TRUE)
  }
  )


setGeneric (
  name = "getc",
  def  = function(.Object){standardGeneric("getc")}
)

setMethod(
  f = "getc", 
  signature = "HirukawaTSKernel",
  definition = function(.Object) {
    .Object@modified
  })

setMethod(
  f = "density", 
  signature = "HirukawaTSKernel",
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
		kernel1 <- new(Class="Chen99Kernel",dataPoints=.Object@dataPoints, b=.Object@b, dataPointsCache=.Object@dataPointsCache,
						modified=.Object@modified)
		kernel2 <- new(Class="Chen99Kernel",dataPoints=.Object@dataPoints, b=.Object@b/.Object@c, dataPointsCache=.Object@dataPointsCache,
						modified=.Object@modified)
				
# 		x.densities[x.indices == 0] <- density(kernel1,x.new,scaled=TRUE)^(1/(1-.Object@c)) * density(kernel2,x.new,scaled=TRUE)^(-.Object@c/(1-.Object@c))
    aux1 <- density(kernel1,x.new,scaled=TRUE)^(1/(1-.Object@c))
    aux2 <- density(kernel2,x.new,scaled=TRUE)^(-.Object@c/(1-.Object@c))
    aux3 <- aux1 * aux2		
    # when the estimated density tends to 0 when estimating with the Chen99Kernel (kernel1 or kernel2), it may be possible 
    # that aux1 -> 0 and aux2 -> inf which results in a NaN. In those cases, the density with HirukawaTSKernel should 
    # tends alsto to 0
    aux3[aux1 == 0] <- 0
    x.densities[x.indices == 0] <- aux3
	}
	else{}
	
	# include the density (density=0) of the out-of-bound x points in the final result
	aux.density <- numeric(numDataPoints)
	aux.density[index.nozero] <- x.densities
	x.densities <- aux.density
    
	#if x is a matrix, we storte the densities as a matrix object
	if(isMatrix.x){
		dim(x.densities) <- dims
	}
		
	#if data are in another scale (not in the [0,1] interval) we should
	# normalize the density by dividing it by the length
	# of the interval so that the density integrates to 1
    domain.length <- .Object@upper.limit - .Object@lower.limit
    if(!scaled){
      x.densities <- x.densities/domain.length
    }    
    return(x.densities)
  }
  )




#####################################
## Constructor functions for users ##
#####################################

hirukawaTSKernel <- function(dataPoints, c, b=length(dataPoints)^(-2/5), dataPointsCache=NULL, modified = FALSE, lower.limit=0,upper.limit=1){
  #cat("~~~~~~ HirukawaTSKernel: constructor ~~~~~~\n")
  
  dataPoints.scaled <- dataPoints
  dataPointsCache.scaled <- dataPointsCache
  if(is.null(dataPointsCache)){
	dataPointsCache.scaled <- seq(0,1,0.01)
  }
  
  if(lower.limit!=0 || upper.limit!=1){
    dataPoints.scaled <- (dataPoints-lower.limit)/(upper.limit-lower.limit)
    if(!is.null(dataPointsCache)){
      dataPointsCache.scaled <- (dataPointsCache-lower.limit)/(upper.limit-lower.limit)
    }    
  }
  
  kernel <- new(Class="HirukawaTSKernel",dataPoints = dataPoints.scaled, b = b, dataPointsCache = dataPointsCache.scaled, modified = modified, c = c, lower.limit=lower.limit,upper.limit=upper.limit)
  setDensityCache(kernel, densityFunction=NULL)
  return(kernel)
}

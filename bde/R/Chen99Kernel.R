#************************************************************************
# Kernel estimator as described in Chen's 99 paper:
#
#	@article{Chen1999,
#	title = {Beta kernel estimators for density functions},
#	author = {Chen, Song Xi},
#	journal = {Computational Statistics \& Data Analysis},
#	year = {1999},
#	volume = {31},
#	pages = {131--145}
#	}
#
#************************************************************************

setClass(
  Class = "Chen99Kernel",
  representation = representation(modified = "logical"),
  prototype = prototype(modified = FALSE),
  contains = "KernelDensity"
  )

setValidity(
  Class = "Chen99Kernel",
  method = function(object) {
    if (length(object@dataPoints) == 0){
      stop("A data set with at least one point is needed")
    }else if (any(object@dataPoints < 0) || any(object@dataPoints > 1)){
      stop("Data points outside the bounds")
    }else if (object@modified & object@b > 0.25){
		cat(c("WARNING: b parameter is too large (b =", object@b, "). If modified=TRUE and b > 0.25 left and right boundary overlap and \n it may result", 
			"in extrange discontinuities in the density function\n"))
	}else{}
    return(TRUE)
  }
  )

setGeneric (
  name = "getmodified",
  def  = function(.Object){standardGeneric("getmodified")}
)

setMethod(
  f = "getmodified", 
  signature = "Chen99Kernel",
  definition = function(.Object) {
    .Object@modified
  })


setMethod(
  f = "density", 
  signature = "Chen99Kernel",
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
	x.indices <- numeric(length(x))
    x.densities <- numeric(length(x))

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
		shape1 <- numeric(x.new.length)
		shape2 <- numeric(x.new.length)
		
		if(.Object@modified){			
			# data points in x.new so that they are < 2*b has the following shape paremeters
			index.cond1 <- which(x.new < 2*.Object@b)
			shape1[index.cond1] <- rho(.Object,x.new[index.cond1])
			shape2[index.cond1] <- (1-x.new[index.cond1]) / .Object@b
			# data points in x.new so that they are ( >= 2*b & <= 1-2*b) has the following shape paremeters
			index.cond2 <- which((x.new >= 2*.Object@b) & (x.new <= 1-2*.Object@b))
			shape1[index.cond2] <- x.new[index.cond2] / .Object@b
			shape2[index.cond2] <- (1-x.new[index.cond2]) / .Object@b
			# the rest of data points in x.new
			index.rest <- which(x.new > 1-2*.Object@b)
			shape1[index.rest] <- (x.new[index.rest] / .Object@b)
			shape2[index.rest] <- rho(.Object,1-x.new[index.rest])				
		}else{
			shape1 <- (x.new / .Object@b) + 1
			shape2 <- ((1-x.new) / .Object@b) + 1			
		}
		
		n <- length(.Object@dataPoints)
		x.aux <- matrix(rep(.Object@dataPoints,x.new.length), nrow=x.new.length, byrow=TRUE)
		aux <- cbind(x.aux,shape1,shape2)
		densities <- apply(aux,MARGIN=1,function(x)               
			dbeta(x[1:n], x[n+1], x[n+2]))
		
		if(!is.matrix(densities)){
 			densities <- matrix(densities, nrow=n)
		}else{}	
		
		x.densities[x.indices == 0] <- colSums(densities) / n	
	}else {}

    # include the density (density=0) of the out-of-bound x points in the final result
	aux.density <- numeric(numDataPoints)
	aux.density[index.nozero] <- x.densities
	x.densities <- aux.density
	
	#if x is a matrix, we store the densities as a matrix object
	if(isMatrix.x){
		dim(x.densities) <- dims
	}
	
    #if data are in another scale (not in the [0,1] interval) we should
	# normalize the density by dividing it by the length
	# of the interval so that the density integrates to 1
    domain.length <- .Object@upper.limit - .Object@lower.limit
    if(!scaled){
      aux.density <- aux.density/domain.length
    }
    return(aux.density)
  }
  )


#######################################################
## Auxiliar function used in Chen99 modified kernels ##
#######################################################

setGeneric (
  name = "rho",
  def  = function(.Object,x){standardGeneric("rho")}
  )

setMethod(
  f = "rho", 
  signature = "Chen99Kernel",
  definition = function(.Object, x){
	2*.Object@b^2 + 2.5 - sqrt(4*.Object@b^4 + 6*.Object@b^2 + 2.25 - x^2 -x/.Object@b)
})


#####################################
## Constructor functions for users ##
#####################################

chen99Kernel <- function(dataPoints, b=length(dataPoints)^(-2/5), dataPointsCache=NULL, modified = FALSE, lower.limit=0,upper.limit=1){
  #cat("~~~~~~ Chen99Kernel: constructor ~~~~~~\n") 
  
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
  
  kernel <- new(Class="Chen99Kernel",dataPoints = dataPoints.scaled, b = b, dataPointsCache = dataPointsCache.scaled, modified = modified, lower.limit=lower.limit,upper.limit=upper.limit)
  setDensityCache(kernel, densityFunction=NULL)
  return(kernel)
}

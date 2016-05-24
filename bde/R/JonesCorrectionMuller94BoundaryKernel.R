#************************************************************************
#  Jone's nonnegative boundary correction for Muller94's Boundary Kernel
#  Density estimator
#
# @article{Jones1996,
#  title = {A simple nonnegative boundary correction method for kernel density estimation},
#  author = {Jones, M.C. and Foster, P.J.},
#	journal = {Statistica Sinica},
#	year = {1996},
#	volume = {6},
#	pages = {1005-1013}
#	}
#************************************************************************

setClass(
  Class = "JonesCorrectionMuller94BoundaryKernel",
  representation = representation(
    normalizedKernel = "NormalizedBoundaryKernel"),
  contains = "Muller94BoundaryKernel"
)


setMethod(
  f = "density", 
  signature = "JonesCorrectionMuller94BoundaryKernel",
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
      # the density value of the Muller91BoundaryKernel #borrar la densityCache antes de llamar a callNextMethod para por si acaso
      f.muller <- callNextMethod(.Object,x.new,scaled=TRUE)
      f.localRenormalized <- density(.Object@normalizedKernel,x.new,scaled=TRUE)
      
      densities <- f.localRenormalized*exp((f.muller/f.localRenormalized)-1)
      # The density when f.localRenormalized = 0 is also 0
      densities[which(f.localRenormalized==0)] =  0
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



#####################################
## Constructor functions for users ##
#####################################

jonesCorrectionMuller94BoundaryKernel <- function(dataPoints, mu=1, b=length(dataPoints)^(-2/5), dataPointsCache=NULL, lower.limit=0,upper.limit=1){
  #cat("~~~~~~ Jones Correction for Muller94 Boundary Kernel: constructor ~~~~~~\n")  
  
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
  
  kernel.localRenormalized <- normalizedBoundaryKernel(dataPoints=dataPoints, mu=mu, b=b, dataPointsCache=dataPointsCache, lower.limit=lower.limit,upper.limit=upper.limit)
  kernel <- new(Class="JonesCorrectionMuller94BoundaryKernel",dataPoints = dataPoints.scaled, b = b, 
                dataPointsCache = dataPointsCache.scaled, mu = mu, normalizedKernel=kernel.localRenormalized,
                lower.limit=lower.limit,upper.limit=upper.limit)
  setDensityCache(kernel, densityFunction=NULL)
  return(kernel)
}
#************************************************************************
# MacroBeta normalized JLN Kernel estimator as described in Hirukawas's 2010 paper:
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
  Class = "MacroBetaHirukawaJLNKernel",
  representation = representation(	
	normalizationConst = "numeric"
  ),  
  contains = "HirukawaJLNKernel"
  )



setGeneric (
  name = "setNormalizationConstant",
  def  = function(.Object){standardGeneric("setNormalizationConstant")}
  )

setMethod(
	f = "setNormalizationConstant",
	signature = "MacroBetaHirukawaJLNKernel",
	definition = function(.Object){
		
		f <- function(x,kernel){
        density(as(kernel,"HirukawaJLNKernel"),x,scaled=TRUE)
     	}
		
		normalizationConst <- integrate(f,lower=0,upper=1,kernel=.Object)$value
		
		# Save the normalization constants in the kernel and export the results to the global object
		# obtain the global name of the variable to modify
		objectGlobalName <- deparse(substitute(.Object))    
		# modify the variable
        .Object@normalizationConst <- normalizationConst
		# assign the local variable to the global variable  
		assign(objectGlobalName,.Object,envir=parent.frame())  	
}
)


setMethod(
  f = "density", 
  signature = "MacroBetaHirukawaJLNKernel",
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
        
    # if the normalization constants have not been calculated yet, we calculate them now
	if(length(.Object@normalizationConst) == 0){
		# obtain the global name of the variable to modify
		objectGlobalName <- deparse(substitute(.Object))    
		# modify the variable
        setNormalizationConstant(.Object)
		# assign the local variable to the global variable  
		assign(objectGlobalName,.Object,envir=parent.frame())  			
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
		#just in case we remove the density cache since it is calculated for a MacroBetaChen99Kernle
		forceDensityCacheTo(.Object,numeric(0))
		x.densities[x.indices == 0] <- callNextMethod(.Object, x.new, scaled=TRUE)/.Object@normalizationConst
	}else{}

	# include the density (density=0) of the out-of-bound x points in the final result
	aux.density <- numeric(numDataPoints)
	aux.density[index.nozero] <- x.densities
	x.densities <- aux.density
	
	#if x is a matrix, we store the densities as a matrix object
	if(isMatrix.x){
		dim(x.densities) <- dims
	}
	
	# if data are in another scale (not in the [0,1] interval) we should
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

macroBetaHirukawaJLNKernel <- function(dataPoints, b=length(dataPoints)^(-2/5), dataPointsCache=NULL, modified = FALSE, lower.limit=0,upper.limit=1){
  #cat("~~~~~~ MacroBetaHirukawaJLNKernel: constructor ~~~~~~\n")
  
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
  
  kernel <- new(Class="MacroBetaHirukawaJLNKernel",dataPoints = dataPoints.scaled, b = b, dataPointsCache = dataPointsCache.scaled, modified = modified, lower.limit=lower.limit,upper.limit=upper.limit)
  
  setNormalizationConstant(kernel)
  setDensityCache(kernel, densityFunction=NULL)  
  return(kernel)
}

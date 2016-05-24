#************************************************************************
# JLN Kernel estimator as described in Hirukawas's 2010 paper:
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
  Class = "HirukawaJLNKernel",
  representation = representation(
	modified = "logical"),
  prototype = prototype(modified = FALSE),
  contains = "Chen99Kernel"
  )


setMethod(
  f = "density", 
  signature = "HirukawaJLNKernel",
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
		#kernel1 <- new(Class="Chen99Kernel",dataPoints=.Object@dataPoints, b=.Object@b, dataPointsCache=.Object@dataPointsCache,
		#				modified=.Object@modified)
		kernel1 <- chen99Kernel(dataPoints=.Object@dataPoints, b=.Object@b, dataPointsCache=.Object@dataPointsCache,
						modified=.Object@modified)
				
				
		#Calculate the normalized term
		# The parameters for the kernel function 
		shape1 <- numeric(0)
		shape2 <- numeric(0)
		
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
		
		# Calculate the numerator in the normalized term
		n <- length(.Object@dataPoints)
		x.aux <- matrix(rep(.Object@dataPoints,x.new.length), nrow=x.new.length, byrow=TRUE)
		aux <- cbind(x.aux,shape1,shape2)
		numerator <- apply(aux,MARGIN=1,function(x)               
			dbeta(x[1:n], x[n+1], x[n+2]))
		
		if(!is.matrix(numerator)){
 			numerator <- matrix(numerator, nrow=n)
		}else{}	
	
		# Calculate the denominator in the normalized term
		den.aux <- density(kernel1,.Object@dataPoints,scaled=TRUE)
		denominator <- matrix(rep(den.aux, x.new.length), nrow=n, byrow=FALSE)
		
		
		normalized.term <- (1/n) * colSums(numerator/denominator)
	
	
		
		
		x.densities[x.indices == 0] <- density(kernel1,x.new) * normalized.term
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

hirukawaJLNKernel <- function(dataPoints, b=length(dataPoints)^(-2/5), dataPointsCache=NULL, modified = FALSE, lower.limit=0,upper.limit=1){
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
  
  kernel <- new(Class="HirukawaJLNKernel",dataPoints = dataPoints.scaled, b = b, dataPointsCache = dataPointsCache.scaled,
                modified = modified, lower.limit=lower.limit, upper.limit=upper.limit)
  setDensityCache(kernel, densityFunction=NULL)
  return(kernel)
}

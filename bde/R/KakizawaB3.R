#************************************************************************
# Kakizawa Bernstein Polynomial approximation as described in 
# Kakizawa's 2004 paper:
#
#	@article{Kakizawa04,
#	title = {Bernsteins polynomial probability density estimation},
#	author = {Kakizawa, Y.},
#	journal = {Journal of Nonparametric Statistics},
#	year = {2004},
#	volume = {16},
#   number = {5},
#	pages = {709-729}
#	}
#
#************************************************************************

setClass(
  Class = "KakizawaB3",
  representation = representation(	
	densityEstimator = "BoundedDensity"),
  contains = "BernsteinPolynomials"
  )

setValidity(
  Class = "KakizawaB3",
  method = function(object) {
    if (length(object@dataPoints) == 0){
      stop("A data set with at least one point is needed")
    }else if (any(object@dataPoints < 0) || any(object@dataPoints > 1)){
      stop("Data points outside the bounds [lower.limit,upper.limit]")
	}else if (object@m < 0){
      stop("The order for Bernstein estimator (m) must be > 0")
    }else{}
    return(TRUE)
  }
  )

setMethod(
  f = "density", 
  signature = "KakizawaB3",
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
	  j <- 0:(.Object@m-1)
	  # aux.pointsForEstimation is a matrix with all the points (j+x)/m where the density is evaluated by the empirical distribution
  	  # (density estimator in this case) for the polinomial approach
	  aux.pointsForEstimation <- (matrix(rep(j, x.new.length), ncol = x.new.length) + matrix(rep(x.new,.Object@m),ncol = x.new.length, byrow=TRUE)) / .Object@m	
	  # The estimated density value at point (j+x)/m
	  f.star <- density(.Object@densityEstimator, aux.pointsForEstimation,scaled=TRUE)
	  # The binomial probability (P_{k,m-a}) 
	  P <- sapply(x.new,FUN = function(x,j,m){
			#dbeta(x,shape1,shape2)
			choose(m,j)*x^j*(1-x)^(m-j)			
		},
		j = j, m = .Object@m-1)
	
	densities <- colSums(f.star * P)
    x.densities[x.indices == 0] <- densities      
    }else{}
    
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


#
#####################################
## Constructor functions for users ##
#####################################

kakizawaB3 <- function(dataPoints, estimator=NULL,m=round(length(dataPoints)^(2/5)),dataPointsCache=NULL, lower.limit=0,upper.limit=1){
  #cat("~~~~~~ KakizawaB3: constructor ~~~~~~\n")
  
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
  
  if(is.null(estimator)){
    estimator <- muller94BoundaryKernel(dataPoints=dataPoints.scaled)
  }
  
  polinomialModel <- new(Class="KakizawaB3",dataPoints = dataPoints.scaled, m = m, dataPointsCache = dataPointsCache.scaled,
                         densityEstimator = estimator, lower.limit=lower.limit,upper.limit=upper.limit)   
  setDensityCache(polinomialModel, densityFunction=NULL)
  return(polinomialModel)
}

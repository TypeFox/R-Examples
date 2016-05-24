#************************************************************************
# Vitale 75 Bernstein Polynomial approximation as described in 
# Leblanc's 09 paper:
#
#	@article{Leblanc2010,
#	title = {A bias-reduced approach to density estimation using 
#            Bernstein polynomials},
#	author = {Leblanc, Alexandre},
#	journal = {Journal of Nonparametric Statistics},
#	year = {2010},
#	volume = {22},
#   number = {4},
#	pages = {459--475}
#	}
#
#   @article{Vitale1975,
#	title = {A Bernstein polynomial approach to density function
#            estimation},
#	author = {Vitale, R. A.},
#	journal = {Statistical Inference and Related Topics},
#	year = {1975},
#	volume = {2},
#	pages = {87--99}
#	}
#
#************************************************************************

setClass(
  Class = "Vitale",
  representation = representation(),
  contains = "BernsteinPolynomials"
  )

setValidity(
  Class = "Vitale",
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
  signature = "Vitale",
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
	  aux <- sapply(0:(.Object@m - 1), FUN = 
		function(k,m,x){
			n <- length(.Object@dataPoints)
			#Fn <- sum((.Object@dataPoints > k/m) & (.Object@dataPoints <= (k+1)/m)) / n
			Fn <- (sum(.Object@dataPoints <= (k+1)/m) - sum(.Object@dataPoints <= k/m)) / n
			return(Fn * dbeta(x,k+1,m-k))
		}, m = .Object@m, x = x.new)

	  # if x.new contains more than one data points, aux is a matrix where the calculated values for each 
	  # data point in x.new are stored by rows. However, if x.new is a single value, aux cointains a vector
 	  # with the calculated values the data point in x.new
	  if(is.matrix(aux)){
		densities <- rowSums(aux)  
	  }else{
		densities <- sum(aux)
	  }
		     
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
	
	##if data are in another scale (not in the [0,1] interval) we should
	## normalize the density by dividing it by the length
	## of the interval so that the density integrates to 1
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

vitale <- function(dataPoints, m=round(length(dataPoints)^(2/5)), dataPointsCache=NULL, lower.limit=0,upper.limit=1){
  #cat("~~~~~~ Vitale: constructor ~~~~~~\n")  
  
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
  
  polinomialModel <- new(Class="Vitale",dataPoints = dataPoints.scaled, m = m, dataPointsCache = dataPointsCache.scaled,
                         lower.limit=lower.limit,upper.limit=upper.limit)
  setDensityCache(polinomialModel, densityFunction=NULL)
  return(polinomialModel)
}

#************************************************************************
# Bias reduced version of Vitale 75 Bernstein Polynomial approximation 
# proposed in Leblanc's 09 paper:
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
  Class = "BrVitale",
  representation = representation(
	M = "numeric"),
  contains = "Vitale"
  )

setValidity(
  Class = "BrVitale",
  method = function(object) {
    if (object@M > object@m){
      stop("The parameter M must be smaller than the order for Bernstein estimator (m) and > 0")
    }else{}
    return(TRUE)
  }
  )

setGeneric (
  name = "getM",
  def  = function(.Object){standardGeneric("getM")}
)

setMethod(
  f = "getM",
  signature = "BernsteinPolynomials", 
  definition = function(.Object) {
    return(.Object@m)
  }
)

setMethod(
  f = "density", 
  signature = "BrVitale",
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
	
	########################################################################
	### auxiliar function to calculate density as in vitale approximation
	 density.vitale <- function(x,dataPoints,m){
		aux <- sapply(0:(m-1), FUN = 
			function(k,m,x){
				n <- length(dataPoints)
				Fn <- sum((dataPoints > k/m) & (dataPoints <= (k+1)/m)) / n
				return(Fn * dbeta(x,k+1,m-k))
		}, m = m, x = x)
	
		# if x.new contains more than one data points, aux is a matrix where the calculated values for each 
		# data point in x.new are stored by rows. However, if x.new is a single value, aux cointains a vector
		# with the calculated values the data point in x.new
		density <- numeric(0)
		if(is.matrix(aux)){
			densities <- rowSums(aux)  
		}else{
			densities <- sum(aux)
		}	
		
		return(densities)
	}	
	### end of auxiliar function to calculate density as in vitale approximation
	##############################################################################
	
		fm <- density.vitale(x.new,.Object@dataPoints,.Object@m)     
		fM <- density.vitale(x.new,.Object@dataPoints,.Object@M)
		densities <- (.Object@m/(.Object@m - .Object@M))*fm - (.Object@M/(.Object@m - .Object@M))*fM   			      
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


#####################################
## Constructor functions for users ##
#####################################

brVitale <- function(dataPoints, m=round(length(dataPoints)^(2/5)), M=NULL, dataPointsCache=NULL, lower.limit=0,upper.limit=1){
  #cat("~~~~~~ Bias Reduction Vitale: constructor ~~~~~~\n")
  
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
  
  if(is.null(M)){
    M <- m/2
  }
  polinomialModel <- new(Class="BrVitale",dataPoints = dataPoints.scaled, m = m, M = M, 
                         dataPointsCache = dataPointsCache.scaled,lower.limit=lower.limit,upper.limit=upper.limit)
  setDensityCache(polinomialModel, densityFunction=NULL)
  return(polinomialModel)
}

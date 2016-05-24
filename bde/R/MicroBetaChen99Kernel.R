#************************************************************************
# Kernel estimator with kernel-wise normalization as described in 
# Chen's 99 paper:
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
  Class = "MicroBetaChen99Kernel",
  representation = representation(	
	normalizationConstants = "numeric"
  ),
  prototype = prototype(
	modified = FALSE,
	normalizationConstants = numeric(0)
  ),
  contains = "Chen99Kernel"
  )

setValidity(
  Class = "MicroBetaChen99Kernel",
  method = function(object) {
    if (any(object@dataPoints == 0) | any(object@dataPoints == 1)){
		stop(paste("MicroBeta normalization presents problems when there are data samples in the boundaries (", object@lower.limit, 
			object@upper.limit,").\n",
 			"Use instead MacroBeta normalization or remove those data samples in the boundaries"))   
    }else{}
    return(TRUE)
  }
  )

setMethod(
  f = "setDataPoints",
  signature = "MicroBetaChen99Kernel", 
  definition = function(.Object,x) {    
    # obtain the global name of the variable to modify
    objectGlobalName <- deparse(substitute(.Object))
    
    # modify the variable localy
    .Object@dataPoints <- x    
    # The kernel model has changed -> the densityCache must be cleaned
    .Object@densityCache <- numeric(0)    
	# The kernel model has changed -> normalization constants also change
	.Object@normalizationConstants <- numeric(0)
    
    # assign the local variable to the global variable  
    assign(objectGlobalName,.Object,envir=parent.frame())      	
  }
  )
	

setGeneric (
  name = "setNormalizationConstants",
  def  = function(.Object){standardGeneric("setNormalizationConstants")}
  )

setMethod(
	f = "setNormalizationConstants",
	signature = "MicroBetaChen99Kernel",
	definition = function(.Object){
		
		##############################
		##### Auxiliar Functions #####
		##############################
		# calculate the density value at a point X_i of the Chen99's kernel function for both modified (K^*_{x,b}(X_i)) 
		# and no modified (K^*_{x/b + 1, (1-x)/b + 1}(X_i)) versions 
		density.kernelFunction <- function(x,Xi,b,modified){		
					
			shape1 <- numeric(0)
			shape2 <- numeric(0)

			if(.Object@modified){
				
				# data points in x so that they are < 2*b has the following shape paremeters
				index.cond1 <- which(x < 2*.Object@b)
				shape1[index.cond1] <- rho(.Object,x[index.cond1])
				shape2[index.cond1] <- (1-x[index.cond1]) / .Object@b
				# data points in x so that they are ( >= 2*b & <= 1-2*b) has the following shape paremeters
				index.cond2 <- which((x >= 2*.Object@b) & (x <= 1-2*.Object@b))
				shape1[index.cond2] <- x[index.cond2] / .Object@b
				shape2[index.cond2] <- (1-x[index.cond2]) / .Object@b
				# the rest of data points in x.new
				index.rest <- which((x > 1-2*.Object@b))
				shape1[index.rest] <- (x[index.rest] / .Object@b)
				shape2[index.rest] <- rho(.Object, 1-x[index.rest])		
			
			}else{
				# Not modified or modified and (x >= 2*.Object@b & x <= (1-2*.Object@b))
				shape1 <- x/b + 1
				shape2 <- (1-x)/b + 1	
			}	
		
			dbeta(Xi,shape1,shape2)					
		}		
		
		##################################
		##### End Auxiliar Functions #####
		##################################
		
		normalizationConstants <- numeric(0)
		
		
## In the boundaries (0 and 1) the integral of the kernel function between points 0 and 1 takes the value 0.
## It cause a wrong normalization of the density function which will not integrate to 1
## This is a problem with the MicroBeta normalization -> it is better to avoid these situations by removing 
## the datasamples in the boundary
boundaryDataPoints.indices <- (.Object@dataPoints == 0 | .Object@dataPoints == 1)
normalizationConstants[boundaryDataPoints.indices] = 0

		
		
		normalizationConstants[!boundaryDataPoints.indices] <- sapply(
			.Object@dataPoints[!boundaryDataPoints.indices], 
			FUN= function(x,b,modified){
				integrate(density.kernelFunction,0,1, Xi=x, b=b, modified=modified)$value
			},
			b = .Object@b, 
			modified = .Object@modified)
		
		
		# Save the normalization constants in the kernel and export the results to the global object
		# obtain the global name of the variable to modify
		objectGlobalName <- deparse(substitute(.Object))    
		# modify the variable
        .Object@normalizationConstants <- normalizationConstants
		# assign the local variable to the global variable  
		assign(objectGlobalName,.Object,envir=parent.frame())  	
})


setMethod(
  f = "density", 
  signature = "MicroBetaChen99Kernel",
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

    # if the normalization constants have not been calculated yet, we calculate them now
	if(length(.Object@normalizationConstants) == 0){
		# obtain the global name of the variable to modify
		objectGlobalName <- deparse(substitute(.Object))    
		# modify the variable
        setNormalizationConstants(.Object)
		# assign the local variable to the global variable  
		assign(objectGlobalName,.Object,envir=parent.frame())  			
	}
	
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
			index.rest <- which((x.new > 1-2*.Object@b))
			shape1[index.rest] <- (x.new[index.rest] / .Object@b)
			shape2[index.rest] <- rho(.Object, 1-x.new[index.rest])				
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
	
		# densities is a matrix with length(.Object@dataPoints) rows and x.new.length columns, to calculate the 
		# final density at the points in x.new we have to normalize and average over the columns.
 		# .Object@normalizationConstants is a vector with length(.Object@dataPoints) values (a normalization constant for
		# the kernel function evaluated at each point in .Object@dataPoints. By doing densities/.Object@normalizationConstants 		
		# we divide each column in densities by the normalization constants in .Object@normalizationConstants element by element.		
		x.densities[x.indices == 0] <- colSums(densities/.Object@normalizationConstants) / n	
	}

	# include the density (density=0) of the out-of-bound x points in the final result
	aux.density <- numeric(numDataPoints)
	aux.density[index.nozero] <- x.densities
	x.densities <- aux.density
    
	#if x is a matrix, we store the densities as a matrix object
	if(isMatrix.x){
		dim(x.densities) <- dims
}

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

microBetaChen99Kernel <- function(dataPoints, b=length(dataPoints)^(-2/5), dataPointsCache=NULL, modified = FALSE, lower.limit=0,upper.limit=1){
  #cat("~~~~~~ MicroBetaChen99Kernel: constructor ~~~~~~\n")  
  
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
  
  kernel <- new(Class="MicroBetaChen99Kernel",dataPoints = dataPoints.scaled, b = b, dataPointsCache = dataPointsCache.scaled, modified = modified, lower.limit=lower.limit,upper.limit=upper.limit)
  setNormalizationConstants(kernel)
	#density(kernel,0.5)
  setDensityCache(kernel, densityFunction=NULL)  
  return(kernel)
}

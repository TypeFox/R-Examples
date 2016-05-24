#************************************************************************
#  Boundary Kernel Density class with normalization as described in 
#  kakizawa's 2004 paper
#
# @article{kakizawa2004,
#	title = {Bernstein polynomial probability density estimation},
#	author = {Kakizawa, Y.},
#	journal = {Journal of nonparametric Statistics},
#	year = {2004},
#	volume = {16},
#	number = {5},
#	pages = {709-729}
#	}
#************************************************************************

setClass(
  Class = "NormalizedBoundaryKernel",
  representation = representation(),
  contains = "BoundaryKernel"
  )

setGeneric (
  name = "integrate.interiorKernelFunction",
  def  = function(.Object,lower,upper){standardGeneric("integrate.interiorKernelFunction")}
  )

setMethod(
	f = "singlePoint.density",
	signature = "NormalizedBoundaryKernel",
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
			#q <- (1-x)/.Object@b
			q <- (x-1)/.Object@b
			sum(rightBoundaryKernelFunction(.Object,q,u)) / (.Object@b*length(.Object@dataPoints))
		} else {}
		
})

setMethod(
	f = "integrate.interiorKernelFunction",
	signature = "NormalizedBoundaryKernel",
	definition = function(.Object,lower,upper){

		############## Auxiliar function #####################################################		
		# implements analytical formulation of the integral for the interior kernel function #	
		int <- function(u,mu){
			if (mu == 0){
				u/2
			}else if (mu == 1){
				(3*u)/4 - u^3/4
			}else if (mu == 2){
				15/16*(u+u^5/5-(2*u^3)/3)
			}else if (mu == 3){
				35/32*(u-u^3+(3*u^5)/5-u^7/7)
			}else {
				stop("Mu can only take values 0,1,2 or 3")
			}
		}	
		######################################################################################
	
	int(upper,.Object@mu) - int(lower,.Object@mu)
})


setMethod(
	f = "leftBoundaryKernelFunction",
	signature = "NormalizedBoundaryKernel",
	definition = function(.Object,q,u){	
		
		norm = integrate.interiorKernelFunction(.Object,lower=-1,upper=q)
		interiorKernelFunction(.Object,u)/norm
})

setMethod(
	f = "interiorKernelFunction",
	signature = "NormalizedBoundaryKernel",
	definition = function(.Object,u){
		if (.Object@mu == 0){
			rep(1/2,times=length(u))	
		}else if (.Object@mu == 1){
			(3/4)*(1-u^2)
		}else if (.Object@mu == 2){
			(15/16)*(1-u^2)^2
		}else if (.Object@mu == 3){
			(35/32)*(1-u^2)^3
		}else {
			stop("Mu can only take values 0,1,2 or 3")
		}	
})


setMethod(
	f = "rightBoundaryKernelFunction",
	signature = "NormalizedBoundaryKernel",
	definition = function(.Object,q,u){						
		norm = integrate.interiorKernelFunction(.Object,lower=q,upper=1)
		interiorKernelFunction(.Object,u)/norm		
})

#####################################
## Constructor functions for users ##
#####################################

normalizedBoundaryKernel <- function(dataPoints, mu=1, b=length(dataPoints)^(-2/5), dataPointsCache=NULL, lower.limit=0,upper.limit=1){
  #cat("~~~~~~ Normalized Boundary Kernel: constructor ~~~~~~\n")  
  
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
  
  kernel <- new(Class="NormalizedBoundaryKernel",dataPoints = dataPoints.scaled, b = b, 
                dataPointsCache = dataPointsCache.scaled, mu = mu, lower.limit=lower.limit,upper.limit=upper.limit)
  setDensityCache(kernel, densityFunction=NULL)
  return(kernel)
}

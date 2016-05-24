#************************************************************************
#  No Boundary Kernel Density class
#
#	This kernel uses the same kernel function for all left boundary
#	interior and right boundary regions
#************************************************************************

setClass(
  Class = "NoBoundaryKernel",
  representation = representation(),
  contains = "BoundaryKernel"
  )

setMethod(
	f = "leftBoundaryKernelFunction",
	signature = "NoBoundaryKernel",
	definition = function(.Object,q,u){	
		interiorKernelFunction(.Object,u)		
})

setMethod(
	f = "interiorKernelFunction",
	signature = "NoBoundaryKernel",
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
	signature = "NoBoundaryKernel",
	definition = function(.Object,q,u){
		interiorKernelFunction(.Object,u)
})

#####################################
## Constructor functions for users ##
#####################################

noBoundaryKernel <- function(dataPoints, mu=1, b=length(dataPoints)^(-2/5), dataPointsCache=NULL, lower.limit=0,upper.limit=1){
  #cat("~~~~~~ No Boundary Kernel: constructor ~~~~~~\n")  
  
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
  
  kernel <- new(Class="NoBoundaryKernel",dataPoints = dataPoints.scaled, b = b, dataPointsCache = dataPointsCache.scaled,
                mu = mu, lower.limit=lower.limit,upper.limit=upper.limit)
  setDensityCache(kernel, densityFunction=NULL)
  return(kernel)
}

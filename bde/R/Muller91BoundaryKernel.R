#************************************************************************
#  Boundary Kernel Density class as described in Muller91's paper
#
# @article{Muller1991,
#	title = {Smooth optimum kernel estimators near endpoints},
#	author = {M\"uller, H.},
#	journal = {Biometrika},
#	year = {1991},
#	volume = {79},
#	number = {3}
#	pages = {521-530}
#	}
#************************************************************************

setClass(
  Class = "Muller91BoundaryKernel",
  representation = representation(),
  contains = "BoundaryKernel"
  )

setMethod(
	f = "leftBoundaryKernelFunction",
	signature = "Muller91BoundaryKernel",
	definition = function(.Object,q,u){		
		if (.Object@mu == 0){
			1/(1+q)*(1 + 3*((1-q)/(1+q))^2 + 6*((1-q)/((1+q)^2))*u)		
		}else if (.Object@mu == 1){
			6*(1+u)*(q-u)*(1/((1+q)^3))*(1 + 5*((1-q)/((1+q))^2) + 10*((1-q)/((1+q)^2))*u)
		}else if (.Object@mu == 2){
			30*(1+u)^2*(q-u)^2*(1/((1+q)^5))*(1 + 7*((1-q)/(1+q))^2 + 14*((1-q)/((1+q)^2))*u)
		}else if (.Object@mu == 3){
			140*(1+u)^3*(q-u)^3*(1/((1+q)^7))*(1 + 9*((1-q)/(1+q))^2 + 18*((1-q)/((1+q)^2))*u)
		}else {
			stop("Mu can only take values 0,1,2 or 3")
		}	
})

setMethod(
	f = "interiorKernelFunction",
	signature = "Muller91BoundaryKernel",
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
	signature = "Muller91BoundaryKernel",
	definition = function(.Object,q,u){
		# the same computation as in the leftBoundaryKernel Function but evaluated at points -u 		
		leftBoundaryKernelFunction(.Object,q,-u)				
})

#####################################
## Constructor functions for users ##
#####################################

muller91BoundaryKernel <- function(dataPoints, mu=1, b=length(dataPoints)^(-2/5),dataPointsCache=NULL, lower.limit=0,upper.limit=1){
  #cat("~~~~~~ Muller91 Boundary Kernel: constructor ~~~~~~\n")  
  
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
  
  kernel <- new(Class="Muller91BoundaryKernel",dataPoints = dataPoints.scaled, b = b, dataPointsCache = dataPointsCache.scaled,
                mu = mu, lower.limit=lower.limit,upper.limit=upper.limit)
  setDensityCache(kernel, densityFunction=NULL)
  return(kernel)
}

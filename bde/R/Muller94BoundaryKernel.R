#************************************************************************
#  Boundary Kernel Density class as described in Muller94's paper
#
# @article{Muller1994,
#	title = {Hazard rate estimation under random censoring with varying
#			kernels and bandwidths},
#	author = {M\"uller, H. and Wang, J.},
#	journal = {Biometrics},
#	year = {1994},
#	volume = {50},
#	pages = {61-76}
#	}
#************************************************************************

setClass(
  Class = "Muller94BoundaryKernel",
  representation = representation(),
  contains = "BoundaryKernel"
  )

setMethod(
	f = "leftBoundaryKernelFunction",
	signature = "Muller94BoundaryKernel",
	definition = function(.Object,q,u){		
		if (.Object@mu == 0){
			2/((1+q)^3)*(3*(1-q)*u+2*(1-q+q^2))
		}else if (.Object@mu == 1){
			12/((1+q)^4)*(u+1)*((1-2*q)*u+(1-2*q+3*q^2)*0.5)
		}else if (.Object@mu == 2){
			15*(1+u)^2*(q-u)*(1/(1+q)^5)*(2*u*(5*((1-q)/(1+q))-1)+3*q-1+(5*(1-q)^2)/(1+q))
		}else if (.Object@mu == 3){
			70*(1+u)^3*(q-u)^2*(1/(1+q)^7)*(2*u*(7*((1-q)/(1+q))-1)+3*q-1+(7*(1-q)^2)/(1+q))
		}else {
			stop("Mu can only take values 0,1,2 or 3")
		}	
})

setMethod(
	f = "interiorKernelFunction",
	signature = "Muller94BoundaryKernel",
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
	signature = "Muller94BoundaryKernel",
	definition = function(.Object,q,u){
		# the same computation as in the leftBoundaryKernel Function but evaluated at points -u 		
		leftBoundaryKernelFunction(.Object,q,-u)	
})

#####################################
## Constructor functions for users ##
#####################################

muller94BoundaryKernel <- function(dataPoints, mu=1, b=length(dataPoints)^(-2/5), dataPointsCache=NULL, lower.limit=0,upper.limit=1){
  #cat("~~~~~~ Muller94 Boundary Kernel: constructor ~~~~~~\n")
  
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
  
  kernel <- new(Class="Muller94BoundaryKernel",dataPoints = dataPoints.scaled, b = b, 
                dataPointsCache = dataPointsCache.scaled, mu = mu, lower.limit=lower.limit,upper.limit=upper.limit)
  setDensityCache(kernel, densityFunction=NULL)
  return(kernel)
}


#************************************************************************
#           Kernel Density class (Virtual Class)
#************************************************************************

setClass(
  Class = "KernelDensity",
  representation = representation(
    dataPoints = "numeric",
    b = "numeric",
    "VIRTUAL"), 
  prototype = prototype(dataPointsCache=seq(0,1,0.01)),
  contains = "BoundedDensity"
  )

setValidity(
  Class = "KernelDensity",
  method = function(object) {
    if (object@b < 0 | object@b > 1){
      stop("The parameter b must be 0<=b<=1")
    }else{}
    return(TRUE)
  }
  )

setMethod(
  f = "print",
  signature = "KernelDensity", 
  definition = function(x,...) {  
	str(x)
  }
  )

setMethod(
  f = "show",
  signature = "KernelDensity", 
  definition = function(object) {    
	str(object)
  }
  )

setGeneric (
  name = "getdataPoints",
  def  = function(.Object){standardGeneric("getdataPoints")}
  )

setMethod(
  f = "getdataPoints",
  signature = "KernelDensity", 
  definition = function(.Object) {
    return(getOriginalPoints(.Object,.Object@dataPoints))
  }
  )

setGeneric (
  name = "getb",
  def  = function(.Object){standardGeneric("getb")}
  )

setMethod(
  f = "getb",
  signature = "KernelDensity", 
  definition = function(.Object) {
    return(.Object@b)
  }
  )

setGeneric (
  name = "setDataPoints",
  def  = function(.Object,x){standardGeneric("setDataPoints")}
  )

setMethod(
  f = "setDataPoints",
  signature = "KernelDensity", 
  definition = function(.Object,x) {  
    
    if(.Object@upper.limit != 1 && .Object@lower.limit != 0){
      # scale the data to the 0-1 interval
      x <- getScaledPoints(.Object,x)    
    }
    
    # obtain the global name of the variable to modify
    objectGlobalName <- deparse(substitute(.Object))
    
    # modify the variable localy
    .Object@dataPoints <- x    
    # The kernel model has changed -> the densityCache must be cleaned
    .Object@densityCache <- numeric(0)    
    
    # assign the local variable to the global variable  
    assign(objectGlobalName,.Object,envir=parent.frame())      	
  }
  )

setGeneric (
  name = "setB",
  def  = function(.Object,b){standardGeneric("setB")}
  )

setMethod(
  f = "setB",
  signature = "KernelDensity", 
  definition = function(.Object,b) {    
    # obtain the global name of the variable to modify
    objectGlobalName <- deparse(substitute(.Object))
    
    # modify the variable localy
    .Object@b <- b    
    # The kernel model has changed -> the densityCache must be cleaned
    .Object@densityCache <- numeric(0)    
    
    # assign the local variable to the global variable  
    assign(objectGlobalName,.Object,envir=parent.frame())    
  }
  )

setMethod(
  f = "setDensityCache",
  signature = "KernelDensity", 
  definition = function(.Object, densityFunction=NULL) {  
    # This method is only for internal use, therefore it is not affected by the use of differente scale
    # (different from 0-1)
    
    # obtain the global name of the variable to modify
    objectGlobalName <- deparse(substitute(.Object))
	
    # modify the variable localy
	if(is.null(densityFunction)){
      .Object@densityCache <- density(.Object, .Object@dataPointsCache, scaled=TRUE)
    }else{
      cat("WARNING: A new density function is provided and it may be different from the kernel\n")
      cat("density function. Use preferably the parameter densityFunction=null (default)")
      .Object@densityCache <- sapply(.Object@dataPointsCache,densityFunction)          
    }
    
    # assign the local variable to the global variable  
    assign(objectGlobalName,.Object,envir=parent.frame())        
  }
  )


setMethod(
  f = "distribution", 
  signature = "KernelDensity",
  definition = function(.Object, x,... , discreteApproximation=TRUE, scaled=FALSE) {
   
  if(!scaled){
      # scale the data to the 0-1 interval
      x <- getScaledPoints(.Object,x)    
  }
    
	cummulative.distribution.x <- numeric(0)
	
	if(discreteApproximation){
 
      if(length(.Object@densityCache) == 0){
		# obtain the global name of the variable to modify
		objectGlobalName <- deparse(substitute(.Object))    
		# modify the variable
        setDensityCache(.Object,densityFunction=NULL)
		# assign the local variable to the global variable  
		assign(objectGlobalName,.Object,envir=parent.frame())  
      }else{
        #do nothing (densityCache is already set)
      }
      
      cummulative.distribution.x <- discreteApproximationToDistribution(.Object,x)
    }else{
      f <- function(x,kernel){
		#cat(length(x),sep="\n")
		density(kernel,x,scaled=TRUE)				
      }
	
		if(is.null(dim(x))){#x is a vector or a number
			cummulative.distribution.x <- sapply(x, function(x){integrate(f, 0, x, kernel=.Object, subdivisions=1000)$value})
		}else{# x is a matrix or other object with 2 dimensions
			cummulative.distribution.x <- apply(x, MARGIN=c(1,2), FUN = function(x) 
				{integrate(f, 0, x, kernel=.Object, subdivisions=10000)$value})
		}  		
    }
    return(cummulative.distribution.x)
  }
  )


setMethod(
  f = "plot",
  signature = "KernelDensity", 
  definition = function(x,y,...) {
    if(length(x@densityCache) == 0){
		# obtain the global name of the variable to modify
		objectGlobalName <- deparse(substitute(x))
    	# modify the variable localy 
      	setDensityCache(x)
		# assign the local variable to the global variable
      	assign(objectGlobalName,x,envir=parent.frame()) 
    }  
    callNextMethod()
	return(invisible())
  }
  )
#************************************************************************
#           Bernsteins Polynomials class (Virtual Class)
#************************************************************************

setClass(
  Class = "BernsteinPolynomials",
  representation = representation(
	dataPoints = "numeric",
    m = "numeric",
    "VIRTUAL"), 
  prototype = prototype(dataPointsCache=seq(0,1,0.01)),
  contains = "BoundedDensity"
  )

setValidity(
  Class = "BernsteinPolynomials",
  method = function(object) {
    isWholeNumber <- function(x){
	dec <- x-round(x)
	if(dec > 0){
		return(FALSE)
	}else{}
	
	return(TRUE)		
	}
	
	if (!isWholeNumber(object@m) | object@m == 0){
      stop("The parameter m must be a integer grater than 0")
    }else{}
    return(TRUE)
  }
)

setMethod(
  f = "print",
  signature = "BernsteinPolynomials", 
  definition = function(x,...) {
    callNextMethod()
#     cat("*** The method also cotains \n")
#     cat("* m = "); print(x@m); cat("\n")     
#     cat("******* End Print (BernsteinPolynomials) ******* \n")
  }
  )

setMethod(
  f = "show",
  signature = "BernsteinPolynomials", 
  definition = function(object) {
    callNextMethod()
#     nDataShown <- min(10,length(object@dataPointsCache))
#     cat("*** The method also cotains \n")
#     cat("* m = "); print(object@b); cat("\n")    
#     cat("******* End show (BernsteinPolynomials) ******* \n")
  }
  )

setGeneric (
  name = "getdataPoints",
  def  = function(.Object){standardGeneric("getdataPoints")}
  )

setMethod(
  f = "getdataPoints",
  signature = "BernsteinPolynomials", 
  definition = function(.Object) {
    return(.Object@dataPoints)
  }
  )

setGeneric (
  name = "getm",
  def  = function(.Object){standardGeneric("getm")}
  )

setMethod(
  f = "getm",
  signature = "BernsteinPolynomials", 
  definition = function(.Object) {
    return(.Object@m)
  }
  )


setMethod(
  f = "setDensityCache",
  signature = "BernsteinPolynomials", 
  definition = function(.Object, densityFunction=NULL) {  
    # obtain the global name of the variable to modify
    objectGlobalName <- deparse(substitute(.Object))
	
    # modify the variable localy
	if(is.null(densityFunction)){
     .Object@densityCache <- density(.Object, .Object@dataPointsCache, scaled=T)
    }else{
      cat("WARNING: A new density function is provided and it may be different from the current Bernstein Polynomial model\n")
      cat("density function. Use preferably the parameter densityFunction=null (default)")
      .Object@densityCache <- sapply(.Object@dataPointsCache,densityFunction)          
    }
    
    # assign the local variable to the global variable  
    assign(objectGlobalName,.Object,envir=parent.frame())        
  }
  )

setMethod(
  f = "distribution", 
  signature = "BernsteinPolynomials",
  definition = function(.Object, x, ..., discreteApproximation=TRUE, scaled=FALSE) {
    
    if(!scaled){
      # scale the data to the 0-1 interval
      x <- getScaledPoints(.Object,x)    
    }

    if(discreteApproximation){
 
      if(length(.Object@densityCache) == 0){
		# obtain the global name of the variable to modify
		objectGlobalName <- deparse(substitute(.Object))    
		# modify the variable
        setDensityCache(.Object,densityFunction=NULL)
		# assign the local variable to the global variable  
		assign(objectGlobalName,.Object,envir=parent.frame())          
      }else{}
      
      cummulativeDistribution.x <- callNextMethod(.Object, x, discreteApproximation=TRUE,scaled=TRUE)
    }else{
      f <- function(x,.Object){
        density(.Object,x,scaled=TRUE)		
      }
	
		if(!is.matrix(x)){
			cummulativeDistribution.x <- sapply(x, function(x) 
				integrate(f,0,x, .Object=.Object)$value)
		}else{
			cummulativeDistribution.x <- apply(x, MARGIN=c(1,2), FUN = function(x) 
				integrate(f,0,x, .Object=.Object)$value)
		}  		
    }
    return(cummulativeDistribution.x)
  }
  )

setMethod(
  f = "plot",
  signature = "BernsteinPolynomials", 
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

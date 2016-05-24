#           Bounded Density class
#************************************************************************

setClass(
  Class = "BoundedDensity",
  representation = representation(
    dataPointsCache = "numeric", 
    densityCache = "numeric",
    distributionCache = "numeric",
    lower.limit = "numeric",
    upper.limit = "numeric"
  )
)


setValidity(
  Class = "BoundedDensity",
  method = function(object) {
    if (length(object@dataPointsCache) == 0){
      stop("ERROR: A data set with at least one point is needed")
    }else if (any(object@dataPointsCache < 0) | any(object@dataPointsCache > 1)){
      stop("ERROR: Data points outside the bounds [lower.limit,upper.limit]")
    }else if ( (length(object@densityCache) != 0) & 
                 (length(object@densityCache) != length(object@dataPointsCache)) ){
      stop("ERROR: dataPointsCache and densityCache must have the same dimension")
    }else {
      
      # data shold be stored in dataPointsCache in ascending order
      if(is.unsorted(object@dataPointsCache)){
        cat("WARNING: data shold be stored in dataPointsCache in ascending order.")
        cat("** dataPointsCache and densityCache have been sorted out")
        dataPoints.order <- order(object@dataPointsCache)
        object@dataPointsCache <- object@dataPointsCache[dataPoints.order]
        # if densities for dataPoints in x are given, they are also sorted      
        if(length(object@dataPointsCache) == length(object@densityCache)){
          object@densityCache <- object@densityCache[dataPoints.order]
        }
      }
      
    }
    return(TRUE)
  }
)

setMethod(
  f = "print",
  signature = "BoundedDensity", 
  definition = function(x,...) {
    str(x)
  }
)

setMethod(
  f = "show",
  signature = "BoundedDensity", 
  definition = function(object) {  
    str(object)
  }
)

setGeneric (
  name = "getScaledPoints",
  def  = function(.Object,x){standardGeneric("getScaledPoints")}
)

setMethod(
  f = "getScaledPoints",
  signature = "BoundedDensity", 
  definition = function(.Object,x) {
	value <- (x-.Object@lower.limit)/(.Object@upper.limit - .Object@lower.limit)
    return(value)
  }
)

setGeneric (
  name = "getOriginalPoints",
  def  = function(.Object,x){standardGeneric("getOriginalPoints")}
)

setMethod(
  f = "getOriginalPoints",
  signature = "BoundedDensity", 
  definition = function(.Object,x) {
    return((.Object@upper.limit - .Object@lower.limit)*x + .Object@lower.limit)
  }
    )

setGeneric (
  name = "getdataPointsCache",
  def  = function(.Object){standardGeneric("getdataPointsCache")}
)

setMethod(
  f = "getdataPointsCache",
  signature = "BoundedDensity", 
  definition = function(.Object) {
	data <- getOriginalPoints(.Object,.Object@dataPointsCache) 
    return(data)
  }
)


setGeneric (
  name = "getdensityCache",
  def  = function(.Object){standardGeneric("getdensityCache")}
)

setMethod(
  f = "getdensityCache",
  signature = "BoundedDensity", 
  definition = function(.Object) {
    domain.length <- .Object@upper.limit - .Object@lower.limit
	return(.Object@densityCache/domain.length)
  }
)

setGeneric (
  name = "getdistributionCache",
  def  = function(.Object){standardGeneric("getdistributionCache")}
)

setMethod(
  f = "getdistributionCache",
  signature = "BoundedDensity", 
  definition = function(.Object) {
	domain.length <- .Object@upper.limit <- .Object@lower.limit
    return(.Object@distributionCache/domain.length)
  }
)

setGeneric (
  name = "setDataPointsCache",
  def  = function(.Object,x){standardGeneric("setDataPointsCache")}
)

setMethod(
  f = "setDataPointsCache",
  signature = "BoundedDensity", 
  definition = function(.Object,x) {
    
    if(.Object@upper.limit != 1 && .Object@lower.limit != 0){
      # scale the data to the 0-1 interval
      x <- getScaledPoints(.Object,x)    
    }
    
    # obtain the global name of the variable to modify
    objectGlobalName <- deparse(substitute(.Object))
    
    # modify the variable localy    
    # dataPointsCache is stored in ascending order
    if(is.unsorted(x)){
      x.order <- order(x)
      .Object@dataPointsCache <- x[x.order]      
      
    }else{
      .Object@dataPointsCache <- x
    }
    
    # clean the densityCache
    .Object@densityCache <- numeric(0)
    
    if(validObject(.Object)){
      # assign the local variable to the global variable  
      assign(objectGlobalName,.Object,envir=parent.frame()) 
    }else{}
  }
)

setGeneric (
  name = "setDensityCache",
  def  = function(.Object,densityFunction=NULL){standardGeneric("setDensityCache")}
)

setMethod(
  f = "setDensityCache",
  signature = "BoundedDensity", 
  definition = function(.Object,densityFunction) { 
    # We need a density function to calculate the densities
    if(!is.function(densityFunction)){
      stop("A density function is needed to calculate the densities of the data points in dataPointsCache\n")
    }
    
    # obtain the global name of the variable to modify
    objectGlobalName <- deparse(substitute(.Object))
    
    # modify the variable localy
    .Object@densityCache <- sapply(.Object@dataPointsCache,densityFunction)    
    
    if(validObject(.Object)){
      # assign the local variable to the global variable  
      assign(objectGlobalName,.Object,envir=parent.frame())  
    }else{}
    
  }
)

setGeneric (
  name = "setDistributionCache",
  def  = function(.Object){standardGeneric("setDistributionCache")}
)

setMethod(
  f = "setDistributionCache",
  signature = "BoundedDensity", 
  definition = function(.Object) { 
    # We need the density cache to calculate the distribution values
    if(length(.Object@dataPointsCache) != length(.Object@densityCache)){
      stop("A density cache is needed to calculate the distribution values of the data points in dataPointsCache\n")
    }
    
    # obtain the global name of the variable to modify
    objectGlobalName <- deparse(substitute(.Object))
    
    # modify the variable localy
    .Object@distributionCache <- discreteApproximationToDistribution(.Object,.Object@dataPointsCache)    
    
    if(validObject(.Object)){
      # assign the local variable to the global variable  
      assign(objectGlobalName,.Object,envir=parent.frame())  
    }else{}
    
  }
)

setGeneric (
  name = "forceDensityCacheTo",
  def  = function(.Object,densities){standardGeneric("forceDensityCacheTo")}
)

setMethod(
  f = "forceDensityCacheTo",
  signature = "BoundedDensity", 
  definition = function(.Object,densities) {
    #cat("WARNING: This is a method that can easily lead to programing errors. Use setDensityCache instead when possible\n")
    
    # obtain the global name of the variable to modify
    objectGlobalName <- deparse(substitute(.Object))
    
    # modify the variable localy
    .Object@densityCache <- densities    
    
    if(validObject(.Object)){
      # assign the local variable to the global variable  
      assign(objectGlobalName,.Object,envir=parent.frame())  
    }else{}
    
  }
)

if(!isGeneric("density")){
setGeneric (
  name = "density",
  def  = function(x,...){standardGeneric("density")}
  #These parameters agree with those from the density method defined in the stats package
 )
}


# Evaluates the boundedDensity in a single point using the information from densityCache
setMethod(
  f = "density",
  signature = "BoundedDensity", 
  definition = function(x,values,scaled = FALSE) {    
    #These parameters agree with those from the density method defined in the stats package
    .Object <- x
    x <- values
    
    isMatrix.x <- is.matrix(x)
    #dims = [nrows,ncols]
    dims <- dim(x)    
   
    
    # the scaled parameter controls the scale of the data:
    #  * if scaled = TRUE we assume that the points in x are scaled between 0 and 1
    #  * if scaled = FALSE se assume that the points in x are scaled between .Object@lower.limit and .Object@upper.limit
    if(!scaled){
      x <- getScaledPoints(.Object,x)
    }
	
	# if any value in x is lower than 0 or grater than 1 its density is 0
	numDataPoints <- length(x)
	index.nozero <- which(x>=0 & x <=1)
	
	x <- x[index.nozero]
	if(length(x) == 0){ # all elements in x where out of bound
		return(rep(0,numDataPoints - length(index.nozero)))
	}
	
    
    density.x <- numeric(length(x))
    
    if(length(.Object@dataPointsCache) == 0){
      stop("No data points in the cache")
      return(numeric(0))
    }
    else if(length(.Object@dataPointsCache) == 1){
      # There is only a point in dataPointsCache: We assume the same density for any other 
      # point in [lower.limit,upper.limit]
      
      density.x <- rep(.Object@densityCache[1],length(x))
    }else{
      # we calculate the density using the slope of the line between the previous and the subsequent points (in the cache) 
      # to the values in x
      a <- numeric(0) # heigh of the triangle formed by the preceding and the subsequent points to the values in x
      b <- numeric(0) # base of the triangle formed by the preceding and the subsequent points to the values in x
      b.x <- numeric(0) # base of the triangle formed by x and it preceding value in the dataPointsCache
      a.x <- numeric(0) # height of the triangle formed by x and it preceding value in the dataPointsCache
      
      # when x is a value in dataPointsCache. We take its density from densityCache
      #  indices.dataPointsCache: The position in the cache where the values x[indices.equal] are stored
      indices.dataPointsCache <- match(x, .Object@dataPointsCache)		
	  # From the previous "match" command, if x is not in the cache, the corresponding index is a "NA" value.
 	  # indices.equal: the values from x which are already calculated in the cache.
      indices.x.equal <- which(!is.na(indices.dataPointsCache))
      # We should remove these NA's from indices.dataPointsCache
      indices.dataPointsCache <- indices.dataPointsCache[!is.na(indices.dataPointsCache)]
     
      density.x[indices.x.equal] <- .Object@densityCache[indices.dataPointsCache]
      # a,b and b.x are set to the following values just to allow a correct calculculation for the density of all the data
      # points in x at the same time. See at the end of the function:
      ##	density.x <- density.x + a.x
      # Thus, for those values of x which are in the dataPointsCache a.x will be 0 and therefore the density of x will be the 
      # density value form the densityCache
      a[indices.x.equal] <- 0
      b[indices.x.equal] <- 1
      b.x[indices.x.equal] <- 0
      
      # when x is smaller than any other value in dataPointsCache. The trend observed in the first two points in 
      # dataPointsCache is used to calculate the slope   
      indices.x.less <- which(x < .Object@dataPointsCache[1])
      
      if(length(indices.x.less) > 0){
        # we multiply by (-1) since we need to substract a quantity to .Object@densityCache[1] in order to obtain the density of x
        # Thus, we will directly obtain a negative value
        a[indices.x.less] <- (-1) * (.Object@densityCache[2] - .Object@densityCache[1])
        b[indices.x.less] <- .Object@dataPointsCache[2] - .Object@dataPointsCache[1]
        b.x[indices.x.less] <- .Object@dataPointsCache[1] - x[indices.x.less]      
        density.x[indices.x.less] <- .Object@densityCache[1]
      }else{}
      
      # when x is greater than any other value in dataPointsCache. The trend observed in the last two points in
      # dataPointsCache is used to calculate the slope
      indices.x.great <- which(x > .Object@dataPointsCache[length(.Object@dataPointsCache)])		
      n <- length(.Object@dataPointsCache)
      
      if(length(indices.x.great) > 0){
        a[indices.x.great] <- .Object@densityCache[n] - .Object@densityCache[n-1]
        b[indices.x.great] <- .Object@dataPointsCache[n] - .Object@dataPointsCache[n-1]
        b.x[indices.x.great] <- x[indices.x.great] - .Object@dataPointsCache[n]
        density.x[indices.x.great] <- .Object@densityCache[n]
      }else{}
      
      # when x is between two values in dataPointsCache. The trend observed among the preceding and the 
      # subsequent values is used to calculate the slope
      indices.x.rest <- 1:length(x)
      indices.x.used <- c(indices.x.equal, indices.x.less, indices.x.great)
      
      if(length(indices.x.used) > 0){
        indices.x.rest <- indices.x.rest[-indices.x.used]
      }else{}
      
      if(length(indices.x.rest) > 0){
        prec.x <- sapply(x[indices.x.rest], function(x)  max(which(.Object@dataPointsCache < x)))
        subseq.x <- prec.x + 1
        a[indices.x.rest] <- .Object@densityCache[subseq.x] - .Object@densityCache[prec.x]
        b[indices.x.rest] <- .Object@dataPointsCache[subseq.x] - .Object@dataPointsCache[prec.x]
        b.x[indices.x.rest] <- x[indices.x.rest] - .Object@dataPointsCache[prec.x]
        density.x[indices.x.rest] <- .Object@densityCache[prec.x]
      }
      
      tan.alpha <- a/b
      a.x <- tan.alpha * b.x
      density.x <- density.x + a.x		
    }
	
	# due to discrete approximation to density function the obtained values may be negative. We substitute negative values
	# by 0
	density.x[density.x < 0] <- 0
	
    aux.density <- numeric(numDataPoints)
	aux.density[index.nozero] <- density.x
	
	#if x is a matrix, we store the densities as a matrix object
    if(isMatrix.x){
        dim(aux.densities) <- dims
    }
  
	#if data are in another scale (not in the [0,1] interval) we should normalize the density by dividing it by the length
 	# of the interval so that the density integrates to 1
	domain.length <- .Object@upper.limit - .Object@lower.limit
    if(!scaled){
      aux.density <- aux.density/domain.length
    }
    
	return(aux.density)     
    
  }
)

setGeneric (
  name = "discreteApproximationToDistribution",
  def  = function(.Object,x, ...){standardGeneric("discreteApproximationToDistribution")}
)

# due to discrete approximation, discreteApproximationToDistribution may return distribution values greater than 1
# specially when there are a few data points in dataPointsCache. 
setMethod(
  f = "discreteApproximationToDistribution", 
  signature = "BoundedDensity",
  definition = function(.Object, x, ...){       
    
	# sort x 
    order.x <- order(x)
	x <- x[order.x]
	
	if(length(.Object@dataPointsCache) == 1){
      
      cummulative.distribution <- x * .Object@densityCache[1] 
    }
	else{
      
      # distribution.for.a.single.point calculates de cummulated distribution value for a single point x, 
      # this is used latter to calculate the distribution for a vector or a matrix of points
      distribution.for.a.single.point <- function(x){			
        # the cummulated distribution is approximated by calculating the ara under the curve
        # described by the points in dataPointsCache
        # First, we approximate the area using bins (w: width, h: height)
        # no. of points in the cache which are smaller than X  
        w <- numeric(0)
        h <- numeric(0)	
        n.smallerThanX <- length(which(.Object@dataPointsCache < x))
        
        if(n.smallerThanX == 0){
		  # add the point x, it is smaller than any other point in the dataPointsCache
          w <- x
          h <- density(.Object,x)
		}
		else if(n.smallerThanX == 1){
          w <- .Object@dataPointsCache[1] #- 0
          h <- .Object@densityCache[1]
          # add the point x
          w <- c(w,(x-.Object@dataPointsCache[1]))
          h <- c(h,density(.Object,x))	
        }
		else{ #(n.smallerThanX > 1)
          w <- .Object@dataPointsCache[1:n.smallerThanX] - c(0,.Object@dataPointsCache[1:n.smallerThanX-1])
          h <- .Object@densityCache[1:n.smallerThanX]
          # add the point x
          w <- c(w,(x-.Object@dataPointsCache[n.smallerThanX]))
          h <- c(h,density(.Object,x))	
        }
        
        cummulative.distribution.pointX <- sum(w*h)
        
        # Finally, the approximation is refined (b:base ,h:height of the triangle)
        b <- w
        if(n.smallerThanX == 0){
		  # the point x is smaller than any other point in the dataPointsCache
		  h <- density(.Object,x) - density(.Object,0)			
		}
		else if(n.smallerThanX == 1){
		  h <- .Object@densityCache[1] - density(.Object,0)
          # add to h the corresponding value for the point x          
          h <- c(h, density(.Object,x) - .Object@densityCache[1])	
		}
        else #(n.smallerThanX > 1)
		{
          h <- .Object@densityCache[1:n.smallerThanX] - c(density(.Object,0),.Object@densityCache[1:n.smallerThanX-1])
          # add to h the corresponding value for the point x
          h <- c(h,density(.Object,x)-.Object@densityCache[n.smallerThanX])
		}
       
				
        
        cummulative.distribution.pointX <- cummulative.distribution.pointX - sum(b*h/2)					
        return(cummulative.distribution.pointX)
      }
      
      cummulative.distribution <- sapply(x,distribution.for.a.single.point)			
    }
	
	# set the original order
	cummulative.distribution[order.x] <- cummulative.distribution
    
    #if x is a matrix, we storte the cummulative.distribution as a matrix object
    if(is.matrix(x)){
      dim(cummulative.distribution) <- c(nrow(x),ncol(x))
    }  
    
    return(cummulative.distribution)
  }
)

if(!isGeneric("distribution")){
setGeneric (
  name = "distribution",
  def  = function(.Object, ...){standardGeneric("distribution")}
)
}

setMethod(
  f = "distribution", 
  signature = "BoundedDensity",
  definition = function(.Object, x, discreteApproximation = TRUE,scaled = FALSE) {     
        
    # the scaled parameter controls the scale of the data:
    #  * if scaled = TRUE we assume that the points in x are scaled between 0 and 1
    #  * if scaled = FALSE se assume that the points in x are scaled between .Object@lower.limit and .Object@upper.limit
    
    #print(c("Estoy en distribution y los valores son discreteApproximation=", discreteApproximation, "y scaled = ",scaled))
    
    if(!scaled){
      x <- getScaledPoints(.Object,x)
    }
        
    cummulative.distribution <- discreteApproximationToDistribution(.Object,x)
    
    # due to discrete approximation, cummulative distribution values may be greater than 1
    # specially when there are a few data points in dataPointsCache. We set those values
    # to 1
    # cummulative.distribution[cummulative.distribution > 1] <- 1	
    if(any(cummulative.distribution > 1)){
      cat("WARNING: due to discrete approximation, distribution value may be greater than 1\n")
    }
    
    return(cummulative.distribution)  	
  }
)

if(!isGeneric("quantile")){
setGeneric (
  name = "quantile",
  def = function(x,...){standardGeneric("quantile")}
  #These parameters agree with those from the density method defined in the stats package
)
}

setMethod(
  f= "quantile", 
  signature = "BoundedDensity",
  definition = function(x, p){  
    #These parameters agree with those from the density method defined in the stats package
    .Object <- x
    
    # if distribution cache is empty, calculate the values
    if(length(.Object@distributionCache) == 0){
      # obtain the global name of the variable to modify
      objectGlobalName <- deparse(substitute(.Object))
      
      # modify the variable localy
      setDistributionCache(.Object)
      
      if(validObject(.Object)){
        # assign the local variable to the global variable  
        assign(objectGlobalName,.Object,envir=parent.frame())  
      }else{}      
    } 
    
	# A bounded Density estimator may not integrate to 1. Therefore, we consider the quantile in a 
	# scale [0,distribution(.Object,1)] instead of in the scale [0,1]
	F.max <- distribution(.Object,1)
	p <- p*F.max
	
	
    index.lowerBound <- sapply(p,
                               FUN=function(x,.Object){sum(.Object@distributionCache <= x)},
                               .Object = .Object)
    index.upperBound <- index.lowerBound + 1
    #print(index.lowerBound)
    # if p is == last element in .Object@distributionCache we use the last element for both
    index.upperBound[index.lowerBound == length(.Object@distributionCache)] <- length(.Object@distributionCache)
    
    # If index.lowerBound contains 0s we may have problems when accessing .Object@dataPointsCache[index.lowerBound]
	dataPoints.lowerBound <- array(0,length(index.lowerBound))
	distribution.lowerBound <- array(0,length(index.lowerBound))
	not.zero <- (index.lowerBound != 0)
	dataPoints.lowerBound[not.zero] <- .Object@dataPointsCache[index.lowerBound]
	distribution.lowerBound[not.zero] <- .Object@distributionCache[index.lowerBound]
		
	inc.dataPoints <- .Object@dataPointsCache[index.upperBound] - dataPoints.lowerBound
    inc.distribution <- .Object@distributionCache[index.upperBound] - distribution.lowerBound
    inc.distributionNeeded <- p - distribution.lowerBound
    inc.dataPointNeeded <-  (inc.dataPoints * inc.distributionNeeded) / inc.distribution
    
    quantile.value <- dataPoints.lowerBound + inc.dataPointNeeded
    # re-scale the quantile values from [0,1] interval to their original scale
	return(getOriginalPoints(.Object,quantile.value))
  })

setGeneric (
  name = "rsample",
  def  = function(.Object,n){standardGeneric("rsample")}
)

setMethod(
  f= "rsample", 
  signature = "BoundedDensity",
  definition = function(.Object, n){    
    s <- runif(n,min=0,max=distribution(.Object,1,discreteApproximation=TRUE))
    quantile(.Object,s)
  }
)


setMethod(
  f = "plot",
  signature = "BoundedDensity", 
  definition = function(x,main="Bounded density",type="l",xlab="X",ylab="Density",...) {  
    plot(x = getdataPointsCache(x), y = getdensityCache(x), main=main,xlab=xlab, ylab=ylab, type=type,...)    
  }
)

setMethod(
  f = "lines",
  signature = "BoundedDensity", 
  definition = function(x,...) {  
    
    lines(x = getdataPointsCache(x), y = getdensityCache(x),, ...)    
  }
)
setGeneric (
  name = "gplot",
  def  = function(.Object,show=F,includePoints=F,lwd=1,alpha=1){standardGeneric("gplot")}
)

setMethod(
  f= "gplot", 
  signature = "list",
  definition = function(.Object,show,includePoints,lwd,alpha){    
    aux<-lapply(1:length(.Object),FUN=function(x){
      pts=getdataPointsCache(.Object[[x]])
      dens=getdensityCache(.Object[[x]])
      lab=rep(names(.Object)[x],length(dens))
      data.frame(X=pts,Density=dens,Label=lab)})
    df<-do.call(rbind,aux)
    X <- df$X
    p<-ggplot(df,aes_string(x="X",y="Density",col="Label")) + geom_line(size=lwd) + scale_color_discrete("")
    if(includePoints){
      wd<-(max(df$Density)-min(df$Density))
      y<-min(df$Density)-0.1*wd
      df2<-data.frame(x=df$X,y=y,col=df$Label)
      p<-p + geom_point(data=df2,mapping=aes_string(x="x",y="y",col="col"),position=position_jitter(height=0.05*wd),alpha=alpha)
    }
    if (show) print(p)
    return(p)
  }
)


setMethod(
  f= "gplot", 
  signature = "BoundedDensity",
  definition = function(.Object,show,includePoints,lwd,alpha){         
    df<-data.frame(X=getdataPointsCache(.Object), Density = getdensityCache(.Object))
    p<- ggplot(df,aes_string(x="X",y="Density")) + geom_line(size=lwd)
    if(includePoints){
      wd<-(max(.Object@densityCache)-min(.Object@densityCache))
      y<-min(.Object@densityCache)-0.1*wd       
      df2<-data.frame(x=getOriginalPoints(.Object,.Object@dataPoints) ,y=y)
      p<-p + geom_point(data=df2,mapping=aes_string(x="x",y="y"),position=position_jitter(height=0.05*wd),alpha=alpha)
    }
    if (show) print(p)
    return(p)
  }
)


# #####################################
# ## Constructor functions for users ##
# #####################################
# 
boundedDensity <- function(x,densities,lower.limit=0,upper.limit=1){
  #cat("~~~~~~ BoundedDensity: constructor ~~~~~~\n")
  x.scaled <- x
  if(lower.limit!=0 || upper.limit!=1){
    x.scaled <- (x-lower.limit)/(upper.limit-lower.limit)
  }
  new(Class="BoundedDensity",dataPointsCache = x.scaled,densityCache = densities, lower.limit=lower.limit,upper.limit=upper.limit)  
}

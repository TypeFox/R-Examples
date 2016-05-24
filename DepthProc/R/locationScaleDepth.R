# All functions realted with lsd package is here
# Some function were rewritten in C++ for greater efficiency
# and are located in src/LocationDepth.cpp

# Functions names from lsd were translated to harmonize name
# conventions in depthproc (lowerCamelCase for functions)
## - sample.max.depth -> lsdSampleDepth (CPP)
## - sample.depth.contours -> lsdSampleDepthContours (R)


#' @title Location-Scale depth class
#' @export
setClass("LSDepth", slots = c(max_depth = "numeric",
                              mu = "numeric",
                              sigma = "numeric"))


#' @title Location-Scale depth contour class
#' @export
setClass("LSDepthContour", slots = c(cont_depth = "numeric", sample = "numeric"), contains = "list")


#' @title Calculates the maximum sample location-scale depth
#' 
#' @description Calculates the maximum Student depth estimator of location and scale for one dimensional data (an alternative for MED and MAD or for the mean and standard deviation).
#' 
#' @param x one dimensional vector with sample
#' @param iter maximum number of iterations in algorith for calculation Location-Scale Depth
#' @param eps tolerance level
#' @param p_length is the maximum length of the precision step at the end
#' 
#' @details
#' Calculations are based on lsdepth algorithm written by Ch. Muller.
#' 
#' @export
#' 
#' @references
#' 
#' Mizera, I., Muller, C. H., 2004. Location-scale depth (with discussion). Journal of the American Statistical Association 99, 949 - 966.
#' 
#' 
#' @examples
#' x = rnorm(100)
#' lsdSampleMaxDepth(x)
#' y = rf(100, 4,10)
#' lsdSampleMaxDepth(y)
lsdSampleMaxDepth = function(x,iter=100,eps=0.0001,p_length=10)
{
  res = sampleMaxLocScaleDepthCPP(ry=as.numeric(x),iter=iter, eps=eps, p_length)
  res = as.numeric(res)
  
  new("LSDepth", max_depth = res[1], mu = res[2], sigma = res[3])
  #names(res)<-c("max.depth","mu","sigma")
  #return(res)
}


#' @title Calculate sample Mizera and Muller Student depth contours
#' 
#' @param x one dimensional vector with sample
#' @param depth depth level for contours
#' @param lengthmu number of points to evalute depth
#' 
#' @export
#' 
#' @details
#' Calculations are based on lsdepth algorithm written by Ch. Muller.
#' 
#' @references
#' 
#' Mizera, I., Muller, C. H., 2004. Location-scale depth (with discussion).  Journal of the American Statistical Association 99, 949-966.
#' 
#' @examples
#' ## EXAMPLE 1 for F-distribution
#' dcont = lsdSampleDepthContours(rf(200,4,7))
#' plot(dcont)
#' 
#' ## EXAMPLE 2 for normal distribution
#' ## - more contours calculated 
#' dcont_norm = lsdSampleDepthContours(rnorm(100),seq(0.05,0.4,0.05))
#' plot(dcont_norm)

lsdSampleDepthContours = function(x, depth = c(0.1,0.2,0.3,0.4), lengthmu=1000)
{
  depth = round(depth*length(x))
  x = sort(x)
  n = length(x)  
  
  dlen = length(depth)
  cont.all = vector("list",length=length(depth))
  
  for(i in 1:dlen)
  {
    d = depth[i]
    mu = getMuLS(x,n,d,lengthmu)
    cont = t(sapply(mu, function(mu) as.numeric(sampleDepthContForMuCPP(d,mu[1],x))))
    colnames(cont) = c("lbound","ubound","tbound","case","M")
    
    if(sum(cont[,"tbound"])>1){
      cont.exist = T
      tbound = as.logical(cont[,"tbound"])
      cont = cont[tbound,]
      lbound = cont[,"lbound"]
      ubound = cont[,"ubound"]
      mubound = mu[tbound]
      cont.all[[i]] = list(depth=d,cont.exist=cont.exist,mubound=mubound,
                           lbound=lbound,ubound=ubound)
    }
    else{
      cont.exist = F
      cat(" No contour for the depth ", d, "\n")
      cont.all[[i]] = list(depth=d,cont.exist=cont.exist)
    }
    
  }
  cont.all
  new("LSDepthContour", cont.all, cont_depth = depth/length(x), sample = x)
}


######################## Utils functions - not exported #####################
getMuLS = function(x,n,d,lengthmu)
{
  if(d>1){
    mu = seq(x[1+d-1],x[n-d+1],length=(lengthmu+1))
  }
  else{
    mu = seq(x[1],x[n],length=(lengthmu+1))
  }
  mu
}

#' @title Get location-scale contour from LSDepthContour object
#' @docType methods
#' @rdname lsdGetContour-methods
#' @export
#' 
#' @param x object of class LSDepthContour
#' @param cont single numeric - depth of contour to return
#' 
#' @details
#' Calculations are based on lsdepth algorithm written by Ch. Muller.
#'  
#' @examples
#' dcont = lsdSampleDepthContours(rf(200,4,7), depth = c(0.1,0.2))
#' 
#' #get contour that is present in dcont object
#' lsdGetContour(dcont,0.1)
#' 
#' # get contour that is not present in dcont
#' # it will be automatically calculated
#' lsdGetContour(dcont,0.3)
setGeneric("lsdGetContour", function(x, cont) standardGeneric("lsdGetContour"))

#' @rdname lsdGetContour-methods
#' @aliases lsdGetContour,LSDepthContour
#' @export
setMethod("lsdGetContour", signature = "LSDepthContour",function(x, cont)
{
  i = which(x@cont_depth == cont)
  if(length(i) > 0) return(x@.Data[[i]])
  lsdSampleDepthContours(x@sample, depth = cont)[[1]] 
})

#' @title Adds location scale depth contour to a plot.
#' @docType methods
#' @rdname lsdAddContour-methods
#' @export
#' 
#' @param x object of class LSDepthContour
#' @param cont depth of contour to plot
#' @param ... other arguments passed to polygon function
#' 
#' @examples
#' smp = rf(100,5,10)
#' x = lsdSampleDepthContours(smp)
#' plot(x)
#' lsdAddContour(x,0.1, col = "grey50")
#' lsdAddContour(x,0.3, col = "grey10", border = "red", lwd = 4)
setGeneric("lsdAddContour", function(x, cont = NULL,...) standardGeneric("lsdAddContour"))
#' @rdname lsdAddContour-methods
#' @aliases lsdAddContour,LSDepthContour
#' @export
setMethod("lsdAddContour", signature = c(x = "LSDepthContour"), function(x, cont = NULL, ...) 
{
  contour = lsdGetContour(x, cont)
  
  if(!contour$cont.exist) cat(" No contour for the depth ", cont, "\n")

  mubound = contour$mubound
  lbound  = contour$lbound
  ubound  = contour$ubound
  d   = contour$d
  lmu = length(mubound)
    
  polygon(c(mubound, rev(mubound)), c(ubound, rev(lbound)), ...)
})

#' @title Plot Location-Scale depth contours.
#' @export
#' 
#' @param x object of class LSDepthContour
#' @param cont plotted contours. Default NULL means that all contours stored in x will be plotted
#' @param ratio ratio
#' @param mu_min mu_min
#' @param mu_max mu_max
#' @param col vectors with area colors passed to polygon function
#' @param border vector with colors for borders
#' @param ... other parameters passed to polygon
#' 
#' @examples
#' 
#' smp = rf(100,5,10)
#' x = lsdSampleDepthContours(smp)
#' plot(x, col = paste0("grey", col = rev(seq(10,40,10))))
#' 
setMethod("plot", signature = c(x = "LSDepthContour"), function(x, cont = NULL, ratio=1,mu_min=NULL,mu_max=NULL, col = NULL, border = NULL,...)
{
  
  if(is.null(cont))   cont   = x@cont_depth
  ## numbers of conturs
  k = length(cont)
  cont = sort(cont)
  
  tmp_cont = lsdGetContour(x,cont[1])
  
  mubound = tmp_cont$mubound
  lbound  = tmp_cont$lbound
  ubound  = tmp_cont$ubound
  
  if(is.null(mu_min)) mu_min = mubound[1]
  if(is.null(mu_max)) mu_max = tail(mubound,1)
  
  if(is.null(col)) col = 0
  if(is.null(border)) border = 1
  if(length(col) != k) col = rep(col,k)
  if(length(border) != k) border = rep(border,k)
  
  plot(mubound,ubound,type="n", ylim=c(0,(1/ratio)*(mu_max-mu_min)),
         xlim=c(mu_min,mu_max), ylab=expression(sigma),xlab=expression(mu))
  sapply(1:length(cont), function(i,...) lsdAddContour(x,cont[i], col = col[i], border = border[i],...), ...)

  
  return(invisible())
})

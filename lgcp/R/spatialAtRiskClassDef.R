###
# spatialAtRisk class definition and functions
###

##' spatialAtRisk function
##' 
##' The methods for this generic function:\link{spatialAtRisk.default}, \link{spatialAtRisk.fromXYZ}, \link{spatialAtRisk.im}, \link{spatialAtRisk.function}, 
##' \link{spatialAtRisk.SpatialGridDataFrame}, \link{spatialAtRisk.SpatialPolygonsDataFrame} and \link{spatialAtRisk.bivden} are used to represent the fixed 
##' spatial component, lambda(s) in the log-Gaussian Cox process model. Typically lambda(s) would be represented as a spatstat object of class im, that encodes 
##' population density information. However, regardless of the physical interpretation of lambda(s), 
##' in lgcp we assume that it integrates to 1 over the observation window. The above methods make sure this condition is satisfied (with the 
##' exception of the method for objects of class function), as well as providing a framework for manipulating these structures. lgcp uses bilinear interpolation 
##' to project a user supplied lambda(s) onto a discrete grid ready for inference via MCMC, this grid can be obtained via the \link{selectObsWindow} function.
##'
##' Generic function used in the construction of spatialAtRisk objects. The class of spatialAtRisk objects provide a framework for describing
##' the spatial inhomogeneity of the at-risk population, lambda(s). This is in contrast to the class of 
##' \link{temporalAtRisk} objects, which describe the global levels of the population at risk, mu(t).
##'
##' Unless the user has specified lambda(s) directly by an R function (a mapping the from the real plane onto the non-negative real numbers, 
##' see ?spatialAtRisk.function), then it is only necessary to describe the population at risk up to a constant of proportionality, as the routines 
##' automatically normalise the lambda provided to integrate to 1.
##'
##' For reference purposes, the following is a mathematical description of a log-Gaussian Cox Process, it is best viewed in the pdf version of the manual.
##'
##' Let \eqn{\mathcal Y(s,t)}{\mathcal Y(s,t)} be a spatiotemporal Gaussian process, \eqn{W\subset R^2}{W\subset R^2} be an 
##' observation window in space and \eqn{T\subset R_{\geq 0}}{T\subset R_{\geq 0}} be an interval of time of interest. 
##' Cases occur at spatio-temporal positions \eqn{(x,t) \in W \times T}{(x,t) \in W \times T} 
##'  according to an inhomogeneous spatio-temporal Cox process,
##' i.e. a Poisson process with a stochastic intensity \eqn{R(x,t)}{R(x,t)},
##'   The number of cases, \eqn{X_{S,[t_1,t_2]}}{X_{S,[t_1,t_2]}}, arising in 
##'   any \eqn{S \subseteq W}{S \subseteq W} during the interval \eqn{[t_1,t_2]\subseteq T}{[t_1,t_2]\subseteq T} is 
##'   then Poisson distributed conditional on \eqn{R(\cdot)}{R(\cdot)},
##' \deqn{X_{S,[t_1,t_2]} \sim \mbox{Poisson}\left\{\int_S\int_{t_1}^{t_2} R(s,t)d sd t\right\}}{X_{S,[t_1,t_2]} \sim \mbox{Poisson}\left\{\int_S\int_{t_1}^{t_2} R(s,t)d sd t\right\}.}
##' Following Brix and Diggle (2001) and Diggle et al (2005), the intensity is decomposed multiplicatively as
##' \deqn{R(s,t) = \lambda(s)\mu(t)\exp\{\mathcal Y(s,t)\}.}{R(s,t) = \lambda(s)\mu(t)Exp\{\mathcal Y(s,t)\}.}
##' In the above, the fixed spatial component, \eqn{\lambda:R^2\mapsto R_{\geq 0}}{\lambda:R^2\mapsto R_{\geq 0}}, 
##' is a known function, proportional to the population at risk at each point in space and scaled so that
##' \deqn{\int_W\lambda(s)d s=1,}{\int_W\lambda(s)d s=1,}
##' whilst the fixed temporal component, 
##'  \eqn{\mu:R_{\geq 0}\mapsto R_{\geq 0}}{\mu:R_{\geq 0}\mapsto R_{\geq 0}}, is also a known function with
##' \deqn{\mu(t) \delta t = E[X_{W,\delta t}],}{\mu(t) \delta t = E[X_{W,\delta t}],}
##' for \eqn{t}{t} in a small interval of time, \eqn{\delta t}{\delta t}, over which the rate of the process over \eqn{W}{W} can be considered constant. 
##'
##' @param X an object
##' @param ... additional arguments
##' @return method spatialAtRisk
##' \enumerate{
##'     \item Brix A, Diggle PJ (2001). Spatiotemporal Prediction for log-Gaussian Cox processes. Journal of the Royal Statistical Society, Series B, 63(4), 823-841.
##'     \item Diggle P, Rowlingson B, Su T (2005). Point Process Methodology for On-line Spatio-temporal Disease Surveillance. Environmetrics, 16(5), 423-434.
##' }
##' @seealso \link{selectObsWindow} \link{lgcpPredict}, link{lgcpSim}, \link{spatialAtRisk.default}, \link{spatialAtRisk.fromXYZ}, \link{spatialAtRisk.im}, \link{spatialAtRisk.function}, \link{spatialAtRisk.SpatialGridDataFrame}, \link{spatialAtRisk.SpatialPolygonsDataFrame}, \link{spatialAtRisk.bivden} 
##' @export

spatialAtRisk <- function(X,...){
    UseMethod("spatialAtRisk")
}



##' spatialAtRisk.default function
##'
##' The default method for creating a spatialAtRisk object, which attempts to extract x, y and Zm values from the object using \code{xvals},
##' \code{yvals} and \code{zvals}.
##'
##' @method spatialAtRisk default
##' @param X an object
##' @param ... additional arguments  
##' @return object of class spatialAtRisk
##' \enumerate{
##'     \item Brix A, Diggle PJ (2001). Spatiotemporal Prediction for log-Gaussian Cox processes. Journal of the Royal Statistical Society, Series B, 63(4), 823-841.
##'     \item Diggle P, Rowlingson B, Su T (2005). Point Process Methodology for On-line Spatio-temporal Disease Surveillance. Environmetrics, 16(5), 423-434.
##' }
##' @seealso \link{lgcpPredict}, link{lgcpSim}, \link{spatialAtRisk.fromXYZ}, \link{spatialAtRisk.im}, \link{spatialAtRisk.function}, \link{spatialAtRisk.SpatialGridDataFrame}, \link{spatialAtRisk.SpatialPolygonsDataFrame}, \link{spatialAtRisk.bivden}, \code{xvals}, \code{yvals}, \code{zvals}
##' @export 

spatialAtRisk.default <- function(X,...){
    return(spatialAtRisk.fromXYZ(X=xvals(X),Y=yvals(X),Zm=zvals(X)))
}



##' spatialAtRisk.fromXYZ function
##'
##' Creates a spatialAtRisk object from a list of X, Y, Zm giving respectively the x and y coordinates of the grid and the 'z' values ie
##' so that Zm[i,j] is proportional to the at-risk population at X[i], Y[j]. 
##'
##' @method spatialAtRisk fromXYZ
##' @param X vector of x-coordinates
##' @param Y vector of y-coordinates
##' @param Zm matrix such that Zm[i,j] = f(x[i],y[j]) for some function f
##' @param ... additional arguments  
##' @return object of class spatialAtRisk
##' \enumerate{
##'     \item Brix A, Diggle PJ (2001). Spatiotemporal Prediction for log-Gaussian Cox processes. Journal of the Royal Statistical Society, Series B, 63(4), 823-841.
##'     \item Diggle P, Rowlingson B, Su T (2005). Point Process Methodology for On-line Spatio-temporal Disease Surveillance. Environmetrics, 16(5), 423-434.
##' }
##' @seealso \link{lgcpPredict}, link{lgcpSim}, \link{spatialAtRisk.default}, \link{spatialAtRisk.im}, \link{spatialAtRisk.function}, \link{spatialAtRisk.SpatialGridDataFrame}, \link{spatialAtRisk.SpatialPolygonsDataFrame}, \link{spatialAtRisk.bivden}
##' @export 

spatialAtRisk.fromXYZ <- function(X,Y,Zm,...){
    X <- sort(X)
    Y <- sort(Y)
    dx <- diff(X)
    dy <- diff(Y)
    if (!isTRUE(all.equal(rep(dx[1],length(dx)),dx)) | !isTRUE(all.equal(rep(dx[1],length(dx)),dx))){
        stop("Both X and Y must be equi-spaced vectors")
    }
    obj <- list()
    obj$X <- X
    obj$Y <- Y
    obj$Zm <- Zm
    NC <- dx[1]*dy[1]*sum(obj$Zm,na.rm=TRUE)
    obj$Zm <- obj$Zm / NC # normalise so lambda(s) integrates to 1
    if (any(obj$Zm<0,na.rm=TRUE)){
        stop("Negative at risk entries in matrix")
    }
    class(obj) <- c("fromXYZ","spatialAtRisk")
    attr(class(obj),"package") <- "lgcp"
    attr(obj,"NC") <- NC # normalising constant
    return(obj)
}

##' spatialAtRisk.im function
##'
##' Creates a spatialAtRisk object from a spatstat pixel image (im) object.
##'
##' @method spatialAtRisk im
##' @param X object of class im
##' @param ... additional arguments  
##' @return object of class spatialAtRisk
##' \enumerate{
##'     \item Brix A, Diggle PJ (2001). Spatiotemporal Prediction for log-Gaussian Cox processes. Journal of the Royal Statistical Society, Series B, 63(4), 823-841.
##'     \item Diggle P, Rowlingson B, Su T (2005). Point Process Methodology for On-line Spatio-temporal Disease Surveillance. Environmetrics, 16(5), 423-434.
##' }
##' @seealso \link{lgcpPredict}, link{lgcpSim}, \link{spatialAtRisk.default}, \link{spatialAtRisk.fromXYZ}, \link{spatialAtRisk.function}, \link{spatialAtRisk.SpatialGridDataFrame}, \link{spatialAtRisk.SpatialPolygonsDataFrame}, \link{spatialAtRisk.bivden}
##' @export 

spatialAtRisk.im <- function(X,...){
    obj <- list()
    obj$X <- sort(X$xcol)
    obj$Y <- sort(X$yrow)
    dx <- diff(obj$X)
    dy <- diff(obj$Y)
    if (!isTRUE(all.equal(rep(dx[1],length(dx)),dx)) | !isTRUE(all.equal(rep(dx[1],length(dx)),dx))){
        stop("Both X and Y must be equi-spaced vectors")
    }
    obj$Zm <- t(X$v)
    if (any(obj$Zm<0,na.rm=TRUE)){
        stop("Negative at risk entries in matrix")
    }
    NC <- dx[1]*dy[1]*sum(obj$Zm,na.rm=TRUE)
    obj$Zm <- obj$Zm / NC # normalise so lambda(s) integrates to 1
    class(obj) <- c("fromXYZ","spatialAtRisk")
    attr(class(obj),"package") <- "lgcp"
    attr(obj,"NC") <- NC # normalising constant
    return(obj)
}



##' spatialAtRisk.function function
##'
##' Creates a spatialAtRisk object from a function mapping R^2 onto the non negative reals. Note that for spatialAtRisk
##' objects defined in this manner, the user is responsible for ensurng that the integral of the function is 1 over the observation
##' window of interest.
##'
##' @method spatialAtRisk function
##' @param X a function with accepts arguments x and y that returns the at risk population at coordinate (x,y), which should be a numeric of length 1
##' @param warn whether to issue a warning or not
##' @param ... additional arguments
##' @return object of class spatialAtRisk
##' NOTE The function provided is assumed to integrate to 1 over the observation window, the user is responsible for ensuring this is the case.
##' \enumerate{
##'     \item Brix A, Diggle PJ (2001). Spatiotemporal Prediction for log-Gaussian Cox processes. Journal of the Royal Statistical Society, Series B, 63(4), 823-841.
##'     \item Diggle P, Rowlingson B, Su T (2005). Point Process Methodology for On-line Spatio-temporal Disease Surveillance. Environmetrics, 16(5), 423-434.
##' }
##' @seealso \link{lgcpPredict}, link{lgcpSim}, \link{spatialAtRisk.default}, \link{spatialAtRisk.fromXYZ}, \link{spatialAtRisk.im}, \link{spatialAtRisk.SpatialGridDataFrame}, \link{spatialAtRisk.SpatialPolygonsDataFrame}, \link{spatialAtRisk.bivden}
##' @export 

spatialAtRisk.function <- function(X,warn=TRUE,...){
    if(warn){
        warning("The function provided is assumed to integrate to 1 over the observation window.")
    }
    obj <- list()
    obj$f <- X
    class(obj) <- c("fromFunction","spatialAtRisk")
    attr(class(obj),"package") <- "lgcp"
    return(obj)
}



##' spatialAtRisk.SpatialGridDataFrame function
##'
##' Creates a spatialAtRisk object from an sp SpatialGridDataFrame object
##'
##' @method spatialAtRisk SpatialGridDataFrame
##' @param X a SpatialGridDataFrame object
##' @param ... additional arguments
##' @return object of class spatialAtRisk
##' \enumerate{
##'     \item Brix A, Diggle PJ (2001). Spatiotemporal Prediction for log-Gaussian Cox processes. Journal of the Royal Statistical Society, Series B, 63(4), 823-841.
##'     \item Diggle P, Rowlingson B, Su T (2005). Point Process Methodology for On-line Spatio-temporal Disease Surveillance. Environmetrics, 16(5), 423-434.
##' }
##' @seealso \link{lgcpPredict}, link{lgcpSim}, \link{spatialAtRisk.default}, \link{spatialAtRisk.fromXYZ}, \link{spatialAtRisk.im}, \link{spatialAtRisk.function}, \link{spatialAtRisk.SpatialPolygonsDataFrame}, \link{spatialAtRisk.bivden}
##' @export 

spatialAtRisk.SpatialGridDataFrame <- function(X,...){
    return(spatialAtRisk(X=xvals(X),Y=yvals(X),Zm=zvals(X)))
}



##' spatialAtRisk.SpatialPolygonsDataFrame function
##'
##' Creates a spatialAtRisk object from a SpatialPolygonsDataFrame object.
##'
##' @method spatialAtRisk SpatialPolygonsDataFrame
##' @param X a SpatialPolygonsDataFrame object; one column of the data frame should have name "atrisk", containing the aggregate population at risk for that region
##' @param ... additional arguments
##' @return object of class spatialAtRisk
##' \enumerate{
##'     \item Brix A, Diggle PJ (2001). Spatiotemporal Prediction for log-Gaussian Cox processes. Journal of the Royal Statistical Society, Series B, 63(4), 823-841.
##'     \item Diggle P, Rowlingson B, Su T (2005). Point Process Methodology for On-line Spatio-temporal Disease Surveillance. Environmetrics, 16(5), 423-434.
##' }
##' @seealso \link{lgcpPredict}, link{lgcpSim}, \link{spatialAtRisk.default}, \link{spatialAtRisk.fromXYZ}, \link{spatialAtRisk.im}, \link{spatialAtRisk.function}, \link{spatialAtRisk.SpatialGridDataFrame}, \link{spatialAtRisk.bivden}
##' @export 

spatialAtRisk.SpatialPolygonsDataFrame <- function(X,...){
    if(!any(names(data.frame(X))=="atrisk")){
        stop("One of the columns of the data frame must be be named 'atrisk', see ?spatialAtRisk.SpatialPolygonsDataFrame")
    }
    obj <- list()
    region.areas <- sapply(X@polygons,function(x){x@area})
    NC <- region.areas*sum(data.frame(X)$atrisk)
    obj$atrisk <- (data.frame(X)$atrisk/NC)  # normalise so lambda(s) integrates to 1
    X$atrisk <- obj$atrisk 
    obj$spdf <- X
    class(obj) <- c("fromSPDF","spatialAtRisk")
    attr(obj,"NC") <- NC # normalising constant
    return(obj)
}



##' spatialAtRisk.bivden function
##'
##' Creates a spatialAtRisk object from a sparr bivden object
##'
##' @method spatialAtRisk bivden
##' @param X a bivden object
##' @param ... additional arguments
##' @return object of class spatialAtRisk
##' \enumerate{
##'     \item Brix A, Diggle PJ (2001). Spatiotemporal Prediction for log-Gaussian Cox processes. Journal of the Royal Statistical Society, Series B, 63(4), 823-841.
##'     \item Diggle P, Rowlingson B, Su T (2005). Point Process Methodology for On-line Spatio-temporal Disease Surveillance. Environmetrics, 16(5), 423-434.
##' }
##' @seealso \link{lgcpPredict}, link{lgcpSim}, \link{spatialAtRisk.default}, \link{spatialAtRisk.fromXYZ}, \link{spatialAtRisk.im}, \link{spatialAtRisk.function}, \link{spatialAtRisk.SpatialGridDataFrame}, \link{spatialAtRisk.SpatialPolygonsDataFrame}
##' @export 

spatialAtRisk.bivden <- function(X,...){
    return(spatialAtRisk.fromXYZ(X=xvals(X),Y=yvals(X),Zm=zvals(X)))
}

##' spatialAtRisk.lgcpgrid function
##'
##' Creates a spatialAtRisk object from an lgcpgrid object
##'
##' @method spatialAtRisk lgcpgrid
##' @param X an lgcpgrid object
##' @param idx in the case that X$grid is a list of length > 1, this argument specifies which element of the list to convert. By default, it is the last.
##' @param ... additional arguments
##' @return object of class spatialAtRisk
##' \enumerate{
##'     \item Brix A, Diggle PJ (2001). Spatiotemporal Prediction for log-Gaussian Cox processes. Journal of the Royal Statistical Society, Series B, 63(4), 823-841.
##'     \item Diggle P, Rowlingson B, Su T (2005). Point Process Methodology for On-line Spatio-temporal Disease Surveillance. Environmetrics, 16(5), 423-434.
##' }
##' @seealso \link{lgcpPredict}, link{lgcpSim}, \link{spatialAtRisk.default}, \link{spatialAtRisk.fromXYZ}, \link{spatialAtRisk.im}, \link{spatialAtRisk.function}, \link{spatialAtRisk.SpatialGridDataFrame}, \link{spatialAtRisk.SpatialPolygonsDataFrame}
##' @export 

spatialAtRisk.lgcpgrid <- function(X,idx=length(X$grid),...){
    return(spatialAtRisk.fromXYZ(X=X$xvals,Y=X$yvals,Zm=X$grid[[idx]]))
}



###
# Conversion functions
###



##' as.im.fromXYZ function
##'
##' Convert an object of class fromXYZ (created by spatialAtRisk for example) into a spatstat im object.
##'
##' @method as.im fromXYZ
##' @param X object of class fromXYZ
##' @param ... additional arguments
##' @return object of class im containing normalised intensities
##' @seealso \link{as.im.fromSPDF}, \link{as.im.fromFunction}, \link{as.fromXYZ}
##' @export

as.im.fromXYZ <- function(X,...){
    mat <- zvals(X)
    return(im(mat=t(mat),xcol=xvals(X),yrow=yvals(X)))
}



##' as.im.fromSPDF function
##'
##' Convert an object of class fromSPDF (created by spatialAtRisk for example) into a spatstat im object.
##'
##' @method as.im fromSPDF
##' @param X an object of class fromSPDF
##' @param ncells number of cells to divide range into; default 100
##' @param ... additional arguments
##' @return object of class im containing normalised intensities
##' @seealso \link{as.im.fromXYZ}, \link{as.im.fromSPDF}, \link{as.im.fromFunction}, \link{as.fromXYZ}
##' @export

as.im.fromSPDF <- function(X,ncells=100,...){
    bbox <- summary(X$spdf)$bbox
    xgrid <- seq(bbox[1,1],bbox[1,2],length.out=ncells)
    ygrid <- seq(bbox[2,1],bbox[2,2],length.out=ncells)
    xyvals <- SpatialPoints(matrix(c(rep(xgrid,ncells),rep(ygrid,each=ncells)),ncells*ncells,2))
    #EJP: interp <- t(matrix(overlay(X$spdf,xyvals)$atrisk,ncells,ncells))
    interp <- t(matrix(over(xyvals, X$spdf)$atrisk,ncells,ncells))
    return(im(mat=interp,xcol=xgrid,yrow=ygrid)) 
}



##' as.im.fromFunction function
##'
##' Convert an object of class fromFunction(created by spatialAtRisk for example) into a spatstat im object.
##'
##' @method as.im fromFunction
##' @param X an object of class fromSPDF
##' @param xyt and objects of class stppp
##' @param M number of cells in x direction
##' @param N number of cells in y direction
##' @param ... additional arguments
##' @return object of class im containing normalised intensities
##' @seealso \link{as.im.fromXYZ}, \link{as.im.fromSPDF}, \link{as.im.fromFunction}, \link{as.fromXYZ}
##' @export

as.im.fromFunction <- function(X,xyt,M=100,N=100,...){
    del1 <- (xyt$window$xrange[2]-xyt$window$xrange[1])/M
	del2 <- (xyt$window$yrange[2]-xyt$window$yrange[1])/N
	xgrid <- xyt$window$xrange[1]+.5*del1+(0:(M-1))*del1
	ygrid <- xyt$window$yrange[1]+.5*del2+(0:(N-1))*del2
    xyvals <- cbind(rep(xgrid,N),rep(ygrid,each=M))
    matvals <- t(matrix(apply(xyvals,1,function(xy){X$f(x=xy[1],y=xy[2])}),M,N))
    return(im(mat=matvals,xcol=xgrid,yrow=ygrid)) 
}



##' as.fromXYZ function
##'
##' Generic function for conversion to a fromXYZ object (eg as would have been produced by spatialAtRisk for example.)
##'
##' @param X an object
##' @param ... additional arguments
##' @return generic function returning method as.fromXYZ
##' @seealso \link{as.im.fromXYZ}, \link{as.im.fromSPDF}, \link{as.im.fromFunction}, \link{as.fromXYZ}
##' @export

as.fromXYZ <- function(X,...){
    UseMethod("as.fromXYZ")
}



##' as.fromXYZ.fromFunction function
##'
##' Method for converting from the fromFunction class of objects to the fromXYZ class of objects. Clearly this requires the
##' user to specify a grid onto which to compute the discretised verion. 
##'
##' @method as.fromXYZ fromFunction
##' @param X an object of class fromFunction
##' @param xyt and objects of class stppp
##' @param M number of cells in x direction
##' @param N number of cells in y direction
##' @param ... additional arguments
##' @return object of class im containing normalised intensities
##' @seealso \link{as.im.fromXYZ}, \link{as.im.fromSPDF}, \link{as.im.fromFunction}, \link{as.fromXYZ}
##' @export

as.fromXYZ.fromFunction <- function(X,xyt,M=100,N=100,...){
    del1 <- (xyt$window$xrange[2]-xyt$window$xrange[1])/M
	del2 <- (xyt$window$yrange[2]-xyt$window$yrange[1])/N
	xgrid <- xyt$window$xrange[1]+.5*del1+(0:(M-1))*del1
	ygrid <- xyt$window$yrange[1]+.5*del2+(0:(N-1))*del2
    xyvals <- cbind(rep(xgrid,N),rep(ygrid,each=M))
    matvals <- matrix(apply(xyvals,1,function(xy){X$f(x=xy[1],y=xy[2])}),M,N)
    return(spatialAtRisk(X=xgrid,Y=ygrid,Zm=matvals)) 
}



##' spatialIntensities function
##'
##' Generic method for extracting spatial intensities.
##'
##' @param X an object
##' @param ... additional arguments
##' @return method spatialintensities
##' @seealso \link{spatialIntensities.fromXYZ}, \link{spatialIntensities.fromSPDF}
##' @export

spatialIntensities <- function(X,...){
    UseMethod("spatialIntensities")
}



##' spatialIntensities.fromXYZ function
##'
##' Extract the spatial intensities from an object of class fromXYZ (as would have been created by spatialAtRisk for example).
##'
##' @method spatialIntensities fromXYZ
##' @param X object of class fromXYZ
##' @param xyt object of class stppp or a list object of numeric vectors with names $x, $y
##' @param ... additional arguments
##' @return normalised spatial intensities
##' @seealso \link{spatialIntensities}, \link{spatialIntensities.fromSPDF}
##' @export

spatialIntensities.fromXYZ <- function(X,xyt,...){ 
    sgdf <- as.SpatialGridDataFrame(X)    
    #EJP: return(overlay(sgdf,SpatialPoints(cbind(xyt$x,xyt$y)))$atrisk)
    return(over(SpatialPoints(cbind(xyt$x,xyt$y)), sgdf)$atrisk)
}



##' spatialIntensities.fromSPDF function
##'
##' Extract the spatial intensities from an object of class fromSPDF (as would have been created by spatialAtRisk.SpatialPolygonsDataFrame for example).
##'
##' @method spatialIntensities fromSPDF
##' @param X an object of class fromSPDF
##' @param xyt object of class stppp or a list object of numeric vectors with names $x, $y
##' @param ... additional arguments
##' @return normalised spatial intensities
##' @seealso \link{spatialIntensities}, \link{spatialIntensities.fromXYZ}
##' @export

spatialIntensities.fromSPDF <- function(X,xyt,...){
    #EJP: return(overlay(X$spdf,SpatialPoints(cbind(xyt$x,xyt$y)))$atrisk)
    return(over(SpatialPoints(cbind(xyt$x,xyt$y)), X$spdf)$atrisk)
}



###
# Printing/plotting functions
###



##' print.fromXYZ function
##'
##' Print method for objects of class fromXYZ.
##'
##' @method print fromXYZ
##' @param x an object of class spatialAtRisk
##' @param ... additional arguments  
##' @return prints the object
##' @export

print.fromXYZ <- function(x,...){
    cat("SpatialAtRisk object\n")
    suppressWarnings(cat(paste("   X range : [",min(x$X),",",max(x$X),"]\n",sep="")))
    suppressWarnings(cat(paste("   Y range : [",min(x$Y),",",max(x$Y),"]\n",sep="")))
    suppressWarnings(cat(paste("   dim     : ",nrow(x$Zm)," x ",ncol(x$Zm),"\n",sep="")))
}



##' print.fromFunction function
##'
##' Print method for objects of class fromFunction.
##'
##' @method print fromFunction
##' @param x an object of class spatialAtRisk
##' @param ... additional arguments  
##' @return prints the object
##' @export

print.fromFunction <- function(x,...){
    cat("SpatialAtRisk object\n")
    print(x$f)
}



##' print.fromSPDF function
##'
##' Print method for objects of class fromSPDF.
##'
##' @method print fromSPDF
##' @param x an object of class spatialAtRisk
##' @param ... additional arguments  
##' @return prints the object
##' @export

print.fromSPDF <- function(x,...){
    NextMethod("print",x)
}



##' plot.fromXYZ function
##'
##' Plot method for objects of class fromXYZ.
##'
##' @method plot fromXYZ
##' @param x object of class spatialAtRisk
##' @param ... additional arguments
##' @return an image plot
##' @export

plot.fromXYZ <- function(x,...){
    image(x$X,x$Y,x$Zm,...)
}



##' plot.fromSPDF function
##'
##' Plot method for objects of class fromSPDF.
##'
##' @method plot fromSPDF
##' @param x an object of class spatialAtRisk
##' @param ... additional arguments  
##' @return prints the object
##' @export

plot.fromSPDF <- function(x,...){
    plot(x$spdf,...)
}

## spplot.fromSPDF function
##
## @name spplot.fromSPDF
## @docType methods
## @rdname spplot.fromSPDF
## method for spplot
#suppressWarnings(setMethod("spplot","fromSPDF", function(obj,...){spplot(obj$spdf,...)}))


###
# functions to extract x values
###



##' xvals function
##'
##' Generic for extractign the 'x values' from an object.
##'
##' @param obj an object of class spatialAtRisk
##' @param ... additional arguments
##' @return the xvals method
##' @seealso \link{yvals}, \link{zvals}, \link{xvals.default}, \link{yvals.default}, \link{zvals.default}, \link{xvals.fromXYZ}, \link{yvals.fromXYZ}, \link{zvals.fromXYZ}, \link{xvals.SpatialGridDataFrame}, \link{yvals.SpatialGridDataFrame}, \link{zvals.SpatialGridDataFrame}
##' @export

xvals <- function(obj,...){
    UseMethod("xvals")
}



##' xvals.default function
##'
##' Default method for extracting 'x values' looks for $X, $x in that order.
##'
##' @method xvals default
##' @param obj an object
##' @return the x values
##' @param ... additional arguments
##' @seealso \link{xvals},  \link{yvals}, \link{zvals}, \link{yvals.default}, \link{zvals.default}, \link{xvals.fromXYZ}, \link{yvals.fromXYZ}, \link{zvals.fromXYZ}, \link{xvals.SpatialGridDataFrame}, \link{yvals.SpatialGridDataFrame}, \link{zvals.SpatialGridDataFrame}
##' @export

xvals.default <- function(obj,...){
    if (any(names(obj)=="X")){
        return(obj$X)
    }
    else if (any(names(obj)=="x")){
        return(obj$x)
    }
    else{
        stop("Cannot extract x values")
    }
}



##' xvals.fromXYZ function
##'
##' Method for extracting 'x values' from an object of class fromXYZ
##'
##' @method xvals fromXYZ
##' @param obj a spatialAtRisk object
##' @param ... additional arguments
##' @return the x values
##' @seealso \link{xvals},  \link{yvals}, \link{zvals}, \link{xvals.default}, \link{yvals.default}, \link{zvals.default}, \link{yvals.fromXYZ}, \link{zvals.fromXYZ}, \link{xvals.SpatialGridDataFrame}, \link{yvals.SpatialGridDataFrame}, \link{zvals.SpatialGridDataFrame}
##' @export

xvals.fromXYZ <- function(obj,...){
    return(xvals.default(obj))
}



##' xvals.SpatialGridDataFrame function
##'
##' Method for extracting 'x values' from an object of class spatialGridDataFrame
##'
##' @method xvals SpatialGridDataFrame
##' @param obj an object
##' @param ... additional arguments
##' @return the x values
##' @seealso \link{xvals},  \link{yvals}, \link{zvals}, \link{xvals.default}, \link{yvals.default}, \link{zvals.default}, \link{xvals.fromXYZ}, \link{yvals.fromXYZ}, \link{zvals.fromXYZ}, \link{yvals.SpatialGridDataFrame}, \link{zvals.SpatialGridDataFrame}
##' @export

xvals.SpatialGridDataFrame <- function(obj,...){
    return(obj@grid@cellcentre.offset[1] + obj@grid@cellsize[1]*(0:(obj@grid@cells.dim[1]-1)))
}



###
# functions to extract y values
###



##' yvals function
##'
##' Generic for extractign the 'y values' from an object.
##'
##' @param obj an object of class spatialAtRisk
##' @param ... additional arguments
##' @return the yvals method
##' @seealso \link{xvals}, \link{zvals}, \link{xvals.default}, \link{yvals.default}, \link{zvals.default}, \link{xvals.fromXYZ}, \link{yvals.fromXYZ}, \link{zvals.fromXYZ}, \link{xvals.SpatialGridDataFrame}, \link{yvals.SpatialGridDataFrame}, \link{zvals.SpatialGridDataFrame}
##' @export

yvals <- function(obj,...){
    UseMethod("yvals")
}



##' yvals.default function
##'
##' Default method for extracting 'y values' looks for $Y, $y in that order.
##'
##' @method yvals default
##' @param obj an object
##' @param ... additional arguments
##' @return the y values
##' @seealso \link{xvals},  \link{yvals}, \link{zvals}, \link{xvals.default}, \link{zvals.default}, \link{xvals.fromXYZ}, \link{yvals.fromXYZ}, \link{zvals.fromXYZ}, \link{xvals.SpatialGridDataFrame}, \link{yvals.SpatialGridDataFrame}, \link{zvals.SpatialGridDataFrame}
##' @export

yvals.default <- function(obj,...){
    if (any(names(obj)=="Y")){
        return(obj$Y)
    }
    else if (any(names(obj)=="y")){
        return(obj$y)
    }
    else{
        stop("Cannot extract y values")
    }
}



##' yvals.fromXYZ function
##'
##' Method for extracting 'y values' from an object of class fromXYZ
##'
##' @method yvals fromXYZ
##' @param obj a spatialAtRisk object
##' @param ... additional arguments
##' @return the y values
##' @seealso \link{xvals},  \link{yvals}, \link{zvals}, \link{xvals.default}, \link{yvals.default}, \link{zvals.default}, \link{xvals.fromXYZ}, \link{zvals.fromXYZ}, \link{xvals.SpatialGridDataFrame}, \link{yvals.SpatialGridDataFrame}, \link{zvals.SpatialGridDataFrame}
##' @export

yvals.fromXYZ <- function(obj,...){
    return(yvals.default(obj))
}



##' yvals.SpatialGridDataFrame function
##'
##' Method for extracting 'y values' from an object of class SpatialGridDataFrame
##'
##' @method yvals SpatialGridDataFrame
##' @param obj an object
##' @param ... additional arguments
##' @return the y values
##' @seealso \link{xvals},  \link{yvals}, \link{zvals}, \link{xvals.default}, \link{yvals.default}, \link{zvals.default}, \link{xvals.fromXYZ}, \link{yvals.fromXYZ}, \link{zvals.fromXYZ}, \link{xvals.SpatialGridDataFrame}, \link{zvals.SpatialGridDataFrame}
##' @export

yvals.SpatialGridDataFrame <- function(obj,...){
    return(obj@grid@cellcentre.offset[2] + obj@grid@cellsize[2]*(0:(obj@grid@cells.dim[2]-1)))
}



###
# functions to extract z values
###



##' zvals function
##'
##' Generic for extractign the 'z values' from an object.
##'
##' @param obj an object 
##' @param ... additional arguments
##' @return the zvals method
##' @seealso \link{xvals},  \link{yvals}, \link{xvals.default}, \link{yvals.default}, \link{zvals.default}, \link{xvals.fromXYZ}, \link{yvals.fromXYZ}, \link{zvals.fromXYZ}, \link{xvals.SpatialGridDataFrame}, \link{yvals.SpatialGridDataFrame}, \link{zvals.SpatialGridDataFrame}
##' @export

zvals <- function(obj,...){
    UseMethod("zvals")
}



##' zvals.default function
##'
##' Default method for extracting 'z values' looks for $Zm, $Z, $z in that order.
##'
##' @method zvals default
##' @param obj an object
##' @param ... additional arguments
##' @return the x values
##' @seealso \link{xvals},  \link{yvals}, \link{zvals}, \link{xvals.default}, \link{yvals.default}, \link{xvals.fromXYZ}, \link{yvals.fromXYZ}, \link{zvals.fromXYZ}, \link{xvals.SpatialGridDataFrame}, \link{yvals.SpatialGridDataFrame}, \link{zvals.SpatialGridDataFrame}
##' @export

zvals.default <- function(obj,...){
    if (any(names(obj)=="Zm")){
        return(obj$Zm)
    }
    if (any(names(obj)=="Z")){
        return(obj$Z)
    }
    else if (any(names(obj)=="z")){
        return(obj$z)
    }
    else{
        stop("Cannot extract z values")
    }
}



##' zvals.fromXYZ function
##'
##' Method for extracting 'z values' from an object of class fromXYZ
##'
##' @method zvals fromXYZ
##' @param obj a spatialAtRisk object
##' @param ... additional arguments
##' @return the z values
##' @seealso \link{xvals},  \link{yvals}, \link{zvals}, \link{xvals.default}, \link{yvals.default}, \link{zvals.default}, \link{xvals.fromXYZ}, \link{yvals.fromXYZ}, \link{xvals.SpatialGridDataFrame}, \link{yvals.SpatialGridDataFrame}, \link{zvals.SpatialGridDataFrame}
##' @export

zvals.fromXYZ <- function(obj,...){
    return(zvals.default(obj))
}



##' zvals.SpatialGridDataFrame function
##'
##' Method for extracting 'z values' from an object of class SpatialGridDataFrame
##'
##' @method zvals SpatialGridDataFrame
##' @param obj an object
##' @param ... additional arguments
##' @return the z values
##' @seealso \link{xvals},  \link{yvals}, \link{zvals}, \link{xvals.default}, \link{yvals.default}, \link{zvals.default}, \link{xvals.fromXYZ}, \link{yvals.fromXYZ}, \link{zvals.fromXYZ}, \link{xvals.SpatialGridDataFrame}, \link{yvals.SpatialGridDataFrame}
##' @export

zvals.SpatialGridDataFrame<- function(obj,...){
    return(matrix(obj@data[,1],length(xvals(obj)),length(yvals(obj)))[,length(xvals(obj)):1])
}



###
# functions to convert between data formats
###



##' as.SpatialGridDataFrame function
##'
##' Generic method for convertign to an object of class SpatialGridDataFrame.
##'
##' @param obj an object
##' @param ... additional arguments
##' @return method as.SpatialGridDataFrame
##' @seealso \link{as.SpatialGridDataFrame.fromXYZ}
##' @export

as.SpatialGridDataFrame <- function(obj,...){
    UseMethod("as.SpatialGridDataFrame")
}


##' as.SpatialGridDataFrame.fromXYZ function
##'
##' Method for converting objects of class fromXYZ into those of class SpatialGridDataFrame
##'
##' @method as.SpatialGridDataFrame fromXYZ
##' @param obj an object of class spatialAtRisk
##' @param ... additional arguments
##' @return an object of class SpatialGridDataFrame
##' @seealso \link{as.SpatialGridDataFrame}
##' @export



as.SpatialGridDataFrame.fromXYZ <- function(obj,...){   
    sg <- SpatialGrid(GridTopology(c(obj$X[1],obj$Y[1]),c(diff(obj$X[1:2]),diff(obj$Y[1:2])),c(length(obj$X),length(obj$Y))))
	sgdf <- SpatialGridDataFrame(grid=sg,data=data.frame(atrisk=as.vector(obj$Zm[,length(obj$Y):1])))
	return(sgdf)
}

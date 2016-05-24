.complement <- function(ivec, imax) {
  result <- rep(TRUE,imax)
  result[ivec] <- FALSE
  return(which(result))
}

# Implimentation of TB algorithm without 'forced' supply points

.tb <- function(d,guess,verbose=FALSE) {
  config <- guess
  n <- ncol(d)
  repeat {
    old.config <- config
    config <- .bestswap(d,config,.complement(config,n))
    if (verbose) {
      cat("Configuration: ",config)
      score <- .dtotal(d,config)
      cat("  Score:",score,"\n")
    }
    if (all(old.config==config)) break 
  }
  return(config)
}

# Implimentation of TB algorithm with 'forced' supply points

.tb2 <- function(d,guess,n_force, verbose=FALSE) {
  config <- guess
  n <- ncol(d)
  repeat {
    old.config <- config
    if (missing(n_force)) {
      config <- .bestswap(d,config,.complement(config,n))
    } else {
      config <- .bestswap2(d,config,.complement(config,n),n_force)
    }
    if (verbose) {
      cat("Configuration: ",config)
      score <- .dtotal(d,config)
      cat("  Score:",score,"\n")
    }
    if (all(old.config==config)) break 
  }
  return(config)
}

##' Teitz-Bart algorithm applied to a 'raw' distance matrix 
##'  
##' @param  d - A distance matrix (not necessarily Euclidean)
##' @param  guess - a guess at the set of p points constituting the \eqn{p}-median
##' @param  verbose - if TRUE print out each swap in the algorithm (default is FALSE)
##' @return Set of point indices for \eqn{p}-median (may be local optimum)
##' 
##' @examples 
##' x1 <- rnorm(100)
##' y1 <- rnorm(100)
##' d <- as.matrix(dist(cbind(x1,y1)))
##' tb.raw(d,c(1,2))
##' @export
##' 


tb.raw <- .tb

##' Euclidean distances from a Spatial* or Spatial*DataFrame object
##'  
##' @param  swdf1 - First Spatial*DataFrame object
##' @param  swdf2 - Second Spatial*DataFrame object (if omitted,  defaults to the same value as \code{swdf1})
##' @param  scale - allows re-scaling eg: value of 1000 means distances in km if coordinates of \code{swdf1}/\code{swdf2} in meters.
##' @return Distance matrix (if \code{swdf1} or \code{swdf2} not SpatialPoints*, distances are based on points obtained from \code{coordinates} function)
##' 
##' @examples 
##' data(meuse)
##' coordinates(meuse) <- ~x+y
##' euc.dists(meuse,scale=1000)
##'
##' @export
##' 


euc.dists <- function(swdf1,swdf2,scale) {
  if (missing(swdf2)) swdf2 <- swdf1
  xy1 <- coordinates(swdf1)
  xy2 <- coordinates(swdf2)
  if (!missing(scale)) {
    xy1 <- xy1/scale
    xy2 <- xy2/scale
  }
  return(.dmat(xy1[,1],xy2[,1],xy1[,2],xy2[,2])) 
}

##' Minkowski distances from a Spatial* or Spatial*DataFrame object
##'  
##' @param  swdf1 - First Spatial*DataFrame object
##' @param  swdf2 - Second Spatial*DataFrame object (if omitted,  defaults to the same value as \code{swdf1})
##' @param  pwr - Minkowski exponent
##' @param  scale - allows re-scaling eg: value of 1000 means distances in km if coordinates of \code{swdf1}/\code{swdf2} in meters.
##' @param  weight - weight for each element in swdf1 (the demand locations)
##' @return Distance matrix (if \code{swdf1} or \code{swdf2} not SpatialPoints*, 
##'   distances are based on points obtained from \code{coordinates} function)
##' 
##' @examples 
##' data(meuse)
##' coordinates(meuse) <- ~x+y
##' d1 <- mink.dists(meuse,pwr=1,scale=1000)   # Taxicab metric
##' d2 <- mink.dists(meuse,pwr=Inf,scale=1000) # Works for limiting case
##'
##' @export
##' 


mink.dists <- function(swdf1,swdf2,pwr,scale,weight) {
  if (missing(swdf2)) swdf2 <- swdf1
  xy1 <- coordinates(swdf1)
  xy2 <- coordinates(swdf2)
  if (!missing(scale)) {
    xy1 <- xy1/scale
    xy2 <- xy2/scale
  }
  temp <- .dmatex(xy1[,1],xy2[,1],xy1[,2],xy2[,2],pwr)
  if (missing(weight)) return(temp)
  sweep(temp, 1, weight, "/")
}



##' Teitz-Bart algorithm applied to Spatial* and Spatial*DataFrame objects 
##' 
##' This reports the \eqn{p}-median set
##'  
##' @param  swdf1 - first Spatial* or Spatial*DataFrame objects - the 'demand' set
##' @param  swdf2 - second Spatial* or Spatial*DataFrame objects - the 'supply' set (if omitted,  defaults to the same value as \code{swdf1})
##' @param  p - either a guess at the initial \eqn{p}-median set of a single integer indicating the size of the set (which is then chosen randomly)
##' @param  metric - the distance matrix (defaults to Euclidean computed via \code{euc.dists(swdf1,swdf2)} if not supplied)
##' @param  verbose - if TRUE print out each swap in the algorithm (default is FALSE)
##' @return Set of point indices for \eqn{p}-median (may be local optimum)
##' 
##' @examples 
##' data(meuse)
##' coordinates(meuse) <- ~x+y
##' tb(meuse,p=5)
##' @export
##' 

tb <- function(swdf1,swdf2,p,metric,verbose=FALSE) {
  if (missing(swdf2))  swdf2 <- swdf1
  n.choices <- nrow(coordinates(swdf2)) 
  if (length(p) == 1) p <- sample(n.choices,p)
  if (missing(metric)) metric <- euc.dists(swdf1,swdf2)
  result <- .tb(metric,p,verbose)
  return(result)
}

##' Creates the lines for a 'star diagram' 
##'  
##' @param  swdf1 - first Spatial* or Spatial*DataFrame objects
##' @param  swdf2 - second Spatial* or Spatial*DataFrame objects (if omitted,  defaults to the same value as \code{swdf1})
##' @param  alloc - a list saying which coordinate in swdf2 is allocated to each point in swdf1 (if ommitted, looks for \code{allocation} column in \code{swdf1})
##' 
##' @examples 
##' data(meuse)
##' coordinates(meuse) <- ~x+y
##' allocations.list <- allocate(meuse,p=5)
##' star.lines <- star.diagram(meuse,alloc=allocations.list)
##' plot(star.lines)
##' 
##' # Acquire allocations from swdf1
##' require(GISTools)
##' set.seed(461976) # Reproducibility
##' data(georgia)
##' georgia3 <- allocations(georgia2,p=8)
##' plot(georgia3,border='grey')
##' plot(star.diagram(georgia3),col='darkblue',lwd=2,add=TRUE)
##' 
##' 
##' @export
##' 


star.diagram <- function(swdf1,swdf2,alloc) {
  if (missing(swdf2))  swdf2 <- swdf1
  if (missing(alloc)) {
    alloc <- swdf1$allocation
  }
  co1 <- coordinates(swdf1)
  co2 <- coordinates(swdf2)
  result <- vector(nrow(co1),mode='list')
  for (i in 1:nrow(co1))  result[[i]] <- Lines(list(Line(cbind(c(co1[i,1],co2[alloc[i],1]),c(co1[i,2],co2[alloc[i],2])))),ID=sprintf("Star%d",i))
  sl <- SpatialLines(result)
  sldf <- SpatialLinesDataFrame(sl,data.frame(allocate=alloc),match.ID=FALSE)
  return(sldf)}




##' Teitz-Bart algorithm applied to Spatial* and Spatial*DataFrame objects 
##' 
##' Return demand Spatial*Dataframe with new columns giving allocation id and distance to supply point
##'  
##' @param  swdf1 - first Spatial* or Spatial*DataFrame objects
##' @param  swdf2 - second Spatial* or Spatial*DataFrame objects (if omitted,  defaults to the same value as \code{swdf1})
##' @param  force - list of supply points or logical vector with length the same as the number of supply points that are forced to be used - eg e
##' @param  p - either a guess at the initial \eqn{p}-median set of a single integer indicating the size of the set (which is then chosen randomly)
##' @param  metric - the distance matrix (defaults to Euclidean computed via \code{euc.dists(swdf1,swdf2)} if not supplied)
##' @param  verbose - if TRUE print out each swap in the algorithm (default is FALSE)
##' @return Copy of swdf1 with extra data columns called \code{allocation} and \code{allocdist}
##' with indices for each element from the \eqn{p}-median set
##' 
##' @examples 
##' 
##' require(RColorBrewer)
##' require(GISTools)
##' data(georgia)
##' georgia3 <- allocations(georgia2,p=5,force=c(1,120,44))
##' col.index <- match(georgia3$allocation,unique(georgia3$allocation))
##' col.alloc <- brewer.pal(5,'Accent')[col.index]
##' par(mfrow=c(1,2))
##' plot(georgia3,col=col.alloc)
##' choropleth(georgia3,georgia3$allocdist)
##' 
##' 
##' # Use in conjunction with rgeos
##' require(rgeos)
##' require(GISTools)
##' georgia3 <- allocations(georgia2,p=5,force=c(1,120,44))
##' georgia4 <- gUnaryUnion(georgia3,georgia3$allocation)
##' plot(georgia4)
##' plot(star.diagram(georgia3),col='darkred',lwd=2,add=TRUE)
##' 
##' @export
##' 

allocations <-function(swdf1,swdf2,force,p,metric,verbose=FALSE) {
  if (missing(swdf2))  swdf2 <- swdf1
  n.choices <- nrow(coordinates(swdf2)) 
  if (missing(metric)) metric <- euc.dists(swdf1,swdf2)
  if (missing(force)) { 
    nni <- allocate(swdf1,swdf2,p=p,metric=metric,verbose=verbose)
  } else {
    nni <- allocate(swdf1,swdf2,force=force,p=p,metric=metric,verbose=verbose)
  }
  n.demand <- nrow(coordinates(swdf1))
  picker <- cbind(1:n.demand,nni)
  dist <- metric[picker]
  swdf1_df <- data.frame(swdf1)
  swdf1_df <- cbind(swdf1_df,data.frame(allocation=nni,allocdist=dist))
  if (class(geometry(swdf1)) == "SpatialPolygons") return(SpatialPolygonsDataFrame(swdf1,swdf1_df))
  if (class(geometry(swdf1)) == "SpatialPoints") return(SpatialPointsDataFrame(swdf1,swdf1_df))
  return(SpatialLinesDataFrame(swdf1,swdf1_df))
}




##' Teitz-Bart algorithm applied to Spatial* and Spatial*DataFrame objects 
##' 
##' This function returns the allocations for each demand point - in terms of the index number of 
##' the record in \code{swdf2} assigned as the supply point.  This version is useful as part of
##' code inside other functions
##'  
##' @param  swdf1 - first Spatial* or Spatial*DataFrame objects
##' @param  swdf2 - second Spatial* or Spatial*DataFrame objects (if omitted,  defaults to the same value as \code{swdf1})
##' @param  force - list of supply points or logical vector with length the same as the number of supply points that are forced to be used - eg existing outlets
##' @param  p - either a guess at the initial \eqn{p}-median set of a single integer indicating the size of the set (which is then chosen randomly)
##' @param  metric - the distance matrix (defaults to Euclidean computed via \code{euc.dists(swdf1,swdf2)} if not supplied)
##' @param  verbose - if TRUE print out each swap in the algorithm (default is FALSE)
##' @return List of nearest neigbour indices for each element from the \eqn{p}-median set
##' 
##' @examples 
##' data(meuse)
##' coordinates(meuse) <- ~x+y
##' allocate(meuse,p=5)
##' 
##'
##' 
##' require(RColorBrewer)
##' require(GISTools)
##' data(georgia)
##' allocations.list <- allocate(georgia2,p=5)
##' zones <- gUnaryUnion(georgia2,allocations.list)
##' plot(zones,col=brewer.pal(5,"Accent"))
##' plot(georgia2,border=rgb(0,0,0,0.1),add=TRUE)
##' points(coordinates(georgia2)[allocations.list,],pch=16,cex=2,col=rgb(1,0.5,0.5,0.1))
##' 
##' @export
##' 

allocate <-function(swdf1,swdf2,force,p,metric,verbose=FALSE) {
  if (missing(swdf2))  swdf2 <- swdf1
  n.choices <- nrow(coordinates(swdf2)) 
  if (length(p) == 1) {
    n.subset <- p
    p <- sample(n.choices,p)
  } else {
    n.subset <- length(p)
  }
  if (missing(metric)) metric <- euc.dists(swdf1,swdf2)
  if (missing(force)) {
    indices <- .tb(metric,p,verbose)
  } else {
    if (is.logical(force)) force <- which(force)
    p <- c(force,setdiff(p,force))[1:n.subset]
    indices <- .tb2(metric,p,as.integer(length(force)), verbose)
  }
  nni <- .rviss(metric,indices)
  return(nni)
}


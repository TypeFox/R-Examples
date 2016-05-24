#' Summary of a prevR object.
#' 
#' Method \code{show} for objects of class \code{\link[=prevR-class]{prevR}}: 
#' shows a summary of the object's characteristics.
#' 
#' @param object object of class \code{\link[=prevR-class]{prevR}}.
#' 
#' @note Exactly the same as \code{\link{print,prevR-method}}.
#' @seealso \code{\link{summary,prevR-method}}.
#' @examples 
#' fdhs
#' \dontrun{
#'  dhs <- rings(fdhs,N=c(100,300,500))
#'  dhs
#' }
#'
#' @aliases show show-methods show,prevR-method

setMethod("show","prevR",
function(object){
  ###############################################################################################
  # Cette fonction affiche un resume des donnees et fourni les parametres disponibles de rings
  #  Exemple 
  #  Number of clusters : 466
  #  Number of valid observations : 466
  #  Are the date weighted ? TRUE
  #  National prevalence : 5.47
  #  Weighted national prevalence : 5.35
  #  Projection: +proj=longlat +ellps=WGS84
  #  
  #  Coordinate range
  #        min      max
  #  x 9.02489 15.44588
  #  y 2.26966 12.77069
  #  
  #  #######################################
  #  ######### rings slot ##################
  #  #######################################
  #  The available (N,R) couples are
  #      N   R
  #  1 250 Inf
  #  2 500 Inf
  #  The available variables are
  #  r.pos, r.n, r.prev, r.radius, r.clusters, r.wpos, r.wn, r.wprev
  ###############################################################################################
 
  clusters               = slot(object,"clusters")
  boundary               = slot(object,"boundary")

  clustersNumber         = nrow(clusters)
  ObservationNumber      = sum(clusters$n)
  PositiveCases          = sum(clusters$pos)
  isWeightedData         = !any(is.na(match(c("wn","wpos"),names(clusters))))
  nationalPrev           = 100*sum(clusters$pos,na.rm=T)/sum(clusters$n,na.rm=T)
  if(isWeightedData){
    weightedNationalPrev = 100*sum(clusters$wpos,na.rm=T)/sum(clusters$wn,na.rm=T)
  }
  proj                   = object@proj@projargs
  coordinatesRange       = rbind(range(clusters$x,na.rm=T),range(clusters$y,na.rm=T))
  dimnames(coordinatesRange) = list(c("x","y"),c("min","max"))
  boundaryCoordinatesRange = NULL
  if(attr(boundary,"valid")){
    boundaryCoordinatesRange = slot(boundary,"bbox")
  }
  
  message("Object of class 'prevR'\n",domain="R-prevR")
  
  message(gettextf("Number of clusters: %i",clustersNumber,domain="R-prevR"))
  message(gettextf("Number of observations: %i",ObservationNumber,domain="R-prevR"))
  message(gettextf("Number of positive cases: %i",PositiveCases,domain="R-prevR"))
  if(isWeightedData){
    message("The dataset is weighted.",domain="R-prevR")
  } else {
    message("The dataset is not weighted.",domain="R-prevR")
  }
  message(gettextf("\nNational prevalence: %.2f%%",nationalPrev,domain="R-prevR"))
  if(isWeightedData){
    message(gettextf("National weighted prevalence: %.2f%%",weightedNationalPrev,domain="R-prevR"))
  }
  message(gettextf("\nProjection used: %s",proj,domain="R-prevR"))
    message("\nCoordinate range",domain="R-prevR")
  print(coordinatesRange)
  if(attr(boundary,"valid")){
    message("\nBoundary coordinate range",domain="R-prevR") 
    print(boundaryCoordinatesRange)
  }
  if(is.prevR(object,"rings")){
    rings = slot(object,"rings")
    N = sapply(rings,function(x) x$N)
    R = sapply(rings,function(x) x$R)
    couples = cbind(N=N,R=R)
    dimnames(couples) = list(1:nrow(couples) ,c("N","R"))
    message("\nAvailable (N,R) couples in the slot 'rings':",domain="R-prevR")
    print(as.data.frame(couples),row.names=FALSE)
  }
}

)
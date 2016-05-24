# ---- roxygen documentation ----
#
#' @title Half-weight Association Index
#'
#' @description
#' This function computes the Half-weight Association Index for examining the presence
#' of dynamic interaction in wildlife telemetry studies. This implementation follows
#' that outlined in the paper Atwood and Weeks (2003).
#'
#' @details
#' This function can be used to test for the presence of dynamic interaction within 
#' the shared area (often termed the overlap zone) of the two animals home ranges. Specifically,
#' HAI is calculated in identical fashion to that for \code{Ca}, but considers only those fixes 
#' in the shared area. Typically, the overlap zone (OZ) is easily obtained by taking the spatial
#' intersection of two polygon home ranges.
#'
#' @param traj1 an object of the class \code{ltraj} which contains the time-stamped
#'    movement fixes of the first object. Note this object must be a \code{type II
#'    ltraj} object. For more information on objects of this type see \code{
#'    help(ltraj)}.
#' @param traj2 same as \code{traj1}.
#' @param OZ spatial polygon associated with the home range (or some other form of) spatial overlap between 
#'    \code{traj1} and \code{traj2}. Required to be an object that coerces to class \code{SpatialPolygons}.
#' @param tc time threshold for determining simultaneous fixes -- see function: \code{GetSimultaneous}.
#' @param dc distance tolerance limit (in appropriate units) for defining when 
#'         two fixes are spatially together.
#'
#' @return
#' This function returns the numeric value of the HAI statistic. Values near 1 indicate
#' attraction within the shared home range area, while values near 0 indicate avoidance
#' within this shared area.
#'
#' @references
#' Atwood, T.C. and Weeks Jr., H.P. (2003) Spatial home-range overlap and temporal
#' interaction in eastern coyotes: The influence of pair types and fragmentation.
#' \emph{Canadian Journal of Zoology}, \bold{81}: 1589-1597.\cr\cr
#'
#' @keywords indices
#' @seealso GetSimultaneous, Ca
#' @examples
#' \dontrun{
#' data(deer)
#' deer37 <- deer[1]
#' deer38 <- deer[2]
#' library(adehabitatHR)
#' library(sp)
#' library(rgeos)
#' #use minimum convex polygon for demonstration...
#' hr37 <- mcp(SpatialPoints(ld(deer37)[,1:2]))
#' hr38 <- mcp(SpatialPoints(ld(deer38)[,1:2]))
#' OZ <- gIntersection(hr37,hr38)
#' #tc = 7.5 minutes, dc = 50 meters
#' HAI(deer37, deer38, OZ=OZ, tc=7.5*60, dc=50)
#' }
#' @export
#
# ---- End of roxygen documentation ----
HAI <- function(traj1, traj2, OZ, tc = 0,dc = 50){
  trajs <- GetSimultaneous(traj1, traj2, tc)
  #convert ltraj objects to dataframes
  tr1 <- ld(trajs[1])
  tr2 <- ld(trajs[2])
  
  n <- nrow(tr1)

  #convert traj to spatial points
  pts1 <- SpatialPoints(tr1[,1:2])
  pts2 <- SpatialPoints(tr2[,1:2])
  #identify only those points in the overlap zone
  ind1 <- gIntersects(pts1,OZ,byid=T)
  ind2 <- gIntersects(pts2,OZ,byid=T)
  tr1 <- tr1[ind1,]
  tr2 <- tr2[ind2,]
  
  #Get only those fixes that are simultaneous
  tr1.t <- tr1[tr1$date %in% tr2$date,]
  tr2.t <- tr2[tr2$date %in% tr1$date,]
  
  euc <- function(x1,y1,x2,y2){sqrt(((x1 - x2)^2) + ((y1 - y2)^2))}
  
  #calculate the count of distances below threshold
  Do <- euc(tr1.t$x,tr1.t$y,tr2.t$x,tr2.t$y)
  n <- length(which(Do < dc))
  #Compute count of all other fixes that are not proximal
  x <- nrow(tr1) - n
  y <- nrow(tr2) - n
  
  hai <- n/(n + 0.5*(x + y))
  #print(paste("A critical distance of ",dc," was used, resulting in HAI = ", hai,".",sep=""))
  return(hai)
}

#======================= End of HAI function ======================================

# ---- roxygen documentation ----
#
#' @title Mapping wildlife contacts
#'
#' @description
#' The function \code{contacts} is a simple function for mapping where wildlife contacts
#' occur on the landscape with wildilfe telemetry data. 
#' 
#' @details
#' The function \code{contacts} can be used to map where contacts occur on the lansdcape, 
#' contacts being defined spatially based on a distance threshold \code{dc} and temporally
#' based on the time threshold \code{tc} -- see the function \code{getsimultaneous}. The location
#' of the contact is calculated as the midpoint between the two fixes that are determined
#' to be a "contact" based on \code{dc} and \code{tc}.  
#'
#' @param traj1 an object of the class \code{ltraj} which contains the time-stamped
#'    movement fixes of the first object. Note this object must be a \code{type II
#'    ltraj} object. For more information on objects of this type see \code{help(ltraj)}.
#' @param traj2 same as \code{traj1}.
#' @param tc time threshold for determining simultaneous fixes -- see function: \code{GetSimultaneous}.
#' @param dc distance tolerance limit (in appropriate units) for defining when 
#'    two fixes are spatially together.
#' @param proj4string a string object containing the projection information to be passed included in the output 
#'    \code{SpatialPolygonsDataFrame} object. For more information see the \code{CRS-class} in the packages
#'    \code{sp} and \code{rgdal}. Default is \code{NA}.
#'
#' @return
#' A \code{SpatialPointsDataFrame} containing the locations of the contacts. The time of the 
#' contact is stored in the attributes of the \code{SpatialPointsDataFrame} object, along with
#' the actual distance between fixes. 
#'
#' @keywords indices
#' @seealso GetSimultaneous, Prox
#' @examples
#' data(deer)
#' deer37 <- deer[1]
#' deer38 <- deer[2]
#' #tc = 7.5 minutes, dc = 50 meters
#' spts <- contacts(deer37, deer38, tc=7.5*60, dc=50)
#' 
#' @export
#

contacts <- function(traj1,traj2,dc=0,tc=0,proj4string=CRS(as.character(NA))){
  #convert ltraj objects to dataframes
  trajs <- GetSimultaneous(traj1, traj2, tc)
  #convert ltraj objects to dataframes
  tr1 <- ld(trajs[1])
  tr2 <- ld(trajs[2])
  trdf <- data.frame(x=(0.5*(tr1$x+tr2$x)),y=(0.5*(tr1$y+tr2$y)),date=tr1$date,prox=sqrt(((tr1$x - tr2$x)^2) + ((tr1$y - tr2$y)^2)))
  ind <- which(trdf$prox <= dc)
  spts <- SpatialPointsDataFrame(trdf[ind,1:2],data=trdf[ind,],proj4string=proj4string)
  return(spts)
}
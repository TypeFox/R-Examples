#' Plot arrows and segments showing the size and direction of currents.
#' 
#' Plot arrows and segments showing the size and direction of currents.
#' 
#' 
#' @param data Data in a list with components \code{lat} and \code{lon} with
#' decimal degrees, and \code{current} with the current magnitude.
#' @param maxsize Maximum current segment size.
#' @param maxn Current given with \code{maxsize}, defaults to
#' \code{max(data$current)}.
#' @param col Color of current arrows and segments.
#' @param lwd Line width of the segments showing current.
#' @param arrowsize Arrow size.
#' @param center Whether or not to center the arrow, defaults to \code{TRUE}.
#' @note Needs further checking and elaboration.
#' @keywords aplot
#' @export currentarrows
currentarrows <-
function(data,maxsize=0.5,maxn,col="blue",lwd=2,arrowsize=0.2,center=T){
  geopar <- getOption("geopar")
  res <- list()
  xsizerat <- geopar$gpar$pin[1]/diff(geopar$origin$lon)
  ysizerat <- geopar$gpar$pin[2]/diff(geopar$origin$lat)
  data$rat <- cos(data$lat*pi/180) 
  if(missing(maxn)) 
    maxn <- max(data$current)
  
  tmp <- data.frame(lat=c(1,1),lon=c(1,1)) 
  for(i in 1:nrow(data)) {
    tmp[1,] <- data[i,c("lat","lon")]
    tmp[2,"lon"] <- tmp[1,"lon"]+maxsize*data$current[i]/maxn*cos(data$angle[i]*pi/180)/data$rat[i]/ysizerat
    tmp[2,"lat"] <- tmp[1,"lat"]+maxsize*data$current[i]/maxn*sin(data$angle[i]*pi/180)/ysizerat
    if(center){ #center the arrow, else start
      dlat <- diff(tmp$lat)
      dlon <- diff(tmp$lon)
      tmp$lat <- tmp$lat -dlat/2
      tmp$lon <- tmp$lon -dlon/2
    }
    res[[i]] <- tmp
    SegmentWithArrow(tmp,lwd=lwd,size=arrowsize,col=col)
  }
  return(invisible(res)) 
}


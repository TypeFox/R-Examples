#' @title Build node data
#' @description Helper funtion to build the origin/destination node data structure. 
#' 
#' @param x              A data.frame contaning node (site) data
#' @param group.ids      Character vector of unique identifer that can be used to join to graph
#' @param from.parms     Character vector of independent "from" variables
#' @param to.parms       Character vector of independent "to" variables. If NULL is the same as from.parms
#' 
#' @return data.frame 
#'
#' @note Unless a different set of parameters will be used as the destination (to) there is no need to define the argument "to.parms" and the "from.parm" will be used to define both set of parameters.  
#' @note The resulting data.frame represents the origin (from) and destination (to) data structure for use in garavity model. This is node structure is also know in the gravity literature as producer (from) and attractor (to). 
#'
#' @author Jeffrey S. Evans  <jeffrey_evans@@tnc.org> and Melanie Murphy <melanie.murphy@@uwyo.edu>
#'
#' @examples 
#' data(ralu.site)
#'
#' # Build from/to site (node) level data structure 
#' site.parms = c("AREA_m2", "PERI_m", "Depth_m", "TDS")
#' site <- build.node.data(ralu.site@data, group.ids = c("SiteName"), 
#'                         from.parms = site.parms ) 
#'  
#' @export
build.node.data <- function(x, group.ids, from.parms, to.parms = NULL) {
    if(!inherits(x, "data.frame")) stop("x is not a data.frame")
	  for(i in from.parms) {
        if (is.na(charmatch(i, names(x))))
          stop(i, " is not a column in the data.frame")
      }
      if(is.null(to.parms)) to.parms = from.parms	
        id.col <- which( names(x) %in% group.ids ) 
	      from.col <- which( names(x) %in% from.parms )
	    to.col <- which( names(x) %in% to.parms )      
      from <- x[,c(id.col,from.col)]
        names(from) <- c(group.ids, paste("from", from.parms, sep="."))   
	  to <- x[,c(id.col,to.col)]
	    names(to) <- c(group.ids, paste("to", to.parms, sep=".")) 
    site <- data.frame(from, to)	
      site <- site[,-(dim(to)[2]+1)]
  return( site )	
}

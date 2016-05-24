#' @title Saturated or K Nearest Neighbour Graph
#' @description Creates a kNN or saturated graph SpatialLinesDataFrame object  
#'
#' @param x              sp SpatialPointsDataFrame object
#' @param row.names      Unique row.names assigned to results  
#' @param k              K nearest neighbours
#' @param max.dist       Maximum length of an edge (used for distance constraint)
#' @param sym            Create symmetrical graph (FALSE/TRUE)
#' @param drop.diag      Drop diag (duplicate edges) (FALSE/TRUE)
#' @param long.lat       Coordinates are longitude-latitude decimal degrees (FALSE/TRUE), in which case distances are measured in kilometers
#' 
#' @return   SpatialLinesDataFrame object with:
#'    i        Name of column in x with FROM (origin) index
#'    j        Name of column in x with TO (destination) index
#'    from_ID     Name of column in x with FROM (origin) region ID
#'    to_ID     Name of column in x with TO (destination) region ID
#'    length   Length of each edge (line) in projection units or kilometers if long.lat = TRUE
#'
#' @note ...
#'
#' @author Jeffrey S. Evans  <jeffrey_evans@@tnc.org> and Melanie Murphy <melanie.murphy@@uwyo.edu>
#'
#' @references
#' Murphy, M. A. & J.S. Evans. (in prep). "GenNetIt: gravity analysis in R for landscape genetics" 
#' Murphy M.A., R. Dezzani, D.S. Pilliod & A.S. Storfer (2010) Landscape genetics of high mountain frog metapopulations. Molecular Ecology 19(17):3634-3649 
#'
#' @examples
#'  library(sp) 
#'    data(ralu.site)
#'
#'  # Saturated spatial graph
#'  sat.graph <- knn.graph(ralu.site, row.names=ralu.site@@data[,"SiteName"])
#'	  head(sat.graph@data)
#'  
#'  # Distanced constrained spatial graph
#'  dist.graph <- knn.graph(ralu.site, row.names=ralu.site@@data[,"SiteName"], max.dist = 5000)
#'	
#'  par(mfrow=c(1,2))	
#'	plot(sat.graph, col="grey")
#'	  points(ralu.site, col="red", pch=20, cex=1.5)
#'      box()
#'      title("Saturated graph")	
#'	plot(dist.graph, col="grey")
#'	  points(ralu.site, col="red", pch=20, cex=1.5)
#'      box()
#'      title("Distance constrained graph")	  
#'		
#' @export			  
knn.graph <- function (x, row.names = NULL, k = NULL, max.dist = NULL, 
                       sym = FALSE, long.lat = FALSE, drop.diag = FALSE) 
   {
    if(is.null(k)) k=(dim(x)[1] - 1)
	  knn <- spdep::knearneigh(sp::coordinates(x), k = k, longlat = long.lat) 
      knn.nb <- spdep::knn2nb(knn, row.names = row.names, sym = sym)
    if(!is.na(sp::proj4string(x))) { 
	  prj <- sp::CRS(sp::proj4string(x))
	} else { 
	  prj <- sp::CRS(as.character(NA)) 
	}  
    if (!is.null(row.names)) {
      if (length(row.names) != knn$np) stop("row.names wrong length")
      if (length(unique(row.names)) != length(row.names)) stop("non-unique row.names given")
    }
    if (knn$np < 1) stop("non-positive number of spatial units")
    if (is.null(row.names)) row.names <- as.character(1:knn$np)
	graph <- spdep::nb2lines(knn.nb, coords = sp::coordinates(x), proj4string = prj)
	graph@data <- data.frame(graph@data[,1:4], length = sp::SpatialLinesLengths(graph, 
	                         longlat = long.lat))
	  names(graph@data)[3:4] <- c("from_ID","to_ID")
	    graph@data[,"from_ID"]<- as.character(graph@data[,"from_ID"])
        graph@data[,"to_ID"]<- as.character(graph@data[,"to_ID"])		
    rm.diag <- function(x) {
      ldiag <- function (x, diag=FALSE) {
            x <- as.matrix(x)
          if (diag) 
              row(x) >= col(x)
          else row(x) > col(x) 
        }  
        ctmx <- table(x@data$i, x@data$j) 	
          ctmx[ldiag(ctmx)] <- 0
          ctmx <- dmatrix.df(as.matrix(ctmx)) 
        ctmx <- data.frame(ij=paste(ctmx[,1], ctmx[,2], sep="."), dup=ctmx[,3])
        x@data <- cbind(ij=paste(x@data[,"i"], x@data[,"j"], sep="."), x@data)
          x <- merge(x, ctmx, by="ij")						 
          x <- x[x$dup == 1,]
        x@data <- x@data[,-which(names(x@data) %in% c("ij","dup"))]
      return(x)
    }
	if(drop.diag == TRUE) { graph <- rm.diag(graph)  }  
	if(!is.null(max.dist)) graph <- graph[graph$length <= max.dist,] 						 
  return( graph )
}

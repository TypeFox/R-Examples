#' @title Distance matrix to data.frame
#' @description Coerces distance matrix to a data.frame object
#' 
#' @param    x        Symmetrical distance matrix
#' @param    rm.diag  (TRUE/FALSE) remove matrix diagonal, self values. 
#' @return   data.frame object representing to and from values
#'
#' @note Function results in data.frame object with "X1" (FROM), "X2" (TO) and "distance" columns. The FROM column represents to origin ID, TO represents destination ID and distance is the associated matrix distance. These results can be joined back to the graph object using either the origin or destination ID's.  
#'
#' @author Jeffrey S. Evans  <jeffrey_evans@@tnc.org> and Melanie Murphy <melanie.murphy@@uwyo.edu>
#'
#' @examples 
#'   library(sp)
#'   pts <- cbind( x=runif(15, 480933, 504250), y=runif(15, 4479433, 4535122))
#'     pts <- SpatialPointsDataFrame(pts, data.frame(ID=paste("ob",1:nrow(pts),sep="")))
#'   
#'   # Create distance matrix  
#'   dm <- spDists(pts, pts)  
#'     colnames(dm) <- pts@@data[,"ID"] 
#'     rownames(dm) <- pts@@data[,"ID"]
#'   
#'   # Coerce to data.frame with TO and FROM ID's and associated distance
#'   dm.df <- dmatrix.df(dm)
#'     head(dm.df)
#'
#' @export
dmatrix.df <- function(x, rm.diag = TRUE) {
  if(rm.diag == TRUE) { diag(x) <- NA }
    varnames = names(dimnames(x))
    values <- as.vector(x)
    dn <- dimnames(x)
    char <- sapply(dn, is.character)
    dn[char] <- lapply(dn[char], utils::type.convert)
    indices <- do.call(expand.grid, dn)
      names(indices) <- varnames
      indices <- data.frame(indices, distance = values)
      indices <- stats::na.omit(indices)	
  return( indices )
}

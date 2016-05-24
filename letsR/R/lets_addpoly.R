#' Add polygon coverage to a PresenceAbscence object
#' 
#' @author Bruno Vilela
#' 
#' @description Add polygon coverage within cells of a PresenceAbsence object.
#'
#' @param x A \code{\link{PresenceAbsence}} object. 
#' @param y Polygon of interest.
#' @param z A character indicating the column name of the polygon containing the attributes to be used.
#' @param onlyvar If \code{TRUE} only the matrix object will be returned.
#' 
#' @return The result is a presence-absence matrix of species with 
#' the polygons' attributes used added as columns at the right-end of the matrix . The Values represent
#' the percentage of the cell covered by the polygon attribute used.   
#'  
#' @seealso \code{\link{lets.presab.birds}}
#' @seealso \code{\link{lets.presab}}
#' @seealso \code{\link{lets.addvar}}
#' 
#' @examples \dontrun{
#' data(PAM)  # Phyllomedusa presence-absence matrix
#' require(maptools)
#' data(wrld_simpl)  # World map
#' Brazil <- wrld_simpl[wrld_simpl$NAME == "Brazil", ]  # Brazil (polygon)
#' 
#' # Check where is the variable name 
#' # (in this case it is in "NAME" which will be my z value)
#' names(Brazil)
#' 
#' PAM_pol <- lets.addpoly(PAM, Brazil, "NAME")
#' }
#' 
#' @export

lets.addpoly <- function(x, y, z, onlyvar = FALSE){
  
  # Get the column position and change the name to a common one
  pos1 <- names(y) == z
  names(y)[pos1] <- "NOME"
  nomes <- y$NOME
  
  # Make the matrix
  valores <- values(x[[2]])
  LenValues <- length(valores)
  n <- nrow(y)
  matriz <- matrix(0, ncol = n, nrow = LenValues)
  colnames(matriz) <- nomes
  
  # Add coordinates
  xy <- xyFromCell(x[[2]], 1:LenValues)
  
  # Calculate the cell area
  areashape <- areaPolygon(y)
  areagrid <- NULL
  global <- all(as.vector(extent(x[[2]])) == c(-180, 180, -90, 90))
  
  if (!global) {
    grid <- rasterToPolygons(x[[2]])
    areagrid <- try(areaPolygon(grid), silent=TRUE)
  }
  
  if (class(areagrid) == "try-error" | global) {
    areagrid <- values(area(x[[2]])) * 1000000
  } 
  
  position <- which(!is.na(valores))
  
  for(i in 1:n){
    
    celu <- extract(x[[2]], y[i, ], cellnumber = TRUE, 
                    small = TRUE, weights = TRUE)
    celu2 <- do.call(rbind.data.frame, celu)
    
    if (!is.null(celu2[, 2])) {
      celu2 <- celu2[!is.na(celu2[, 2]), , drop = FALSE]
    }
    
    calc1 <- (celu2[, 3] * areashape[i])
    calc2 <- areagrid[position %in% celu2[, 1]]
    prop <- round(calc1 / calc2, 2)
    prop1 <- prop > 1
    
    if (any(prop1)) {
      prop[prop1] <- 1
    }
    matriz[celu2[, 1], i] <- prop
  }
  
  r <- rasterize(xy, x[[2]], matriz)
  rExtract <- extract(r, x[[1]][, 1:2])
  rExtract <- as.matrix(rExtract)
  colnames(rExtract) <- nomes
  
  if (onlyvar) {
    return(rExtract) 
  } else {
    resultado <- cbind(x[[1]], rExtract)
    return(resultado)
  }
}

# shp1.path="Desktop/shp/voro.shp"
# shp2.path="Desktop/shp/simulationFileUnit.shp"
# shp1.fieldID = "ID"
# shp2.fieldID = "ID"
# shp1.fieldOS = "OS"
# shp2.fieldOS = "TYPE"

#' Calculates vectorial kappa between two shapes
#' 
#' This function calculates the vectorial Kappa, the associated 
#' standard deviation and also returns the confusion matrix calculated 
#' between two vectorial landscapes.
#' @export KappaV
#' @param shp1.path \code{character}. The path to the first shape file.
#' @param shp2.path \code{character}. The path to the second shape file.
#' @param shp1.fieldID \code{character}. The column name in the .dbf file 
#' to indicate the polygons IDs.
#' @param shp2.fieldID \code{character}. The column name in the .dbf file 
#' to indicate the polygons IDs.
#' @param shp1.fieldOS \code{character}. The column name in the .dbf file 
#' to indicate the nominal variable of interest.
#' @param shp2.fieldOS \code{character}. The column name in the .dbf file
#'  to indicate the nominal variable of interest.
#' @param plot \code{logical}. Whether to plot the two landscapes.
#' @details If not specified the default parameters are
#'  \code{shp1.fieldID = "ID"}, \code{shp2.fieldID = shp1.fieldID}, 
#'  \code{shp1.fieldOS = "OS"},\code{shp2.fieldOS = shp1.fieldOS} 
#'  and \code{plot = FALSE}
#' @return A list with two components:
#' \code{confusion.matrix}: the confusion matrix with the corresponding areas.
#' and \code{$kappa.v}: a numeric with two values: \code{$Kappa} (the value of 
#' vectorial Kappa) and \code{$Kappa.sd }(the associated standard deviation)
#' @references Paper submitted.
#' @seealso \link{Kappa} in the \code{PresenceAbsence} package that 
#' handles the Kappa and its SD calculation from the confusion matrix.
#' 
#' Have a look to the package's vignette: 
#' \url{http://www.vincentbonhomme.fr/KappaV}
#' @examples
#' # Have a look to package's vignette above.
KappaV <-
  function(shp1.path, shp2.path,
           shp1.fieldID = "ID", shp2.fieldID = shp1.fieldID,
           shp1.fieldOS = "OS", shp2.fieldOS = shp1.fieldOS, plot = FALSE) {
    #we import the shp files
    shp1       <- readShapeSpatial(shp1.path)
    shp2       <- readShapeSpatial(shp2.path)
    
    # we create the list of poly(shp1) overlapping those of the shp2
    ov <- over(shp1, shp2, returnList=TRUE)
    
    # we retrieve the OS levels in the two shp files
    shp1.lev <- sort(unique(shp1@data[, shp1.fieldOS]))
    shp2.lev <- sort(unique(shp2@data[, shp2.fieldOS]))
    
    # if required, we plot the landscapes side by side and restore the graphical window
    if (plot) {
      layout(matrix(1:2, ncol=2))
      plot(shp1, col=terrain.colors(length(unique(shp1@data[, shp1.fieldOS])))[shp1@data[, shp1.fieldOS]+1])
      plot(shp2, col=terrain.colors(length(unique(shp2@data[, shp1.fieldOS])))[shp2@data[, shp2.fieldOS]+1])
      layout(matrix(1))}
    
    # full list of OS
    crn <- unique(c(shp1.lev, shp2.lev))
    
    # we prepare the confusion matrix
    res <- matrix(0, nrow=length(crn), ncol=length(crn),
                  dimnames=list(crn, crn))
    
    # we retrieve the nb of polygons in the two landscapes
    n1 <- length(shp1@polygons)
    n2 <- length(shp2@polygons)
    
    # we add another column to data 
    shp1.fieldID2 <- paste0(shp1.fieldID, "2")
    shp2.fieldID2 <- paste0(shp2.fieldID, "2")
    shp1@data[, shp1.fieldID2] <- 1:n1
    shp2@data[, shp2.fieldID2] <- 1:n2
    
    # we loop from every polygon on the first landscape
    for (i in 1:n1) {
      xi <- SpatialPolygons(list(shp1@polygons[[i]]))
      xi <- as(xi, "gpc.poly")
      # on those of the second landscapes
      for (j in seq(along=ov[[i]])) {
        yj     <- SpatialPolygons(list(shp2@polygons[[ov[[i]][j]]]))
        yj     <- as(yj, "gpc.poly")
        int.ij <- area.poly(intersect(xi, yj))
        ri <- which(shp1.lev == shp1@data[shp1@data[ , shp1.fieldID2]==i, shp1.fieldOS])
        ci <- which(shp2.lev == shp2@data[shp2@data[ , shp2.fieldID2]==ov[[i]][j], shp2.fieldOS])
        res[ri, ci] <- res[ri, ci] + int.ij
      }
    }
    #print(res)
    # we finally return both the matrix and the Kappa
    return(list(confusion.matrix=res, kappa.v=Kappa(res)))}

#' KappaV
#' 
#'
#' Calculates vectorial Kappa, i.e. congruence between two vectorial landscapes.
#' 
#' @seealso
#' KappaV' homepage : \url{http://www.vincentbonhomme.fr/KappaV} with tutorials
#' and hotline.
#' 
#' KappaV' GitHub repo : \url{https://github.com/vbonhomme/KappaV} to contribute,
#' among other things.
#' 
#' @references Paper submitted.
#' @importFrom rgeos area.poly
#' @importFrom maptools readShapeSpatial
#' @importFrom sp SpatialPolygons
#' @importFrom PresenceAbsence Kappa
#' @docType package
#' @name KappaV
#' @keywords Abtract
NULL

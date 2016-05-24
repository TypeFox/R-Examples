#' Create a presence-absence matrix based on species' point occurrences
#' 
#' @author Bruno Vilela & Fabricio Villalobos
#' 
#' @description Convert species' occurrences into a presence-absence matrix based on a user-defined grid.
#'
#' 
#' @param xy A matrix with geographic coordinates of species occurrences, first column is the longitude (or x), and the second latitude (or y).
#' @param species Character vector with species names, in the same order as the coordinates.
#' @param xmx Maximun longitude used to construct the grid in which the matrix will be based (i.e. the [gridded] geographic domain of interest)
#' @param xmn Minimun longitude used to construct the grid in which the matrix will be based (i.e. the [gridded] geographic domain of interest)
#' @param ymx Maximun latitude used to construct the grid in which the matrix will be based (i.e. the [gridded] geographic domain of interest)
#' @param ymn Minimun latitude used to construct the grid in which the matrix will be based (i.e. the [gridded] geographic domain of interest)
#' @param resol Numeric vector of length 1 or 2 to set the grid resolution.
#' @param remove.cells Logical, if \code{TRUE} the final matrix will not contain cells in the grid with a value of zero (i.e. sites with no species present).
#' @param remove.sp Logical, if \code{TRUE} the final matrix will not contain species that do not match any cell in the grid.
#' @param show.matrix Logical, if \code{TRUE} only the presence-absence matrix will be returned.
#' @param crs Character representign the PROJ.4 type description of a Coordinate 
#' Reference System (map projection) of the points.
#' @param count Logical, if \code{TRUE} a counting window will open.
#' 
#' 
#' @return The result is an object of class \code{\link{PresenceAbsence}} with the following objects:
#' @return \strong{Presence-Absence Matrix}: A matrix of species' presence(1) and absence(0) information. The first two columns contain the longitude (x) and latitude (y) of the cells' centroid (from the gridded domain used);
#' @return \strong{Richness Raster}: A raster containing species richness data;
#' @return \strong{Species name}: A character vector with species' names contained in the matrix.
#' @return *But see the optional argument \code{show.matrix}.
#' 
#'  
#' @details The function creates the presence-absence matrix based on a raster file. Depending on the cell size, extension used and number of species it may require a lot of memory, 
#'and may take some time to process it. Thus, during the process, if \code{count} argument is set \code{TRUE}, a counting window will open so you can see the progress (i.e. in what polygon the function is working). Note that the number of 
#'polygons is not the same as the number of species that you have (i.e. a species may have more than one polygon/shapefiles).
#' 
#' @seealso \code{\link{plot.PresenceAbsence}}
#' @seealso \code{\link{lets.presab.birds}}
#' @seealso \code{\link{lets.presab}}
#' @seealso \code{\link{lets.shFilter}}
#' 
#' @examples \dontrun{
#' species <- c(rep("sp1", 100), rep("sp2", 100),
#'              rep("sp3", 100), rep("sp4", 100))
#' x <- runif(400, min = -69, max = -51)
#' y <- runif(400, min = -23, max = -4)
#' xy <- cbind(x, y)
#' PAM <- lets.presab.points(xy, species, xmn = -93, xmx = -29, 
#'                           ymn = -57, ymx = 15)
#' summary(PAM)
#' # Species richness map
#' plot(PAM, xlab = "Longitude", ylab = "Latitude",
#'      main = "Species richness map (simulated)")
#' 
#' # Map of the specific species  
#' plot(PAM, name = "sp1")  
#' }
#' 
#' 
#' @export



lets.presab.points <- function(xy, species, xmn = -180, xmx = 180, ymn = -90, 
                               ymx = 90, resol = 1, remove.cells = TRUE, 
                               remove.sp = TRUE, show.matrix = FALSE, 
                               crs = CRS("+proj=longlat +datum=WGS84"), 
                               count = FALSE) {
  
  # Get species name
  if (is.factor(species)) {
    nomes <- levels(species)
    nomes <- nomes[nomes %in% species]
  } else {
    nomes <- levels(factor(species))
  }
  
  # Raster creation
  ras <- raster(xmn = xmn,
                xmx = xmx,
                ymn = ymn,
                ymx = ymx,
                crs = crs)
  res(ras) <- resol
  values(ras) <- 1
  
  # Coordinates xy
  l.values <- length(values(ras))
  coord <- xyFromCell(ras, 1:l.values)
  colnames(coord) <- c("Longitude(x)", "Latitude(y)")
  
  # Matrix creation
  n <- length(nomes)
  matriz <- matrix(0, ncol = n, nrow = l.values)
  colnames(matriz) <- nomes
  
  # With count window
  if (count) {
    
    # Do not set a new device in rstudio to avoid warnings()
    if (!"tools:rstudio" %in% search()) {
      dev.new(width = 2, height = 2, pointsize = 12)
      par(mar = c(0, 0, 0, 0))
    }
    
    # Loop start, running repetitions for the number of species (n) 
    for(i in 1:n){
      plot.new()
      text(0.5, 0.5, paste(paste("Total:", n, "\n","Species to go: ", (n - i))))
      celulas2 <- .extractpos.points(species, nomes[i], xy, ras)
      matriz[celulas2, i] <- 1
    }
    dev.off()
  }
  
  
  if (!count) {    
    for(i in 1:n) {
      celulas2 <- .extractpos.points(species, nomes[i], xy, ras)
      matriz[celulas2, i] <- 1
    }
  }  
  
  
  Resultado <- cbind(coord, matriz)
  
  if (remove.cells) {
    Resultado <- .removeCells(Resultado)
  }
  if (remove.sp) {
    Resultado <- .removeSp(Resultado)
  }
  
  if (show.matrix) {
    return(Resultado)
  } else {
    values(ras) <- rowSums(matriz)
    final <- list("Presence_and_Absence_Matrix" = Resultado,
                  "Richness_Raster"= ras, 
                  "Species_name"= colnames(Resultado)[-(1:2)])
    class(final) <- "PresenceAbsence"
    return(final)
  }
}


# Axuliar function to avoid code repetition inside the loop <<<<<<<<<

# Function to extract cell positions in the raster

.extractpos.points <- function(species, nomesi, xy, ras) {
  pos <- species == nomesi
  xy2 <- xy[pos, , drop = FALSE]
  celulas <- extract(ras, xy2, cellnumbers = T) [, 1]
  return(celulas)
}

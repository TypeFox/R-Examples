#' Create a presence-absence matrix of species' geographic ranges within a
#' grid for the Birdlife spatial data
#' 
#' @author Bruno Vilela & Fabricio Villalobos
#' 
#' @description Convert species' ranges (in shapefile format and stored in particular folders) 
#' into a presence-absence matrix based on a user-defined grid. This function is specially 
#' designed to work with BirdLife Intl. shapefiles (\url{http://www.birdlife.org}).
#'
#' 
#' @param path Path location of folders with one or more species' individual shapefiles. 
#' Shapefiles with more than one species will not work on this function. To use multi-species shapefiles see 
#' \code{\link{lets.presab}}.
#' @param xmx Maximun longitude used to construct the grid in which the matrix will be based 
#' (i.e. the [gridded] geographic domain of interest).
#' @param xmn Minimun longitude used to construct the grid in which the matrix will be based 
#' (i.e. the [gridded] geographic domain of interest).
#' @param ymx Maximun latitude used to construct the grid in which the matrix will be based 
#' (i.e. the [gridded] geographic domain of interest).
#' @param ymn Minimun latitude used to construct the grid in which the matrix will be based 
#' (i.e. the [gridded] geographic domain of interest)
#' @param resol Numeric vector of length 1 or 2 to set the grid resolution.
#' @param remove.cells Logical, if \code{TRUE} the final matrix will not contain cells in the 
#' grid with a value of zero (i.e. sites with no species present).
#' @param remove.sp Logical, if \code{TRUE} the final matrix will not contain species that do 
#' not match any cell in the grid.
#' @param show.matrix Logical, if \code{TRUE} only the presence-absence matrix will be returned.
#' @param crs Character representign the PROJ.4 type description of a Coordinate Reference 
#' System (map projection) of the original polygons.
#' @param crs.grid Character representign the PROJ.4 type description of
#' a Coordinate Reference System (map projection) for the grid. 
#' Note that when you change this options you may probably change 
#' the extent coordinates and the resolution.
#' @param cover Porcentage of the cell covered by the shapefile that will be considered for 
#' presence (values between 0 and 1). This option is only available when the coordinates 
#' are in degrees (longitude/latitude).
#' @param presence A vector with the code numbers for the presence type to be considered 
#' in the process (for IUCN spatial data \url{http://www.iucnredlist.org/technical-documents/spatial-data},
#' see metadata). 
#' @param origin A vector with the code numbers for the origin type to be considered in the process 
#' (for IUCN spatial data).
#' @param seasonal A vector with the code numbers for the seasonal type to be considered in the process 
#' (for IUCN spatial data).
#' @param count Logical, if \code{TRUE} a counting window will open.
#' 
#' @return The result is an object of class \code{\link{PresenceAbsence}} with the following objects:
#' @return \strong{Presence-Absence Matrix}: A matrix of species' presence(1) and absence(0) 
#' information. The first two columns contain the longitude (x) and latitude (y) of the cells'
#' centroid (from the gridded domain used);
#' @return \strong{Richness Raster}: A raster containing species richness data;
#' @return \strong{Species name}: A vector with species' names contained in the matrix.
#' @return *But see the optional argument \code{show.matrix}.
#'  
#' @details The function creates the presence-absence matrix based on a raster file. 
#' Depending on the cell size, extension used and number of species it may require a 
#' lot of memory, and may take some time to process it. Thus, during the process, if 
#' \code{count} argument is set \code{TRUE}, a counting window will open so you can 
#' see the progress (i.e. in what polygon the function is working).
#' Note that the number of polygons is not the same as the number of species that 
#' you have (i.e. a species may have more than one polygon/shapefiles).
#' 
#' @seealso \code{\link{plot.PresenceAbsence}}
#' @seealso \code{\link{lets.presab}}
#' @seealso \code{\link{lets.shFilter}}
#' 
#' @examples \dontrun{
#' # Constructing a Presence/Absence matrix for birds 
#' # Attention: For your own files, omit the 'system.file' 
#' # and 'package="letsR"', these are just to get the path
#' # to files installed with the package. 
#' path.Ramphastos <- system.file("extdata", package = "letsR")
#' PAM <- lets.presab.birds(path.Ramphastos, xmn = -93, xmx = -29, 
#'                          ymn = -57, ymx = 25)
#'                          
#' # Species richness map
#' plot(PAM, xlab = "Longitude", ylab = "Latitude",
#'      main = "Ramphastos species Richness")  
#' 
#' }
#' 
#' @export



lets.presab.birds <- function(path, xmn = -180, xmx = 180, ymn = -90, 
                              ymx = 90, resol = 1, remove.cells = TRUE,
                              remove.sp = TRUE, show.matrix = FALSE, 
                              crs = CRS("+proj=longlat +datum=WGS84"),
                              crs.grid = crs, cover = 0, presence = NULL,
                              origin = NULL, seasonal = NULL, count = FALSE) {
  
  # Shapefile list
  shapes <- list.files(path, pattern = ".shp$", 
                       full.names = T, recursive = T)
  
  # Raster creation
  r <- raster(xmn = xmn,
              xmx = xmx,
              ymn = ymn,
              ymx = ymx,
              crs = crs.grid)
  res(r) <- resol
  values(r) <- 0
  
  # Corrdinates saved
  valores <- values(r)
  xy <- xyFromCell(r, 1:length(valores))
  
  # Generating the empty matrix
  n <- length(shapes)
  nomes <- numeric()
  matriz <- matrix(0, nrow = nrow(xy), ncol = n)
  matriz <- cbind(xy, matriz)
  
  # Control for "after filtering none species distribution left" (see below) 
  k <- 0
  
  # Cell area calculation for cover metrics
  areashape <- NULL
  areagrid <- NULL
  
  if (cover > 0) {
    global <- xmn == -180 & xmx == 180 & ymn == -90 & ymx == 90    
    if (!global) {
      grid <- rasterToPolygons(r)
      areagrid <- try(areaPolygon(grid), silent=TRUE)
    } 
    
    if (class(areagrid) == "try-error" | global) {
      areagrid <- values(area(r)) * 1000000
    }  
  }
  
  
  
  # With count window  
  if (count) {
    
    # Do not set a new device in rstudio to avoid warnings()
    if (!"tools:rstudio" %in% search()) {
      dev.new(width = 2, height = 2, pointsize = 12)
      par(mar = c(0, 0, 0, 0))
    }
    
    # Loop start, running repetitions for the number of polygons (n) 
    for (j in 1:n) {    
      
      # Count window
      plot.new()
      text(0.5, 0.5, paste(paste("Total:", n, "\n",
                                 "Polygons to go: ",
                                 (n - j))))
      
      # Getting species position in the matrix and set to 1
      par.re <- .extractpos.birds(valores,  shapes[j], 
                                  k, r, areashape, areagrid,
                                  cover, presence, origin, 
                                  seasonal, crs, crs.grid)
      
      matriz[, (j + 2)] <- par.re[[1]]
      nomes[j] <- par.re[[2]]
      k <- par.re[[3]]
    }  
    dev.off()
  }
  
  # Without count window  
  if (!count) {
    
    for(j in 1:n){
      
      # Getting species position in the matrix and set to 1
      par.re <- .extractpos.birds(valores,  shapes[j], 
                                  k, r, areashape, areagrid,
                                  cover, presence, origin, 
                                  seasonal, crs, crs.grid)      
      matriz[, (j + 2)] <- par.re[[1]]
      nomes[j] <- par.re[[2]]
      k <- par.re[[3]]
    }  
  }
  
  if (k == 0) {
    stop("after filtering none species distribution left")
  }
  
  colnames(matriz) <- c("Longitude(x)", "Latitude(y)", nomes)
  
  riqueza <- rowSums(matriz[, -c(1, 2), drop = FALSE])  
  
  # Remove cells without presence if TRUE
  if (remove.cells) {
    matriz <- .removeCells(matriz)
  }
  
  # Remove species without presence if TRUE
  if (remove.sp) {
    matriz <- .removeSp(matriz)
  }
  
  # Putting together species with more than one shapefile
  matriz <- .unicas(matriz)
  
  # Return result (depending on what the user wants)
  if (show.matrix) {
    return(matriz)
  } else {
    values(r) <- riqueza
    final <- list("Presence_and_Absence_Matrix" = matriz, "Richness_Raster" = r, 
                  "Species_name" = (colnames(matriz)[-(1:2)]))
    class(final) <- "PresenceAbsence"
    return(final)
  }
}





# Axuliar function to avoid code repetition inside the loop <<<<<<<<<

# Function to extract cell positions in the raster

.extractpos.birds <- function(valores, shapej, k, r,
                              areashape, areagrid, cover,
                              presence, origin, seasonal, crs,
                              crs.grid) {
  
  # Vector to be filled
  valores2 <- valores
  
  # Read species shapefile, get its name and filter it
  shp <- readShapePoly(shapej, delete_null_obj = TRUE,
                       force_ring = TRUE, proj4string = crs)
  shp <- spTransform(shp, crs.grid)
  nomesj <- levels(shp$SCINAME)[1]
  shp <- lets.shFilter(shp, 
                       presence = presence,
                       origin = origin,
                       seasonal = seasonal)
  
  # Just run if after filtering the species has any polygon
  if (!is.null(shp)) {  
    k <- k + 1 # for later error control (see below)
    
    #  Extract cell occurrence positions
    cell <- extract(r, shp, cellnumber = T, 
                    small = T, weights = T)
    
    # Remove null cells
    cell <- cell[!sapply(cell, is.null)]
    
    # Changing colnames
    l.cell <- length(cell)
    
    if(l.cell > 0){
      .rename <- function(x) {
        colnames(x) <- 1:3
        return(x) 
      }          
      cell <- lapply(cell, .rename)
    }
    
    # Getting row positions
    cell2 <- do.call(rbind.data.frame, cell)    
    
    # Correcting presence based on the cover
    
    if (cover > 0) {
      areashape <- areaPolygon(shp)
      prop <- numeric()
      
      for (k1 in 1:l.cell) {
        calc <- cell[[k1]][, 3] * areashape[k1] / areagrid[cell[[k1]][, 1]]
        prop <- c(prop, calc)
      }
      prop.1 <- prop > 1
      
      if (any(prop.1)) {
        prop[prop.1] <- 1
      }
      
      cell2 <- cell2[prop >= cover, , drop = FALSE]
    }
    
    valores2[cell2[, 1]] <- 1
  }
  
  return(list(valores2, nomesj, k))
}

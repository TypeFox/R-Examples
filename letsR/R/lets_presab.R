#' Create a presence-absence matrix of species' geographic ranges within a grid
#' 
#' @author Bruno Vilela & Fabricio Villalobos
#' 
#' @description Convert species' ranges (in shapefile format) into a presence-absence matrix based on a user-defined grid system
#'
#' @param shapes Object of class \code{SpatialPolygonsDataFrame} (see function \code{\link{readShapePoly}} 
#' to open these files) containing the distribution of one or more species.
#' Species names should be in a column (within the .DBF table of the shapefile) 
#' called BINOMIAL/binomial or SCINAME/sciname.
#' @param xmx Maximun longitude used to construct the grid in which the matrix will be based 
#' (i.e. the [gridded] geographic domain of interest)
#' @param xmn Minimun longitude used to construct the grid in which the matrix will be based 
#' (i.e. the [gridded] geographic domain of interest)
#' @param ymx Maximun latitude used to construct the grid in which the matrix will be based 
#' (i.e. the [gridded] geographic domain of interest)
#' @param ymn Minimun latitude used to construct the grid in which the matrix will be based 
#' (i.e. the [gridded] geographic domain of interest)
#' @param resol Numeric vector of length 1 or 2 to set the grid resolution.
#' @param remove.cells Logical, if \code{TRUE} the final matrix will not contain cells in the 
#' grid with a value of zero (i.e. sites with no species present).
#' @param remove.sp Logical, if \code{TRUE} the final matrix will not contain species that 
#' do not match any cell in the grid.
#' @param show.matrix Logical, if \code{TRUE} only the presence-absence matrix will be returned.
#' @param crs Character representign the PROJ.4 type description of
#' a Coordinate Reference System (map projection) of the polygons.
#' @param crs.grid Character representign the PROJ.4 type description of
#' a Coordinate Reference System (map projection) for the grid. 
#' Note that when you change this options you may probably change 
#' the extent coordinates and the resolution.
#' @param cover Porcentage of the cell covered by the shapefile that will be considered for presence 
#' (values between 0 and 1). This option is only available when the coordinates 
#' are in degrees (longitude/latitude).
#' @param presence A vector with the code numbers for the presence type to be considered in the process 
#' (for IUCN spatial data \url{http://www.iucnredlist.org/technical-documents/spatial-data}, see metadata). 
#' @param origin A vector with the code numbers for the origin type to be considered in the process 
#' (for IUCN spatial data).
#' @param seasonal A vector with the code numbers for the seasonal type to be considered in the process 
#' (for IUCN spatial data).
#' @param count Logical, if \code{TRUE} a counting window will open.
#' 
#' 
#' @return The result is a list object of class \code{\link{PresenceAbsence}} with the following objects:
#' @return \strong{Presence-Absence Matrix}: A matrix of species' presence(1) and absence(0) information. 
#' The first two columns contain the longitude (x) and latitude (y) of the cells' centroid 
#' (from the gridded domain used);
#' @return \strong{Richness Raster}: A raster containing species richness data;
#' @return \strong{Species name}: A character vector with species' names contained in the matrix.
#' @return *But see the optional argument \code{show.matrix}.
#' 
#'  
#' @details The function creates the presence-absence matrix based on a raster object. 
#' Depending on the cell size, extension used and number of species it may require a lot of memory, 
#'and may take some time to process it. Thus, during the process, if \code{count} argument is 
#' set \code{TRUE}, a counting window will open so you can see the progress 
#' (i.e. in what polygon/shapefile the function is working). Note that the number of 
#'polygons is not the same as the number of species that you have 
#' (i.e. a species may have more than one polygon/shapefiles).
#' 
#' @seealso \code{\link{plot.PresenceAbsence}}
#' @seealso \code{\link{lets.presab.birds}}
#' @seealso \code{\link{lets.shFilter}} 
#' 
#' @examples \dontrun{
#' # Spatial distribution polygons of south american frogs 
#' # of genus Phyllomedusa. 
#' data(Phyllomedusa)
#' PAM <- lets.presab(Phyllomedusa, xmn = -93, xmx = -29,
#'                    ymn = -57, ymx = 15)
#' summary(PAM)
#' # Species richness map
#' plot(PAM, xlab = "Longitude", ylab = "Latitude",
#'      main = "Phyllomedusa species richness")
#' # Map of the specific species      
#' plot(PAM, name = "Phyllomedusa nordestina")
#' }
#' 
#' 
#' @import raster
#' @import maptools
#' @import maps
#' @import sp
#' @import rgdal 
#' 
#' @export


lets.presab <- function(shapes, xmn = -180, xmx = 180, ymn = -90, 
                        ymx = 90, resol = 1, remove.cells = TRUE, 
                        remove.sp = TRUE, show.matrix = FALSE, 
                        crs = CRS("+proj=longlat +datum=WGS84"),
                        crs.grid = crs, cover = 0, presence = NULL,
                        origin = NULL, seasonal = NULL, count = FALSE) {
  
  # Projection set for spatial polygons
  proj4string(shapes) <- crs
  if (as.character(crs) != as.character(crs.grid)) {
    shapes <- spTransform(shapes, crs.grid)
  }
  # Filter the species range distribution
  if (!all(is.null(presence), is.null(origin), is.null(seasonal))) {
    shapes <- lets.shFilter(shapes = shapes, 
                            presence = presence, 
                            origin = origin, 
                            seasonal = seasonal)
  }
  
  # Error control for no shapes after filtering
  if (is.null(shapes)) {
    stop("After filtering no species distributions left")
  }
  
  # Raster creation
  ras <- raster(xmn = xmn,
                xmx = xmx,
                ymn = ymn,
                ymx = ymx,
                crs = crs.grid)
  res(ras) <- resol
  values(ras) <- 1
  
  # Corrdinates saved
  ncellras <- ncell(ras)
  coord <- xyFromCell(object = ras, cell = 1:ncellras)
  colnames(coord) <- c("Longitude(x)", "Latitude(y)")
  
  # Cell area calculation for cover metrics
  areashape <- NULL
  areagrid <- NULL
  
  if (cover > 0) {
    areashape <- areaPolygon(shapes)
    global <- xmn == -180 & xmx == 180 & ymn == -90 & ymx == 90
    if (!global) {
      grid <- rasterToPolygons(ras)      
      areagrid <- try(areaPolygon(grid), silent = TRUE)
    }
    if (class(areagrid) == "try-error" | global) {
      areagrid <- values(area(ras)) * 1000000
    }
  }
  
  # Getting species name 
  names(shapes) <- toupper(names(shapes))
  names(shapes)[names(shapes) %in% "SCINAME"] <- "BINOMIAL" 
  nomes <- levels(shapes$BINOMIAL)
  n <- length(shapes$BINOMIAL)
  nomes2 <- shapes$BINOMIAL
  nomes <- nomes[nomes %in% nomes2]
  
  # Generating the empty matrix
  matriz <- matrix(0, ncol = length(nomes), nrow = ncellras)
  colnames(matriz) <- nomes
  
  # With count window
  if (count) {
    # Do not set a new device in rstudio to avoid warnings()
    if (!"tools:rstudio" %in% search()) {
      dev.new(width = 2, height = 2, pointsize = 12)
      par(mar = c(0, 0, 0, 0))
    }
    # Loop start, running repetitions for the number of polygons (n) 
    for (i in 1:n) {
      
      # Count window
      plot.new()
      text(0.5, 0.5, paste(paste("Total:", n, "\n",
                                 "Polygons to go: ",
                                 (n - i))))
      
      # Getting species position in the matrix and set to 1
      pospos2 <- .extractpos(ras, shapes@polygons[[i]], nomes, nomes2,
                             cover, areashape, areagrid, i)      
      matriz[pospos2$pos2[, 1], pospos2$pos] <- 1      
    }
    dev.off()
  }
  
  # Wihout count window
  if (!count) {
    
    # Loop start, running repetitions for the number of polygons (n) 
    for(i in 1:n){
      # Getting species position in the matrix and set to 1
      pospos2 <- .extractpos(ras, shapes@polygons[[i]], nomes, nomes2,
                             cover, areashape, areagrid, i)      
      matriz[pospos2$pos2[, 1], pospos2$pos] <- 1            
    }
  }  
  
  # Adding the coordinates to the matrix
  Resultado <- cbind(coord, matriz)
  
  # Remove cells without presence if TRUE
  if (remove.cells) {
    Resultado <- .removeCells(Resultado)
  }
  
  # Remove species without presence if TRUE
  if (remove.sp) {
    Resultado <- .removeSp(Resultado)
  }
  
  # Return result (depending on what the user wants)
  if (show.matrix) {
    return(Resultado)
  } else {
    values(ras) <- rowSums(matriz)
    final <- list("Presence_and_Absence_Matrix" = Resultado,
                  "Richness_Raster" = ras, 
                  "Species_name" = colnames(Resultado)[-(1:2)])
    class(final) <- "PresenceAbsence"
    return(final)
  }
}




# Axuliar function to avoid code repetition inside the loop <<<<<<<<<

# Function to extract cell positions in the raster
.extractpos <- function(ras, shapepol, nomes, nomes2, cover, 
                        areashape, areagrid, i) {
  
  # Try the extraction of cell occurrence positions
  celulas <- try(celulas <- extract(ras, 
                                    SpatialPolygons(list(shapepol)),
                                    cellnumbers = TRUE,
                                    weights = TRUE,
                                    small = TRUE), 
                 silent=T)
  
  # Handle the awkward error that can appear with weights and small = TRUE 
  if (class(celulas) == "try-error") {
    celulas <- extract(ras, SpatialPolygons(list(shapepol)),
                       cellnumbers = TRUE)
    nen <- sapply(celulas, nrow)
    for (ky in 1:length(nen)) {
      weigth <- rep(0, nen[ky]) 
      celulas[[ky]] <- cbind(celulas[[ky]], weigth)  
    } 
  }
  
  # Removing null cells 
  celulas <- celulas[!sapply(celulas, is.null)]
  
  # Changing colnames
  if (length(celulas) > 0) {
    .rename <- function(x) {
      colnames(x) <- 1:3
      return(x) 
    }
    celulas <- lapply(celulas, .rename)
  }
  
  # Getting column positions
  pos <- which(nomes2[i] == nomes)
  
  # Getting row positions
  pos2 <- do.call(rbind.data.frame, celulas)
  
  # Correcting presence based on the cover
  if (cover > 0) {  
    prop <- round((pos2[, 3] * areashape[i]) / areagrid[pos2[, 1]], 2)
    if (any(prop > 1)) {
      prop[prop > 1] <- 1
    }
    pos2 <- pos2[which(prop >= cover), , drop = FALSE]
  }
  
  # Return row and column position
  return (list("pos" = pos, "pos2" = pos2))
}

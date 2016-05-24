ssr <-
function(DIST=25, shp=data$p4no, colours=c("LightGreen","Tan"), PLOT=TRUE) {

  #==================================================
  #
  #  FILENAME:   ShrinkShape.r
  #  NAME:       ssr()
  #  AUTHOR:     Tarmo K. Remmel
  #  DATE:       19 January 2016
  #  NOTES:      Perform ShrinkShape fully in R
  #  NEEDS:      A shapefile (polygons), properly built
  #              libraries: sp, rgeos, rgdal
  #
  #  NOTES:      Shapefile can have only a single shape
  #              Shapes with holes allowed
  #              Projection units must be m, not dd or other angular unit
  #              Shapefile must be imported already (this tool does not do that)
  #
  #              To read Shapefile into object 'data'
  #              data <- readOGR(dsn="./Pond4IslandsPrj.shp", layer="Pond4IslandsPrj")
  #              Call ShrinkShape
  #              data$out <- ssr(DIST=25, shp=data)
  #    
  # CURRENTLY WITH FIXED INTERNAL BUFFER WIDTH AND SHAPEFILE NAME
  # FIX MAY BE WITH POLYGON HOLE COMMENT FUNCTIONS (rgeos)
  # SHAPEFILE MUST NOT BE 'DIRTY'. BE SURE TO SAVE PROPERLY CLEANED AND COMPLETE 
  # SHAPEFILES THAT HAVE CRS INFORMATION (PROJECTION DATA) ATTACHED
  # FULLY IMPLEMENTED IN R (NO SAGA REQUIRED) AND INCLUDES A PLOT OF THE SHRINKING
  # NOTE THAT PLOTTING GREATLY SLOWS THE ALGORITHM         
  #==================================================

  # LOAD REQUIRED LIBRARIES
  #library(rgeos)
  #library(sp)
  #library(rgdal)
  
  # ESTIMATE THE NUMBER OF ITERATIONS REQUIRED BASED ON SHAPEFILE BOUNDING BOX COORDINATES
  noiter <- ceiling(max(abs(shp@bbox[,1] - shp@bbox[,2])) / DIST) + 1
  
  # SET PLOTTING OPTIONS
  if(PLOT) {
    par(mfrow=c(2,3), pty="s")
    fill <- rep(colours,noiter)
    plot(shp, xlim=shp@bbox[1,], ylim=shp@bbox[2,], col=fill[1])
    title("Original Polygon")
    for (a in 1:noiter) {
      plot(gBuffer(shp, width=(-1*a*DIST)), xlim=shp@bbox[1,], ylim=shp@bbox[2,], xlab="X-Coordinate", ylab="Y-Coordinate", col=fill[a+1])
      par(new=TRUE)
    }
    title("Shrinking Phases")
  }
    
  # NOW BUILD VECTORS FOR ACCUMULATING AREA, PERIMETER, AND PARTS
  area <- rep(0,noiter)
  perim <- rep(0,noiter)
  parts <- rep(0,noiter)

  # NOW ACCUMULATE AREA, PERIMETER AND PARTS FOR FULL SHAPE
  area[1] <- gArea(shp)
  perim[1] <- gLength(shp)
  parts[1] <- length(shp@polygons[[1]]@Polygons)
  
  # PERFORM SHRINKING AND ACCUMULATE AREA, PERIMETER, AND PARTS FOR SHRUNKEN SHAPE
  for (a in 2:noiter) {
  	cat(a, ".", sep="")
    # PROCESS BUFFERING ONLY IF IT DOES NOT CAUSE THE SHAPE TO DISAPPEAR
    if(gIsEmpty(gBuffer(shp, width=(-1 * a * DIST))) == FALSE) {
      temp <- gBuffer(shp, width=(-1 * a * DIST))
      area[a] <- gArea(temp)
      perim[a] <- gLength(temp)
      parts[a] <- length(temp@polygons[[1]]@Polygons)
    }
  }

  # BIND CUMULATIVE SHRINKING DISTANCE, AREA, PERIMETER, AND PARTS INTO MATRIX
  # NEED TO ADJUST THE NEXT LINE TO PROPERLY ADJUST THE CUMULATIVE SHRINKING DISTANCE
  tab <- as.data.frame(cbind(area, perim, parts))
  # TRIM tab TO REMOVE EXCESS ZEROS
  cutrow <- min(which(tab[,2] == 0))
  tab <- tab[1:cutrow,]
  # ADD THE CUMULATIVE SHRINKING DISTANCE ATTRIBUTE COLUMN
  tab <- cbind(1:dim(tab)[1], tab)
  # ADD THE SHRINKING PHASE ATTRIBUTE COLUMN
  tab <- cbind(seq(0,noiter*DIST,DIST)[1:dim(tab)[1]], tab)
  # MAKE THE COLUMN NAMES NICE
  names(tab) <- c("CumShrink", "Phase", "Area", "Perimeter", "NumParts")

  if(PLOT) {
  # PLOT AREA, PERIMETER, AND NUMBER OF PARTS RESULTS
    par(new=FALSE)
    plot(tab$Area/10000, type="o", xlab="Shrinking Phase", ylab="Area (ha)")
    title("Area")
    par(new=FALSE)
    plot(tab$Perimeter, type="o", xlab="Shrinking Phase", ylab="Perimeter Length (m)")
    title("Perimeter")
    par(new=FALSE)
    plot(tab$NumParts, type="o", xlab="Shrinking Phase", ylab="Number of Parts")
    title("Parts")
  }

  # SEND RESULTS TABLE BACK AS FUNCTION RETURN
  return(tab)

}

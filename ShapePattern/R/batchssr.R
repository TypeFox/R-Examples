batchssr <-
function(outfile="out.csv", START=1) {
  
  #==================================================
  # 
  #  FILENAME:   ShrinkShape.r
  #  NAME:       batchssr()
  #  AUTHOR:     Tarmo K. Remmel
  #  DATE:       1 March 2016
  #  NOTES:      Calls ssr() -- Shrinkshape (for R) -- on each shapefile output by shpsplitter()
  #  NEEDS:      Shapefiles, each with one polygon, properly built
  #  NOTES:      Projection units must be m, not dd or other angular unit
  #    
  #==================================================

  #library(sp)
  #library(rgdal)
  #library(rgeos)

  files <- list.files(path="./", pattern=".shp")
  cat("Preparing to process ", length(files), " files...\n", sep="")

  # WRITE OUTPUT TABLE HEADERS TO .CSV FILE
  headers <- as.data.frame(t(c("Serial", "Shapefile", "Event", "Resolution", "UniqueID", "CumDist", "Iteration", "Area", "Perimeter", "NParts")))
  write.table(headers, file=outfile, append=TRUE, sep=",", col.names=FALSE, row.names=FALSE)

  for(i in START:length(files)) {

    # PROVIDE SCREEN FEEDBACK OF WHERE WE ARE IN THE LOOPING THROUGH OF FILES
    cat("  File: ", i, " of ", length(files), " --", files[i], "\n", sep="")
    
    # EXTRACT THE RESOLUTION FROM THE FILENAME
    startloc <- gregexpr(pattern="_", files[i])[[1]][1] + 2
    endloc <- gregexpr(pattern="_", files[i])[[1]][2] - 1
    res <- as.integer(substr(files[i], startloc, endloc))
    cat("  Shrinking distance: ", res, "\n", sep="")
    
    # EXTRACT THE SHAPEFILE NAME FROM THE FILENAME
    shpfilename <- strsplit(files[i], "[.]")[[1]][1]
    
    # EXTRACT THE UNIQUE SHAPE ID VALUE
    startloc <- gregexpr(pattern="_", files[i])[[1]][2] + 1
    endloc <- nchar(shpfilename)
    uqid <- as.integer(substr(files[i], startloc, endloc))
    
    # READ SHAPEFILE INTO R OBJECT 
    data <- readOGR(files[i], shpfilename)
    
    # CALL SHRINKSHAPE
    temptab <- ssr(DIST=res, shp=data, PLOT=FALSE)
    temptab <- as.data.frame(cbind(i, shpfilename, substr(files[i], 1, 3), res, uqid, temptab))
    write.table(temptab, file=outfile, append=TRUE, sep=",", col.names=FALSE, row.names=FALSE)
  } # END FOR: i


}

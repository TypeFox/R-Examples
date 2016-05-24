shpsplitter <-
function(SHAPEFILE="f09uqresids64", CODE="f09_R64_") {

  #==================================================
  #
  #  FILENAME:   ShrinkShape.r
  #  NAME:       shpsplitter()
  #  AUTHOR:     Tarmo K. Remmel
  #  DATE:       29 February 2016
  #  NOTES:      Split shapefile into multiple shapefiles, one shape in each file
  # 			 This is necessary for ssr2() to run
  #  NEEDS:      A shapefile (polygons), properly built
  #              libraries: sp, rgdal
  #
  #              shpsplitter(SHAPEFILE="f09uqresids64", CODE="f09_R64_")
  #
  #  NOTES:      Projection units must be m, not dd or other angular unit
  #              Unique shapes must be identified with integer code in attribute GRIDCODE
  #              CODE allows for a unique naming convention to be applied to related parts
  #              of a multi-part shapefile
  #    
  #==================================================

  # LOAD REQUIRED LIBRARIES
  #library(sp)
  #library(rgdal)

  
  filename <- paste("./", SHAPEFILE, ".shp", sep="")
  data <- readOGR(filename, SHAPEFILE)
  unique <- unique(data@data$GRIDCODE)

  for (i in 1:length(unique)) {
    tmp <- data[data$GRIDCODE == unique[i], ]
    writeOGR(tmp, dsn=getwd(), paste(CODE, unique[i], sep=""), driver="ESRI Shapefile", check_exists=TRUE, overwrite_layer=TRUE)
  } # END FOR: i

}

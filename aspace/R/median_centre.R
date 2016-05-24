"median_centre" <-
function(id=1, filename="median_centre_Output.txt", points=activities) {

  #=======================================================
  #
  #  TITLE:     MEDIAN CENTRE CALCULATOR
  #  FUNCTION:  median_centre()
  #  AUTHOR:    RANDY BUI, RON BULIUNG, TARMO K. REMMEL
  #  DATE:      March 28, 2011
  #  NOTES:     USE THE id PARAMETER TO SPECIFY A UNIQUE IDENTIFIER FOR
  #             THE MEDIAN CENTRE; THIS VALUE IS ADDED TO THE OUTPUT filename
  #             AS AN IDENTIFIER THAT CAN BE USED TO EXTRACT RECORDS WHEN 
  #             A USER EMBEDDS THE FUNCTION IN A LOOP TO GENERATE
  #             MULTIPLE MEDIAN CENTRES TO THE SAME FILE.
  #             THE filename PARAMETER CONTROLS WHERE THE COORDINATE INFORMATION 
  #             IS WRITTEN TO. USE medianloc (coordinates) and medianatt (attributes) 
  #             TO PRODUCE SHAPEFILES USING THE CONVERT.TO.SHAPEFILE AND WRITE.SHAPEFILE 
  #             FUNCTIONS FROM THE SHAPEFILES LIBRARY.
  #
  #  OUTPUT:	
  #     ID  		UNIQUE MEDIAN CENTRE IDENTIFIER
  #		median.x	X-COORDINATE OF THE MEDIAN CENTRE
  #		median.y	Y-COORDINATE OF THE MEDIAN CENTRE
  #		medianatt	ATTRIBUTES ABOVE WRITTEN TO DATAFRAME FOR POST-PROCESSING AS SHAPEFILE
  #		medianloc	UNIQUE ID AND X,Y COORDINATES OF THE POINT FOR POST-PROCESSING AS SHAPEFILE
  #
  #=======================================================
  
    # COMPUTE THE MEDIAN CENTRE
	median.x <- median(points[,1])
	median.y <- median(points[,2])
  
	# STORE COORDINATES OF THE MEDIAN CENTRE      
	coordsMC <- cbind(median.x, median.y)
    
    # CREATE ASCII OUTPUT FOR SHAPEFILE CREATION
    medianloc <- cbind(id, coordsMC)
	colnames(medianloc)=c("id","x","y")
    write.table(medianloc, sep=",", file=filename, col.names=FALSE)

	# DATA FRAME WITH COLUMNS IN ORDER ID, X-COORD, Y-COORD FOR CONVERT.TO.SHAPEFILE FUNCTION
	assign("medianloc", medianloc, pos=1)	
	
	# STORE RESULTS INTO A LIST (REQUIRED FOR PLOT FUNCTION)
	r.median <- list(id = id, points = points, median.x = median.x, median.y = median.y) 
	assign("r.median", r.median, pos=1)
    
    # STORE MEDIAN CENTRE ATTRIBUTES INTO A DATA FRAME AND PRINTS RESULTS
    result.median <- list("id"=id, "median.x"=median.x, "median.y"=median.y)
	result.median<-as.data.frame(result.median)
	print(result.median)  
	
	# DATA FRAME OF ATTRIBUTES WITH FIRST COLUMN NAME "ID" FOR CONVERT.TO.SHAPEFILE FUNCTION
	assign("medianatt", result.median, pos=1) 
}
"mean_centre" <-
function(id=1, filename="mean_centre_Output.txt", weighted=FALSE, weights=NULL, points=activities) {

  #=======================================================
  #
  #  TITLE:     MEAN CENTRE CALCULATOR
  #  FUNCTION:  mean_centre()
  #  AUTHOR:    RANDY BUI, RON BULIUNG, TARMO K. REMMEL
  #  DATE:      March 28, 2011
  #  NOTES:     USE THE id PARAMETER TO SPECIFY A UNIQUE IDENTIFIER FOR
  #             THE MEAN CENTRE; THIS VALUE IS ADDED TO THE OUTPUT filename
  #             AS AN IDENTIFIER THAT CAN BE USED TO EXTRACT RECORDS WHEN 
  #             A USER EMBEDDS THE FUNCTION IN A LOOP TO GENERATE
  #             MULTIPLE MEAN CENTRES TO THE SAME FILE.
  #             THE filename PARAMETER CONTROLS WHERE THE COORDINATE INFORMATION 
  #             IS WRITTEN TO. USE meanloc (coordinates) and meanatt (attributes) 
  #             TO PRODUCE SHAPEFILES USING THE CONVERT.TO.SHAPEFILE AND WRITE.SHAPEFILE 
  #             FUNCTIONS FROM THE SHAPEFILES LIBRARY.
  #
  #  OUTPUT:	
  #     ID  		UNIQUE MEAN CENTRE IDENTIFIER
  #		weighted    T|F SHOULD THE CENTRE BE WEIGHTED (WEIGHTED MEAN CENTER)
  #		weights		WEIGHTS (IF APPLIED)
  #		CENTRE.x	X-COORDINATE OF THE CENTRE
  #		CENTRE.y	Y-COORDINATE OF THE CENTRE
  #		meanatt		ATTRIBUTES ABOVE WRITTEN TO DATAFRAME FOR POST-PROCESSING AS SHAPEFILE
  #		meanloc		UNIQUE ID AND X,Y COORDINATES OF THE POINT FOR POST-PROCESSING AS SHAPEFILE
  #
  #=======================================================

  # STORE THE COUNT OF POINTS
  n <- dim(points)[1]
  
  # INITIALIZE centre.xy
  centre.xy <- NULL
  
  if(weighted) {
	# WEIGHT THE POINTS
	wt.x <- points[,1] * weights 
	wt.y <- points[,2] * weights
	
	# COMPUTE AND USE WEIGHTED MEAN CENTRE RATHER THAN USE SPECIFIED LOCATION AS CENTRE (WEIGHTED)
	WMC.x <- c( sum(wt.x) / sum(weights) )
	WMC.y <- c( sum(wt.y) / sum(weights) )    
	centre.xy[1] <- WMC.x
	centre.xy[2] <- WMC.y
   }
  else {
	# COMPUTE AND USE MEAN CENTRE RATHER THAN USE SPECIFIED LOCATION AS CENTRE (NON-WEIGHTED)
	meanx <- sum(points[,1])/n
	meany <- sum(points[,2])/n
	centre.xy[1] <- meanx
	centre.xy[2] <- meany
  } 

	# STORE COORDINATES OF THE MEAN CENTRE      
	coordsMC <- cbind(centre.xy[1], centre.xy[2])
    
    # CREATE ASCII OUTPUT FOR SHAPEFILE CREATION
    meanloc <- cbind(id, coordsMC)
	colnames(meanloc)=c("id","x","y")
    write.table(meanloc, sep=",", file=filename, col.names=FALSE)

	# DATA FRAME WITH COLUMNS IN ORDER ID, X-COORD, Y-COORD FOR CONVERT.TO.SHAPEFILE FUNCTION
	assign("meanloc", meanloc, pos=1)	
	
	# STORE RESULTS INTO A LIST (REQUIRED FOR PLOT FUNCTION)
	r.mean <- list(id = id, points = points, weighted = weighted, weights = weights, 
	              CENTRE.x = centre.xy[1], CENTRE.y = centre.xy[2]) 
	assign("r.mean", r.mean, pos=1)
    
    # STORE MEAN CENTRE ATTRIBUTES INTO A DATA FRAME AND PRINTS RESULTS
    result.mean <- list("id"=id, "CENTRE.x"=centre.xy[1], "CENTRE.y"=centre.xy[2])
	result.mean<-as.data.frame(result.mean)
	print(result.mean)  
	
	# DATA FRAME OF ATTRIBUTES WITH FIRST COLUMN NAME "ID" FOR CONVERT.TO.SHAPEFILE FUNCTION
	assign("meanatt", result.mean, pos=1)
}


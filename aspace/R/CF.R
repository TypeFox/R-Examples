"CF" <-
function(id=1, filename="CF_Output.txt", points=activities) {

  #=======================================================
  #
  #  TITLE:     CENTRAL FEATURE (CF) CALCULATOR 
  #  FUNCTION:  CF()
  #  AUTHOR:    RANDY BUI, RON BULIUNG
  #  DATE:      March 28, 2011
  #  CALLS:     distances()
  #  NOTES:     USE THE id PARAMETER TO SPECIFY A UNIQUE IDENTIFIER FOR
  #             THE CF; THIS VALUE IS ADDED TO THE OUTPUT filename
  #             AS AN IDENTIFIER THAT CAN BE USED TO EXTRACT RECORDS WHEN 
  #             A USER EMBEDDS THE FUNCTION IN A LOOP TO GENERATE
  #             MULTIPLE CF POINTS TO THE SAME FILE.
  #             THE filename PARAMETER CONTROLS WHERE THE COORDINATE INFORMATION 
  #             IS WRITTEN TO. USE cfloc (coordinates) and cfatt (attributes) 
  #             TO PRODUCE SHAPEFILES USING THE CONVERT.TO.SHAPEFILE AND WRITE.SHAPEFILE 
  #             FUNCTIONS FROM THE SHAPEFILES LIBRARY.
  #
  #  OUTPUT:	
  #     ID  	UNIQUE CF IDENTIFIER
  #		CF.x	X-COORDINATE OF THE CENTRAL FEATURE
  #		CF.y	Y-COORDINATE OF THE CENTRAL FEATURE
  #		cfatt	ATTRIBUTES ABOVE WRITTEN TO DATAFRAME FOR POST-PROCESSING AS SHAPEFILE
  #		cfloc	UNIQUE ID AND X,Y COORDINATES OF THE POINT FOR POST-PROCESSING AS SHAPEFILE
  #
  #=======================================================
  
  # DETERMINE THE CENTRAL FEATURE
  count.CF <- length(points[,1])	
  M.CF <- matrix(0,nrow=count.CF,ncol=3)				

  	for(i in 1:count.CF) {
		row.CF <- points[i,]
		coord.CF <- c(row.CF[,1],row.CF[,2])

		dist.CF <- distances(centre.xy=coord.CF, points, verbose=FALSE)
		sum.dist.CF <- sum(dist.CF)
					
		M.CF[i,1] <- sum.dist.CF
		M.CF[i,2] <- coord.CF[1]
		M.CF[i,3] <- coord.CF[2]
	}

  order.CF <- M.CF[order(M.CF[,1]),]
  first.row.CF <- order.CF[1,]
	
  x.CF <- first.row.CF[2]
  y.CF <- first.row.CF[3]

	# STORE COORDINATES OF THE CF     
	coordsCF <- cbind(x.CF, y.CF)
    
    # CREATE ASCII OUTPUT FOR SHAPEFILE CREATION
    cfloc <- as.data.frame(cbind(id, coordsCF))
	colnames(cfloc)=c("id","x","y")
    write.table(cfloc, sep=",", file=filename, col.names=FALSE)
	
	# DATA FRAME WITH COLUMNS IN ORDER ID, X-COORD, Y-COORD FOR CONVERT.TO.SHAPEFILE FUNCTION
	assign("cfloc", cfloc, pos=1)	
	
	# STORE RESULTS INTO A LIST (REQUIRED FOR PLOT FUNCTION)
	r.CF <- list(id = id, points = points, CF.x = x.CF, CF.y = y.CF) 
	assign("r.CF", r.CF, pos=1)
    
    # STORE CF ATTRIBUTES INTO A DATA FRAME AND PRINTS RESULTS
    result.CF <- list("id"=id, "CF.x"=x.CF, "CF.y"=y.CF)
	result.CF<-as.data.frame(result.CF)
	print(result.CF)	
	
	# DATA FRAME OF ATTRIBUTES WITH FIRST COLUMN NAME "ID" FOR CONVERT.TO.SHAPEFILE FUNCTION
	assign("cfatt", result.CF, pos=1)	
}
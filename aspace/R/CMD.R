"CMD" <-
function(id=1, filename="CMD_Output.txt", dist=100, points=activities) {

  #=======================================================
  #
  #  TITLE:     CENTRE OF MINIMUM DISTANCE (CMD) CALCULATOR 
  #  FUNCTION:  CMD() 
  #  AUTHOR:    RON BULIUNG, RANDY BUI
  #  DATE:      March 28, 2011
  #  CALLS:     distances(), gridpts()
  #  NEEDS:     LIBRARIES: splancs
  #  NOTES:     USE THE id PARAMETER TO SPECIFY A UNIQUE IDENTIFIER FOR
  #             THE CMD; THIS VALUE IS ADDED TO THE OUTPUT filename
  #             AS AN IDENTIFIER THAT CAN BE USED TO EXTRACT RECORDS WHEN 
  #             A USER EMBEDDS THE FUNCTION IN A LOOP TO GENERATE
  #             MULTIPLE CMD POINTS TO THE SAME FILE.
  #             THE filename PARAMETER CONTROLS WHERE THE COORDINATE INFORMATION 
  #             IS WRITTEN TO. USE cmdloc (coordinates) and cmdatt (attributes) 
  #             TO PRODUCE SHAPEFILES USING THE CONVERT.TO.SHAPEFILE AND WRITE.SHAPEFILE 
  #             FUNCTIONS FROM THE SHAPEFILES LIBRARY.
  #
  #  OUTPUT:	
  #     ID  				UNIQUE CMD IDENTIFIER
  #		CMD.x				X-COORDINATE OF THE CENTRE OF MINIMUM DISTANCE
  #		CMD.y				Y-COORDINATE OF THE CENTRE OF MINIMUM DISTANCE
  #		Distance			HOLD DISTANCE VALUE BETWEEN I AND ITH ITERATIONS (metres)
  #		Number of Cells		HOLD NUMBER OF CELLS IN EACH GRID CREATED FOR EACH ITERATION
  #		cmdatt				ATTRIBUTES ABOVE WRITTEN TO DATAFRAME FOR POST-PROCESSING AS SHAPEFILE
  #		cmdloc				UNIQUE ID AND X,Y COORDINATES OF THE POINT FOR POST-PROCESSING AS SHAPEFILE
  #
  #=======================================================

	# SET DEPENDENCIES
	require(splancs)

	# CREATE EMPTY OBJECTS
	x<-c() #HOLD X-COORD OF CMD FOR EACH ITERATION
	y<-c() #HOLD Y-COORD OF CMD FOR EACH ITERATION
	d<-c() #HOLD DISTANCE VALUE BETWEEN I AND ITH ITERATIONS
	n<-c() #HOLD ITERATION NUMBER
	cells<-c() #HOLD NUMBER OF CELLS IN EACH GRID CREATED FOR EACH ITERATION

	# INITIALIZE OBJECTS, COUNTERS, GRID SIZE
	i<-1
	x[i]<-0 
	y[i]<-0
	d[i]<-dist 
	n[i]<-0 
	cells[i]<-0 
	dx<-1 #INITIALIZE GRID SPACING IN X
	dy<-1 #INITIALIZE GRID SPACING IN Y, LARGER NUMBER = MORE CELLS, HERE IT IS SET TO 1 CELLS IN X,Y

	# GENERATE MCP
	hpts<-chull(points);MCP<-cbind(points[hpts,1],points[hpts,2])

	while (d[i] >= dist) {
	
	grid<-gridpts(MCP,dx,dy)  
	M.CMD<-matrix(0,nrow=nrow(grid),ncol=3)

		for(j in 1:nrow(grid)) {
				
				coord.CMD<-grid[j,]

				sumdist.CMD<-sum(distances(centre.xy=coord.CMD, points))
				
				M.CMD[j,1]<-sumdist.CMD
				M.CMD[j,2]<-coord.CMD[1]
				M.CMD[j,3]<-coord.CMD[2]		
		}

		if (i >= 1) { 
				order.CMD<-M.CMD[order(M.CMD[,1]),] 
				CMD<-order.CMD[1,]
		
		}else (CMD<-M.CMD[1,])

				
		#DUMP CMD FOR EACH ITERATION
		
		x[i+1]<-CMD[2]
		y[i+1]<-CMD[3]

		d[i+1]<-sqrt((x[i+1] - x[i])^2 + (y[i+1] - y[i])^2) #ESTIMATE DISTANCE BETWEEN CURRENT AND PREVIOUS CMD
		
		n[i+1]<-dx
		
		cells[i+1]<-nrow(grid) 

		rm(grid)
		i<-i+1
		dx<-dx+1
		dy<-dy+1
	}

	#PROCESS RESULTS, CMD HOLDS THE MINIMUM DISTANCE COORDINATES
	result<-cbind(n,round(x,2),round(y,2),round(d,2),cells)
	result<-as.data.frame(result[3:nrow(result),])
	result[,1]<-seq(1,nrow(result),1)
	colnames(result)<-c("n","X","Y","Dist","Cells")

	#CENTRE OF MINIMUM DISTANCE returned from simulation
	CMD<-cbind(result[nrow(result),2:5])
	colnames(CMD)<-c("X","Y","Dist","Cells")

	# STORE COORDINATES OF THE CMD     
	coordsCMD <- cbind(CMD[,1], CMD[,2])

    # CREATE ASCII OUTPUT FOR SHAPEFILE CREATION
    cmdloc <- as.data.frame(cbind(id, coordsCMD))
	colnames(cmdloc)=c("id","x","y")
    write.table(cmdloc, sep=",", file=filename, col.names=FALSE)

	# DATA FRAME WITH COLUMNS IN ORDER ID, X-COORD, Y-COORD FOR CONVERT.TO.SHAPEFILE FUNCTION
	assign("cmdloc", cmdloc, pos=1)	
	
	# STORE RESULTS INTO A LIST (REQUIRED FOR PLOT FUNCTION)
	r.CMD <- list(id = id, points = points, CMD.x = CMD[,1], CMD.y = CMD[,2], result = result) 
	assign("r.CMD", r.CMD, pos=1)
    
    # STORE CMD ATTRIBUTES INTO A DATA FRAME AND PRINTS RESULTS
    result.CMD <- list("id"=id, "CMD.x"=CMD[,1], "CMD.y"=CMD[,2], "Distance"=CMD[,3], "Number of Cells"=CMD[,4])
	result.CMD<-as.data.frame(result.CMD)
	print(result.CMD)

	# DATA FRAME OF ATTRIBUTES WITH FIRST COLUMN NAME "ID" FOR CONVERT.TO.SHAPEFILE FUNCTION
	assign("cmdatt", result.CMD, pos=1)	
}
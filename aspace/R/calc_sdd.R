"calc_sdd" <-
function(id=1, filename="SDD_Output.txt", centre.xy=NULL, calccentre=TRUE, weighted=FALSE, weights=NULL, points=activities, verbose=FALSE) {

  #=======================================================
  #
  #  TITLE:     STANDARD DEVIATION DISTANCE (SDD) CALCULATOR
  #  FUNCTION:  calc_sdd()
  #  AUTHOR:    RANDY BUI, RON BULIUNG, TARMO K. REMMEL
  #  DATE:      March 28, 2011
  #  CALLS:     distances(), as_radians()
  #  NOTES:     USE THE id PARAMETER TO SPECIFY A UNIQUE IDENTIFIER FOR
  #             THE SDD CIRCLE; THIS VALUE IS ADDED TO THE OUTPUT filename
  #             AS AN IDENTIFIER THAT CAN BE USED TO EXTRACT RECORDS WHEN 
  #             A USER EMBEDDS THE FUNCTION IN A LOOP TO GENERATE
  #             MULTIPLE SDD CIRCLES TO THE SAME FILE.
  #             THE filename PARAMETER CONTROLS WHERE THE COORDINATE INFORMATION 
  #             IS WRITTEN TO. USE sddloc (coordinates) and sddatt (attributes) 
  #             TO PRODUCE SHAPEFILES USING THE CONVERT.TO.SHAPEFILE AND WRITE.SHAPEFILE 
  #             FUNCTIONS FROM THE SHAPEFILES LIBRARY.
  #
  #  ERROR:     1000  NO ERRORS DETECTED
  #               25  TOO FEW ACTIVITIES, NEED >= 3
  #               21  INVALID COMBINATION: calccentre=TRUE and centre.xy!=NULL
  #
  #  OUTPUT:	
  #     ID  		UNIQUE SDD IDENTIFIER
  #		calccentre	T|F SHOULD THE MEAN CENTRE BE USED
  #		weighted    T|F SHOULD THE CENTRE BE WEIGHTED (WEIGHTED MEAN CENTER)
  #		CENTRE.x	X-COORDINATE OF THE CENTRE
  #		CENTRE.y	Y-COORDINATE OF THE CENTRE
  #		SDD.radius	RADIUS OF SDD
  #		SDD.area	AREA OF SDD
  #		sddatt		ATTRIBUTES ABOVE WRITTEN TO DATAFRAME FOR POST-PROCESSING AS SHAPEFILE
  #		sddloc		UNIQUE ID AND X,Y COORDINATES OF VERTICES FOR POST-PROCESSING INTO SDD SHAPEFILE
  #=======================================================
  
  # INITIALIZE ERROR CODE TO NO ERROR
  errorcode <- 1000

  # STORE THE COUNT OF POINTS
  n <- dim(points)[1]
  
  if(calccentre) {
    if(length(centre.xy) == 2) {
	  # ERROR: INVALID COMBINATION: calccentre=TRUE AND centre.xy!=NULL
      # SET DESCRIPTIVE ERROR CODE AND GIVE WARNING
      errorcode <- 21
      cat("\n\nWARNING: Invalid combination: calccentre=TRUE and centre.xy!=NULL")
	  cat("\nERROR CODE: ", errorcode, "\n\n", sep="")
      return("ERROR")
	}
	else {
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
    }
  }
  
  # INITIALIZE FUNCTION VARIABLE WITH PARAMETER VALUE
  dist <- distances(centre.xy, points)
  
  # TEST WHETHER A SUFFICIENT NUMBER OF POINTS WERE SUPPLIED
  if(length(dist) >= 3) {
	
	if(weighted) {		
	#PERFORM THE WEIGHTED STANDARD DEVIATION DISTANCE COMPUTATION (WEIGHTED SDD)
	SDD <- sqrt(sum((weights*dist^2)/((sum(weights)) - 2) ) )
	}
	else {
	# PERFORM THE STANDARD DEVIATION DISTANCE COMPUTATION (UNWEIGHTED SDD)
	SDD <- sqrt(sum(dist^2/(length(dist) - 2) ) )
	}
	
    # COMPUTE SDD AREA
    sddarea <- pi * SDD^2
  
    # COMPUTE AND STORE COORDINATES FOR PLOTTING THE SDD CIRCLE (BASED ON ELLIPSEPOINTS FUNCTION FROM SFSMISC LIBRARY)      
    B <- min(SDD, SDD)
    A <- max(SDD, SDD)
    d2 <- (A - B) * (A + B)
    phi <- 2 * pi * seq(0, 1, len = 360)
    sp <- sin(phi)
    cp <- cos(phi)
    r <- SDD * SDD/sqrt(B^2 + d2 * sp^2)
	xy <- r * cbind(cp, sp)
    al <- 0 * pi/180
    ca <- cos(al)
    sa <- sin(al)
    coordsSDD <- xy %*% rbind(c(ca, sa), c(-sa, ca)) + cbind(rep(centre.xy[1], 360), rep(centre.xy[2], 360))
    
    # CREATE ASCII OUTPUT FOR SHAPEFILE CREATION
	sddloc <- as.data.frame(cbind(id, coordsSDD))
	colnames(sddloc)=c("id","x","y")
    write.table(sddloc, sep=",", file=filename, col.names=FALSE)
	
	# DATA FRAME WITH COLUMNS IN ORDER ID, X-COORD, Y-COORD FOR CONVERT.TO.SHAPEFILE FUNCTION
	assign("sddloc", sddloc, pos=1)	

	# STORE RESULTS INTO A LIST (REQUIRED FOR PLOT FUNCTION)
	r.SDD <- list(id = id, points = points, coordsSDD = coordsSDD, SDD = SDD, calccentre = calccentre, weighted = weighted, weights = weights, 
	                   CENTRE.x = centre.xy[1], CENTRE.y = centre.xy[2], SDD.area = sddarea) 
	assign("r.SDD", r.SDD, pos=1)
    
    # STORE SDD ATTRIBUTES INTO A DATA FRAME AND PRINTS RESULTS
    result.sdd <- list("id"=id, "calccentre"=calccentre, "weighted"=weighted, "CENTRE.x"=centre.xy[1], "CENTRE.y"=centre.xy[2],
					   "SDD.radius"=SDD, "SDD.area"=sddarea)
	print(result.sdd)
	result.sdd<-as.data.frame(result.sdd)

	# DATA FRAME OF ATTRIBUTES WITH FIRST COLUMN NAME "ID" FOR CONVERT.TO.SHAPEFILE FUNCTION
	assign("sddatt", result.sdd, pos=1)	
  }
  else {
    # ERROR: TOO FEW POINTS: NEED >= 3
    # SET DESCRIPTIVE ERROR CODE AND GIVE WARNING
    errorcode <- 25
    if(verbose) {
      cat("\n\nWARNING: Not enough values to compute SDD.")
      cat("\nERROR CODE: ", errorcode, "\n\n", sep="")
    }
    return("ERROR")
  }
}


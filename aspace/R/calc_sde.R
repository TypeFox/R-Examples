"calc_sde" <-
function(id=1, filename="SDE_Output.txt", centre.xy=NULL, calccentre=TRUE, weighted=FALSE, weights=NULL, points=activities, verbose=FALSE) {

  #=======================================================
  #
  #  TITLE:     STANDARD DEVIATION ELLIPSE (SDE) CALCULATOR
  #  FUNCTION:  calc_sde()
  #  AUTHOR:    RANDY BUI, RON BULIUNG, TARMO REMMEL
  #  DATE:      August 08, 2012
  #  CALLS:     atan_d(), sin_d(), cos_d(), atan_d()
  #  NOTES:     NOTE THAT R TRIGONOMETRIC FUNCTIONS ARE IN RADIANS NOT DEGREES.
  #             WMC:  WEIGHTED MEAN CENTRE
  #				USE THE id PARAMETER TO SPECIFY A UNIQUE IDENTIFIER FOR
  #             THE SDE; THIS VALUE IS ADDED TO THE OUTPUT filename
  #             AS AN IDENTIFIER THAT CAN BE USED TO EXTRACT RECORDS WHEN 
  #             A USER EMBEDDS THE FUNCTION IN A LOOP TO GENERATE
  #             MULTIPLE SDEs TO THE SAME FILE.
  #             THE filename PARAMETER CONTROLS WHERE THE COORDINATE INFORMATION 
  #             IS WRITTEN TO. USE sdeloc (coordinates) and sdeatt (attributes) 
  #             TO PRODUCE SHAPEFILES USING THE CONVERT.TO.SHAPEFILE AND WRITE.SHAPEFILE 
  #             FUNCTIONS FROM THE SHAPEFILES LIBRARY.
  #
  #  ERROR:     1000  NO ERRORS DETECTED
  #               60  TOO MANY COLUMNS IN DESTINATION MATRIX
  #               61  TOO FEW COLUMNS IN DESTINATION MATRIX
  #               21  INVALID COMBINATION: calcentre=TRUE and centre.xy!=NULL
  #
  #  OUTPUT:	
  #     ID  			   UNIQUE ELLIPSE IDENTIFIER
  #		CALCENTRE		   T|F SHOULD THE MEAN CENTRE BE USED
  #		Weighted     	   T|F SHOULD THE CENTRE BE WEIGHTED (WMC)
  #		CENTRE.x		   X-COORDINATE AT CENTRE
  #		CENTRE.y		   Y-COORDINATE AT CENTRE
  #		Sigma.x			   SIGMA X
  #		Sigma.y			   SIGMA Y
  #		Major			   WHICH AXIS IS THE MAJOR AXIS
  #		Minor			   WHICH AXIS IS THE MINOR AXIS
  #		Theta			   CLOCKWISE ROTATION ANGLE FROM NORTH TO SIGMA Y
  #		Eccentricity	   MEASURE OF SDE ELLIPSE ECCENTRICITY
  #		Area.sde		   AREA OF THE SDE ELLIPSE
  #		TanTheta		   TRIGONOMETRIC RESULT
  #		SinTheta		   TRIGONOMETRIC RESULT
  #		CosTheta		   TRIGONOMETRIC RESULT
  #		SinThetaCosTheta   TRIGONOMETRIC RESULT
  #		Sin2Theta		   TRIGONOMETRIC RESULT
  #		Cos2Theta		   TRIGONOMETRIC RESULT
  #		ThetaCorr		   CLOCKWISE ROTATION ANGLE FROM NORTH TO MAJOR AXIS
  #		sdeatt			   ATTRIBUTES ABOVE WRITTEN TO DATAFRAME FOR POST-PROCESSING AS SHAPEFILE
  #		sdeloc			   UNIQUE ID AND X,Y COORDINATES OF VERTICES FOR POST-PROCESSING INTO SDE SHAPEFILE
  #=======================================================
  
  # INITIALIZE ERROR CODE TO NO ERROR
  errorcode <- 1000
    
  if(length(dim(points)) != 2) {
    # ERROR: TOO FEW COLUMNS IN POINTS COORDINATE MATRIX
    # SET DESCRIPTIVE ERROR CODE AND GIVE WARNING
    errorcode <- 61
    cat("\n\nWARNING: Provided points input matrix has fewer than 2 columns.")
    cat("\nERROR CODE: ", errorcode, "\n\n", sep="")
    return("ERROR")    
  }
  
  if(dim(points)[2] != 2) {
    # ERROR: TOO MANY COLUMNS IN POINTS COORDINATE MATRIX
    # SET DESCRIPTIVE ERROR CODE AND GIVE WARNING
    errorcode <- 60
    cat("\n\nWARNING: Provided points input matrix has too many columns, only 2 are allowed.")
    cat("\nERROR CODE: ", errorcode, "\n\n", sep="")
    return("ERROR")
  }
  else {
    
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
		
		# COMPUTE AND USE WEIGHTED MEAN CENTRE RATHER THAN USER SPECIFIED LOCATION AS CENTRE (WEIGHTED)
        WMC.x <- c( sum(wt.x) / sum(weights) )
        WMC.y <- c( sum(wt.y) / sum(weights) )    
        centre.xy[1] <- WMC.x
        centre.xy[2] <- WMC.y
       }
	  else {
        # COMPUTE AND USE MEAN CENTRE RATHER THAN USER SPECIFIED LOCATION AS CENTRE (NON-WEIGHTED)
        meanx <- sum(points[,1])/n
        meany <- sum(points[,2])/n
        centre.xy[1] <- meanx
        centre.xy[2] <- meany
	  } 
	}
  }
}
  
  # ADD COLUMNS TO points FOR SQUARED x,y TERMS
  points <- cbind(points, points[,1]^2, points[,2]^2)
  
  # ADD COLUMNS TO points FOR TRANSPOSED TERMS
  points <- cbind(points, points[,1]-centre.xy[1], points[,2]-centre.xy[2])
  
  # ADD COLUMNS FOR SQUARED TRANSPOSED TERMS AND PRODUCT OF TRANSPOSED TERMS
  points <- cbind(points, points[,5]^2, points[,6]^2, points[,5]*points[,6])
    
  # ADD COLUMN NAMES TO points
  names(points) <- c("x", "y", "x2", "y2", "x'", "y'", "x'2", "y'2", "x'y'")

  if(weighted) {
  # COMPUTE WEIGHTED THETA (as in EBDON, 1985)
  top1 <- sum(weights*points[,7]) - sum(weights*points[,8])
  top2 <- sqrt( (sum(weights*points[,7])-sum(weights*points[,8]))^2 + 4*(sum(weights*points[,9]))^2 )
  bottom <- (2 * sum(weights*points[,9]))
  tantheta <- (top1 + top2 ) / bottom
  }
  else {
  # COMPUTE THETA (as in EBDON, 1985)
  top1 <- sum(points[,7]) - sum(points[,8])
  top2 <- sqrt( (sum(points[,7])-sum(points[,8]))^2 + 4*(sum(points[,9]))^2 )
  bottom <- (2 * sum(points[,9]))
  tantheta <- (top1 + top2 ) / bottom
  }  
  
  # IF tantheta IS NEGATIVE, ADD 180 TO INVERSE TANGENT OF TANTHETA IN DEGREES
  # TO OBTAIN THE PROPER CLOCKWISE ROTATION ANGLE FROM THE TRANSPOSED AXES
  if(tantheta < 0) {
    theta <- 180 + (atan_d(tantheta)) 
  }
  else {
    theta <- atan_d(tantheta) 
  }
    
  # COMPUTE OTHER TRIGONOMETRIC VALUES
  sintheta <- sin_d(theta)
  costheta <- cos_d(theta)
  sin2theta <- sintheta^2
  cos2theta <- costheta^2
  sinthetacostheta <- sintheta * costheta
  
  # COMPUTE THE STANDARD DEVIATIONS FOR THE TWO NEW AXES OF THE ELLIPSE
  # NOTE THAT EQUATIONS ARE FROM CRIMESTAT (CHAPTER 4) RATHER THAN EBDON (1985)
  # TO ENSURE THE PROPER SIZE FOR THE ELLIPSE

  if(weighted) {
  #PERFORM THE WEIGHTED STANDARD DEVIATIONAL ELLIPSE COMPUTATION (WEIGHTED SDE)
  sigmax <- sqrt(2) * sqrt( ( (sum(weights*points[,7]))*(cos2theta) - 2*(sum(weights*points[,9]))*(sinthetacostheta) + (sum(weights*points[,8]))*(sin2theta) ) / ((sum(weights)) - 2) )
  sigmay <- sqrt(2) * sqrt( ( (sum(weights*points[,7]))*(sin2theta) + 2*(sum(weights*points[,9]))*(sinthetacostheta) + (sum(weights*points[,8]))*(cos2theta) ) / ((sum(weights)) - 2) )
  }
  else {
  #PERFORM THE STANDARD DEVIATIONAL ELLIPSE COMPUTATION (UNWEIGHTED SDE)
  sigmax <- sqrt(2) * sqrt( ( (sum(points[,7]))*(cos2theta) - 2*(sum(points[,9]))*(sinthetacostheta) + (sum(points[,8]))*(sin2theta) ) / (n - 2) )
  sigmay <- sqrt(2) * sqrt( ( (sum(points[,7]))*(sin2theta) + 2*(sum(points[,9]))*(sinthetacostheta) + (sum(points[,8]))*(cos2theta) ) / (n - 2) )
  }  
 
  # CREATE VARIABLES TO HOLD THE IDENTIES OF THE MAJOR AND MINOR AXES FOR OUTPUT
  if(sigmax > sigmay) {
    Major <- "SigmaX"
    Minor <- "SigmaY"
  }
  else {
    Major <- "SigmaY"
    Minor <- "SigmaX"
  }
  
  # STORE THE MAJOR AND MINOR AXIS LENGTHS
  lengthsigmax <- 2 * sigmax
  lengthsigmay <- 2 * sigmay
  
  # STORE THE AREA OF THE SDE
  areaSDE <- pi * sigmax * sigmay
     
  # STORE THE ECCENTRICITY OF THE SDE
  eccentricity <- sqrt(1 - ((min(sigmax,sigmay)^2)/(max(sigmax,sigmay)^2)))

  
  # COMPUTE AND STORE COORDINATES FOR PLOTTING THE SDE ELLIPSE (BASED ON ELLIPSEPOINTS FUNCTION FROM SFSMISC LIBRARY) 
    B <- min(sigmax, sigmay)
    A <- max(sigmax, sigmay)
    d2 <- (A - B) * (A + B)
    phi <- 2 * pi * seq(0, 1, len = 360)
    sp <- sin(phi)
    cp <- cos(phi)
    r <- sigmax * sigmay/sqrt(B^2 + d2 * sp^2)
	xy <- r * cbind(cp, sp)
    al <- (90-theta) * pi/180
    ca <- cos(al)
    sa <- sin(al)
    coordsSDE <- xy %*% rbind(c(ca, sa), c(-sa, ca)) + cbind(rep(centre.xy[1], 360), rep(centre.xy[2], 360))  
      
  if(verbose) {
    cat("\n----------------------------------------------")
    cat("\nCoordinates of centre (x): ", centre.xy[1], sep="")
    cat("\nCoordinates of centre (y): ", centre.xy[2], sep="")
    cat("\nAngle of rotation: ", round(theta, 2), " clockwise degrees", sep="")
    cat("\nLength of X axis: ", round(lengthsigmax, 2), sep="")      
    cat("\nLength of Y axis: ", round(lengthsigmay, 2), sep="")      
    cat("\nArea of SDE ellipse: ", round(areaSDE, 2), sep="")
    cat("\ntantheta: ", tantheta, sep="")      
    cat("\ntheta: ", theta, sep="")      
    cat("\nsintheta: ", sintheta, sep="")      
    cat("\ncostheta: ", costheta, sep="")      
    cat("\nsinthetacostheta: ", sinthetacostheta, sep="")      
    cat("\nsin2theta: ", sin2theta, sep="")      
    cat("\ncos2theta: ", cos2theta, sep="")                                          
    cat("\nsigmax: ", sigmax, sep="")                                          
    cat("\nsigmay: ", sigmay, sep="") 
    cat("\neccentricity: ", eccentricity, sep="")        
    cat("\n----------------------------------------------\n\n")
  }  
    
  # COMPUTE Theta.Corr WHICH IS A CLOCKWISE ANGLE OF ROTATION FOR THE MAJOR AXIS
  # RATHER THAN THE ANGLE TO SIGMA.Y WHICH IS CURRENTLY COMPUTED
  if(sigmax < sigmay) {
    Theta.Corr <- theta
  }
  else {
    Theta.Corr <- theta + 90
  }

  # STORE RESULTS INTO A LIST (REQUIRED FOR PLOT FUNCTION)
  r.SDE <- list(id = id, points = points, coordsSDE = coordsSDE, calccentre = calccentre, CENTRE.x = centre.xy[1], CENTRE.y = centre.xy[2], 
                Major = Major, Minor = Minor, theta = theta, Sigma.x = sigmax, Sigma.y = sigmay, Eccentricity = eccentricity, 
				Area.sde = areaSDE, TanTheta = tantheta, SinTheta = sintheta, CosTheta = costheta, SinThetaCosTheta = sinthetacostheta, 
				Sin2Theta = sin2theta, Cos2Theta = cos2theta, ThetaCorr = Theta.Corr, weighted = weighted, weights = weights)
  assign("r.SDE", r.SDE, pos=1)
	
  # STORE SDE ATTRIBUTES INTO A DATA FRAME AND PRINTS RESULTS
  result.sde <- list("id"=id, "CALCCENTRE"=calccentre, "weighted"=weighted,
                     "CENTRE.x"=centre.xy[1], "CENTRE.y"=centre.xy[2], "Sigma.x"=sigmax, "Sigma.y"=sigmay,
                     "Major"=Major, "Minor"=Minor, "Theta"=theta, "Eccentricity"=eccentricity, "Area.sde"=areaSDE, 
                     "TanTheta"=tantheta, "SinTheta"=sintheta, "CosTheta"=costheta, "SinThetaCosTheta"=sinthetacostheta,
                     "Sin2Theta"=sin2theta, "Cos2Theta"=cos2theta, "ThetaCorr"=Theta.Corr)
  print(result.sde)
  result.sde<-as.data.frame(result.sde)

  # DATA FRAME OF ATTRIBUTES WITH FIRST COLUMN NAME "ID" FOR CONVERT.TO.SHAPEFILE FUNCTION
  assign("sdeatt", result.sde, pos=1)	

  # CREATE ASCII OUTPUT FOR SHAPEFILE CREATION
  sdeloc <- as.data.frame(cbind(id, coordsSDE))
  colnames(sdeloc)=c("id","x","y")
  write.table(sdeloc, sep=",", file=filename, col.names=FALSE)
  
  # DATA FRAME WITH COLUMNS IN ORDER ID, X-COORD, Y-COORD FOR CONVERT.TO.SHAPEFILE FUNCTION
  assign("sdeloc", sdeloc, pos=1)	
}


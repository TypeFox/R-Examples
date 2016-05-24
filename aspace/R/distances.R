"distances" <-
function(centre.xy=centre, destmat=activities, verbose=FALSE) {
  
  #=======================================================
  #
  #  TITLE:     MULTIPLE EUCLIDEAN DISTANCE CALCULATOR
  #  FUNCTION:  distances()
  #  AUTHOR:    TARMO K. REMMEL 
  #  DATE:      27 MARCH 2006
  #  CALLS:     NA
  #  NEEDS:     A VECTOR INPUT HOLDING THE EASTING AND
  #             NORTHING COORDINATES OF A CENTRE
  #  NOTES:     RETURNS A VECTOR OF DISTANCES TO EACH
  #               ACTIVITY LOCATION FROM THE SPECIFIED CENTRE
  #             COLUMN 1 IS EASTING (x) IN UTM COORDINATES
  #             COLUMN 2 is NORTHING (y) IN UTM COORDINATES
  #  ERROR:     1000  NO ERRORS DETECTED
  #               50  TOO MANY COLUMNS IN DESTINATION MATRIX
  #               51  TOO FEW COLUMNS IN DESTINATION MATRIX
  #
  #=======================================================

  # GIVEN A MATRIX OF EASTING AND NORTHING UTM COORDINATES FOR A SET OF POINT
  # LOCATIONS, COMPUTE THE EUCLIDEAN DISTANCE FROM A SELECTED 
  # CENTRE TO EACH POINT LOCATION INDEPENDENTLY - IMPORTANT NOTE: THE FUNCTION ASSUMES PLANAR COORDINATES (UNIVERSE TRANSVERSE MERCATOR)
  
  # INITIALIZE ERROR CODE TO NO ERROR
  errorcode <- 1000
  
  if(length(dim(destmat)) != 2) {
    # ERROR: TOO FEW COLUMNS IN DESTINATION COORDINATE MATRIX
    # SET DESCRIPTIVE ERROR CODE AND GIVE WARNING
    errorcode <- 51
    cat("\n\nWARNING: Provided activity locations input matrix has fewer than 2 columns.")
    cat("\nERROR CODE: ", errorcode, "\n\n", sep="")
    return("ERROR")    
  }
  
  if(dim(destmat)[2] != 2) {
    # ERROR: TOO MANY COLUMNS IN DESTINATION COORDINATE MATRIX
    # SET DESCRIPTIVE ERROR CODE AND GIVE WARNING
    errorcode <- 50
    cat("\n\nWARNING: Provided activity locations input matrix has too many columns, only 2 are allowed.")
    cat("\nERROR CODE: ", errorcode, "\n\n", sep="")
    return("ERROR")
  }
  else {
    # COMPUTE THE EUCLIDEAN DISTANCE
    di <- sqrt( (destmat[,1] - centre.xy[1])^2 + (destmat[,2] - centre.xy[2])^2 ) 
  }
  
  # PROVIDE THE VECTOR OF DISTANCES, di, AS A RETURN PARAMETER TO THE CALLING FUNCTION
  return(di)
  
}


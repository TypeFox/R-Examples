"atan_d" <-
function(theta=0) {

  #=======================================================
  #
  #  TITLE:     COMPUTE ARCTANGENT IN DEGREES
  #  FUNCTION:  atan_d()
  #  AUTHOR:    TARMO K. REMMEL 
  #  DATE:      16 JANUARY 2006
  #  CALLS:     NA
  #  NEEDS:     
  #  NOTES:     TO SIMPLIFY THE USE OF RADIANS IN R
  #
  #=======================================================

  return(atan(theta)*180/pi)
}


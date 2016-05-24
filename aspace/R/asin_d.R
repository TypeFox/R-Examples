"asin_d" <-
function(theta=0) {

  #=======================================================
  #
  #  TITLE:     COMPUTE ARCSINE IN DEGREES
  #  FUNCTION:  asin_d()
  #  AUTHOR:    TARMO K. REMMEL 
  #  DATE:      16 JANUARY 2006
  #  CALLS:     NA
  #  NEEDS:     
  #  NOTES:     TO SIMPLIFY THE USE OF RADIANS IN R
  #
  #=======================================================

  return(asin(theta)*180/pi)
}


"acos_d" <-
function(theta=0) {

  #=======================================================
  #
  #  TITLE:     COMPUTE ARCCOSINE IN DEGREES
  #  FUNCTION:  acos_d()
  #  AUTHOR:    TARMO K. REMMEL 
  #  DATE:      16 JANUARY 2006
  #  CALLS:     NA
  #  NEEDS:     
  #  NOTES:     TO SIMPLIFY THE USE OF RADIANS IN R
  #
  #=======================================================

  return(acos(theta)*180/pi)
}


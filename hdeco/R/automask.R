"t.automask" <-
function(IMG=NULL, minvalue=1) {

  #############################################################
  # 
  # TITLE:  t.automask()
  # AUTHOR: TARMO REMMEL
  # DATE:   25 NOV, 2003
  # CALLS:  
  # NEEDS:  AN INPUT IMAGE - RETURNS A MASK IMAGE [0,1]
  # NOTES:  BUILD MASK TO TO EXCLUDE BAD EOSD CATEGORIES
  #			PROCESS UNDER 1s, MASK OUT UNDER 0s
  #			THIS FUNCTION MUST BE CALLED BEFORE THE HDECO SERIES
  #         IF MINVALUE=2, REMOVES ORIGINAL 0s - AS ALL CATEGORIES
  #         WILL BE INCREASED BY ONE TO CONVERT 0s TO 1s.
  #
  #############################################################

  txtstring <- paste("\nBuilding auto-mask (Blocking out values < ", minvalue, ").", sep="")
  cat(txtstring)
  
  # CONSTRUCT A MASK OF 1s THE SAME SIZE AS THE INPUT IMAGE
  MSK <- (IMG * 0) + 1
  
  # IF MINVALUE=3, MASK OUT EOSD VALUES 1 AND 2 (SHADOW, SNOW/ICE WHICH DO NOT EXIST IN NFI)
  # ALSO MASKS OUT ZEROS AND NEGATIVE VALUES THAT CAUSE PROBLEMS
  MSK[IMG < minvalue] <- 0
 
  return(MSK)
}


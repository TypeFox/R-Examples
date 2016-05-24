"bemar" <-
function (ik=2, BE=.QKEP) {

  #############################################################
  # 
  # TITLE:  	bemar()
  # AUTHOR: 	SANDOR KABOS, TARMO REMMEL
  # DATE:   	20 AUG, 2003
  # CALLS: 		 
  # CALLED BY:	tesz()
  # NEEDS:  	
  # NOTES:  	SUM OVER ALL VARIABLES BY NOT THE ik-TH ONE
  #             ONE SUM FOR EACH LEVEL OF THE VARIABLES INDICATED
  #
  #############################################################
  
  # STORE THE NUMBER OF VARIABLES IN .QKEP
  NC <- length(dim(BE))

  # IF ik IS EQUAL TO NC, RETURN THE ENTIRE .QKEP MATRIX
  if(length(ik) == NC) {
    return(BE)
  }
 
  # SUM .QKEP MATRIX OVER ALL VARIABLES BASED ON GIVEN ik
  A <- c(1:NC)
  KI <- apply(BE, A[ik], sum)
  return(KI)
  
}


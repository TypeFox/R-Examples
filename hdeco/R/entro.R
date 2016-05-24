"entro" <-
function (BE=.QKEP) {

  #############################################################
  # 
  # TITLE:  entro()
  # AUTHOR: SANDOR KABOS (MODIFIED BY TARMO REMMEL)
  # DATE:   20 AUG, 2003
  # CALLS:  
  # NEEDS:  .QKEP
  # NOTES:  COMPUTE ENTROPY USING BASE 2 LOGARITHM
  #
  #############################################################

  # COMPUTE ONLY FOR POSITIVE VALUES IN THE INPUT BE
  BE <- as.vector(BE)
  HOL <- BE > 0
  BE <- BE[HOL]

  # COMPUTE ENTROPY
  KI <- -sum(BE * log(BE, 2))

  # RETURN THE ENTROPY VALUE
  return(KI)

}

"kv" <-
function (BE=c(3,4,1)) {

  #############################################################
  # 
  # TITLE:  	kv()
  # AUTHOR: 	SANDOR KABOS, MODIFIED BY TARMO REMMEL
  # DATE:   	26 NOV, 2003
  # CALLS:  	
  # CALLED BY:	benya()
  # NEEDS:  	
  # NOTES:  	
  #
  #############################################################

  Q <- c(0,0,1,1,0,1,0,1)
  dim(Q) <- c(4,2)

  H <- length(BE)
  KI <- Q[BE[1],]
  if(H != 1){
    for (h  in 2:H){
      KI <- 2 * KI + Q[BE[h],]
    }
  }
  return(KI + 1)
}


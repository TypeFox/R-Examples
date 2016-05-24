"benya" <-
function (BE=KEP1) {

  #############################################################
  # 
  # TITLE:  	benya()
  # AUTHOR: 	SANDOR KABOS, MODIFIED BY TARMO REMMEL
  # DATE:   	26 NOV, 2003
  # CALLS:  	kv(), nalap()
  # CALLED BY:	benc4(), bend4(), bendManyX4()
  # NEEDS:  	
  # NOTES:  	CONVERTS AN IMAGE MATRIX INTO COLUMNAR FORM
  #	        SUCH THAT QKEP MATRIX CAN BE BUILT
  #
  #############################################################
  
  H <- dim(BE)[1]
  L <- trunc(log(H, 2))

  if(H!=dim(BE)[2]) {
    return(cat("\nImage error: Dimensions do not match."))
  }

  H2 <- H^2
  KI <- matrix(nrow=H2, ncol=L+1)

  for(h in 1:H2) {
    N <- nalap(h-1, 4, L) + 1
    KV <- kv(N)
    K <- KV[1]
    V <- KV[2]
    KI[h,] <- c(N,BE[K,V])
  }
  return(KI)
}


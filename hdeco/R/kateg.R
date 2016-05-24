"kateg" <-
function (BE=BEMM) {

  #############################################################
  # 
  # TITLE:  kateg()
  # AUTHOR: SANDOR KABOS (MODIFIED BY TARMO REMMEL)
  # DATE:   2003
  # CALLS:  kossze(), kvissza(), kegy()
  # NEEDS:  
  # NOTES:  USED TO CREATE NESTED CFS CODING STRUCTURE
  #		THIS IS LATER ADDED TO, CREATING .LUT
  #
  #############################################################

  B <- BE
  N <- dim(BE)[1]
  H <- dim(BE)[2]
  KI <- matrix(0, nrow=N, ncol=H, dimnames=list(NULL, paste("KAT", 1:H, sep="")))

  for(h in 1:(H-1)){
    B <- kossze(h, BE)
    KI[,h] <- kegy(B[,h+1])[kvissza(BE[,h])]
  }

  KI[,H] <- kvissza(BE[,H])
  assign(".KATEG", KI, pos=1)
}

"benci" <-
function (BE=MASZK) {
  
  #############################################################
  # 
  # TITLE:  benci()
  # AUTHOR: SANDOR KABOS
  # DATE:   15 AUG, 2003
  # CALLS:  
  # NEEDS:  MASK INPUT AS MATRIX {0,1} - PROCESS UNDER 1s
  # NOTES:  READS THE MASK IMAGE INTO A MULTIDIMENSIONAL ARRAY
  #
  #############################################################
  
  MIK <- names(table(BE))

  if(!prod(MIK == c("0","1"))) {
    cat("\nWarning: MASK might not contain 0s or 1s.")
  }

  KI <- benya(BE)
  NC <- dim(KI)[2]
  return(KI[,NC])
}


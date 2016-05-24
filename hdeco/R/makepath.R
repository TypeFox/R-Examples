"makepath" <-
function(NX=2, NY=6, NZ=2, STEPS=7, TITLE="Custom Decomposition Path") {

  #############################################################
  # 
  # TITLE:  makepath()
  # AUTHOR: TARMO REMMEL
  # DATE:   27 NOV, 2003
  # CALLS:  
  # NEEDS:  
  # NOTES:  BUILDS VFONAL, A CUSTOM DECOMPOSITION PATH MATRIX
  #	      PARAMETERS ARE THE NUMBER OF X, Y, Z VARIABLES
  #	      AND THE NUMBER OF DECOMPOSITION STEPS ALONG WITH 
  #	      A TITLE FOR THE DECOMPOSITION MATRIX
  #	      WORKS IN INTERACTIVE MODE
  #
  #############################################################
  
  # BUILD THE RAW STRUCTURE OF THE DECOMPOSITION PATH MATRIX AND FILL WITH -1 IN EACH CELL
  columns <- c(paste("X", 1:NX, sep=""), paste("Y", 1:NY, sep=""), paste("Z", 1:NZ, sep=""))
  VFONAL <- matrix(-1, nrow=STEPS, ncol=(NX+NY+NZ), dimnames=list(NULL,columns))

  # STARTS AN INTERACTIVE SPREADSHEET-STYLE EDITOR TO FILL CELL VALUES
  fix(VFONAL)
  
  # APPLY DEFAULT TITLE
  attr(VFONAL,"cim") <- TITLE
}


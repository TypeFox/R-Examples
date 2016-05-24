"vfon" <-
function (BE=.VFONAL) {
  
  #############################################################
  # 
  # TITLE:  vfon()
  # AUTHOR: SANDOR KABOS, MODIFIED BY TARMO REMMEL
  # DATE:   26 NOV, 2003
  # CALLS:  
  # CALLED BY:	vfonal()
  # NEEDS:  	.VFONAL (THE DECOMPOSITION PATH)
  # NOTES:  	CREATES .AHIPO, .NHIPO, AND .KIVALO OBJECTS
  #               FROM THE DECOMPOSITION PATH MATRIX
  #
  #############################################################

  NEVEK <- dimnames(BE)[[2]]
  NN <- dim(BE)[1]
  HH <- length(NEVEK)
  HHH <- 1:HH
  names(HHH) <- NEVEK
  NHIPO <- list()
  AHIPO <- list()
  KIVALO <- list()

  for(n in 1:NN){
    BEE <- BE[n,]
    AHIP <- HHH[BEE > 0]
    NHIP <- HHH[BEE == 0]
    KIVA <- HHH[BEE >= 2]
    
    if(length(AHIP) == 0){
      cat("\n\nThere must be at least one 1 in PathDecoMatrix (Step ",n,")\n", sep="")
    }
    if(length(NHIP) == 0){
      cat("\n\nThere must be at least one 0 in PathDecoMatrix (Step ",n,")\n", sep="")
    }
    if(length(KIVA) == 0){
      KIVA <- NULL
    }
    
    NHIPO <- c(NHIPO,list(NHIP))
    AHIPO <- c(AHIPO,list(AHIP))
    KIVALO <- c(KIVALO,list(KIVA))
  }

  # CREATE LOCAL COPIES OF ALTERNATE AND NULL HYPOTHESIS MATRICES
  assign(".AHIPO", AHIPO, pos=1)
  assign(".NHIPO", NHIPO, pos=1)
  assign(".KIVALO", KIVALO, pos=1)
  return(cat("."))
}


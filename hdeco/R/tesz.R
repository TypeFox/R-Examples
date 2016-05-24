"tesz" <-
function (ik=3, QKEP=.QKEP) {

  #############################################################
  # 
  # TITLE:      tesz()
  # AUTHOR:     SANDOR KABOS, MODIFIED BY TARMO REMMEL
  # DATE:   	27 NOV, 2003
  # CALLS:  	entro(), bemar()
  # CALLED BY: 	teszt()
  # NEEDS:  
  # NOTES:  
  #
  #############################################################

  # STORE .QND TO MODELL AND THE COLUMN NAMES TO NEVEK
  MODELL <- .QND
  NEVEK <- names(MODELL)

  # REPLACE MODELL VALUES WITH 1 TO THE LENGTH OF .QND AND ADD THE NAMES BACK
  MODELL <- 1:length(NEVEK)
  names(MODELL) <- NEVEK

  # READ LOCAL VARIABLES FOR GIVEN STEP
  AHIPO <- .AHIPO[[ik]]
  NHIPO <- .NHIPO[[ik]]
  KIVALO <- .KIVALO[[ik]]

  # ERROR CHECKING FOR ALTERNATE HYPOTHESES
  if(prod(AHIPO %in% MODELL)!=1) {
    cat("AHIPO error: Step ",ik,"\n")
    return()
  }

  # ERROR CHECKING FOR NULL HYPOTHESES
  if(prod(NHIPO %in% MODELL)!=1) {
    cat("NHIPO error: Step ",ik,"\n")
    return()
  }

  ALAPF <- sort(union(AHIPO,NHIPO))
  names(ALAPF) <- NEVEK[ALAPF]
  BEVALO.ALAPF <- sort(setdiff(ALAPF,KIVALO))
  BEVALO.AHIPO <- sort(setdiff(AHIPO,KIVALO))
  names(BEVALO.ALAPF) <- NEVEK[BEVALO.ALAPF]
  names(BEVALO.AHIPO) <- NEVEK[BEVALO.AHIPO]

  # COMPUTE ENTROPY OF NULL, JOINT, AND ALTERNATE HYPOTHESES
  HNULLHIPO <- entro(bemar(NHIPO))
  HJOINT <- entro(bemar(ALAPF))
  HAHIPO <- entro(bemar(AHIPO))
  HBEVALO.JOINT <- entro(bemar(BEVALO.ALAPF))
  HBEVALO.AHIPO <- entro(bemar(BEVALO.AHIPO))
  assign(".BASE",ALAPF,pos=1)
  KI <- list(HALAPF=HJOINT,HBEVALO.ALAPF=HBEVALO.JOINT,HHIPO=HNULLHIPO,HAHIPO=HAHIPO,HBEVALO.AHIPO=HBEVALO.AHIPO,BEVALO.AHIPO=BEVALO.AHIPO,ALAPF=ALAPF,BEVALO.ALAPF=BEVALO.ALAPF)
  return(KI)
}


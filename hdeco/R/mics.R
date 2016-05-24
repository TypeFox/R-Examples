"mics" <-
function () {

  #########################################################################
  #
  # TITLE:     mics()
  # AUTHOR:    SANDOR KABOS
  # DATE:      ?
  # CALLS:     
  # NEEDS:     
  # NOTES:     AUTOMATIC CREATION OF A DECOMPOSITION PATH MATRIX
  #
  # SAMPLE:    
  #
  #########################################################################

  KI <- list()
  DIM <- .QND
  NEVEK <- names(DIM)
  EGYKE <- DIM[1]
  
  if(EGYKE!=1) EXEK <- 1 else EXEK <- 0
 
  YEK <- "Y"==substring(NEVEK,1,1)
  NY <- sum(YEK)
  VN <- "Y"
  
  for (i in 1:NY){
    MIK <- i
    KUKK <- list(VARIAB=VN,WHICH=MIK,EXE=EXEK)
    KI <- c(KI,KUKK)
  }

  attr(KI,"cim") <- "Default Decomposition Path"
  return(KI)

}

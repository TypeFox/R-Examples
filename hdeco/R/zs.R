"zs" <-
function (KIK=c("Z"==substring(names(.QND),1,1)),EGYKE=1) {

  # BUILDS .ZS OBJECT
  # CALLS zsz(), entro(), bemar()
  # KIK IS A VECTOR T/F INDICATING WHICH .QND VARIABLES ARE Z VARIABLES
  
  cat("\nIn zs()")
  # STORE THE NUMBER OF ALTERNATE HYPOTHESES
  H <- length(.AHIPO)
  # STORE THE NUMBER OF LEVELS IN THE Z VARIABLES
  N <- prod(.QND[KIK])
  KI <- list()

  # DEPENDING ON THE NUMBER OF IMAGES...
  if(EGYKE == 1) {
    cat("\nIn zs() 1 image")  
    for(h in 1:H) {
      EZ <- zsz(h)

      if(is.null(EZ)) {
        KI[[h]] <- NULL
      }
      else {
        EZ[,1] <- EZ[,1] - EZ[,2]
        KI[[h]] <- EZ[,-2]
      }

    }
  } 
  else {
    cat("\nIn zs() 2 images")  
    # CREATE VECTOR T/F OF X VARIABLES
    XEK <- ("X"==substring(names(.QND),1,1))
    # STORE THE NUMBER OF X VARIABLES
    KEK <- (1:length(XEK))[XEK]
    # ENTROPY OF X VARIABLES
    XENTRO <- entro(bemar(KEK))

    # REPEAT FOR EACH ALTERNATE HYPOTHESIS
    for(h in 1:H) {  
      cat("\nIn Loop calling zsz()")
      # PROBLEM WITH NEXT LINE ***
      EZ <- zsz(h)
      AZ <- zsz(h,T)
      
      cat("\nHere.")

      if(is.null(EZ)|is.null(AZ)) {
        KI[[h]] <- NA
      }
      else {
        EZ[,1] <- AZ[,2]-EZ[,2] + XENTRO
        KI[[h]] <- EZ[,-2]
      }

    }
  }
  assign(".ZS",KI,pos=1)
  cat("\nLeaving zs()")
}


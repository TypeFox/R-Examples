"nalap" <-
function (BE=57,alap=4,BALO=5) {

  KI <- NULL
 
  if(BE < alap){
    KI <- BE
    if(BALO > 1){
      KI <- rep(0,BALO)
      KI[BALO] <- BE
    } 
    return(KI)
  }

  for (i in 1:999){
    M <- BE%%alap
    KI <- c(KI,M)
    BE <- BE%/%alap
    if(BE < alap){
      KI<-c(BE,KI)
      BH <- BALO - length(KI)
      if(BH >0)
        KI <- c(rep(0,BH),KI) 
    return(KI)
    }
  }
}


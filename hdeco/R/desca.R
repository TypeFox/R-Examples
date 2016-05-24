"desca" <-
function (MIT=c(0,1,2),MIHEZ=MIT,HANYSZOR=2) {

  if(is.null(MIHEZ)) {
    return(MIT)
  }

  if(length(dim(MIHEZ))<2) {
    MIHEZ <- matrix(MIHEZ,ncol=1)
  }

  HANYSZOR <- HANYSZOR - 1
  KI <- MIHEZ
  for(h in 1:HANYSZOR) {
    N <- dim(KI)[1]
    K <- NULL
    for (m in MIT){
      K <- rbind(K,cbind(KI,rep(m,N)))
    }
    KI <- K
  }

  NC <- dim(KI)[2]
  KI <- KI[,rev(1:NC)]
  
  return(KI)

}


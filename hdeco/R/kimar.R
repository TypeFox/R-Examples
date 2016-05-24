"kimar" <-
function (ik=2,BE=.QKEP) {

  #SUM OVER THE ik-TH VARIABLES

  NC <- length(dim(BE))
  A <- c(1:NC)
  KI <- apply(BE,A[-ik],sum)
  return(KI)
}


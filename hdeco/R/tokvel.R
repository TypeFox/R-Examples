"tokvel" <-
function (level=6, NSZIN=4, epa=FALSE) {
  M <- 2^level
  N <- M^2
  KI <- sample(1:NSZIN, N, T)
  KI <- matrix(KI, nrow=M, ncol=M)
  attr(KI,"cim") <- "tokveletlen kep"
  if(epa) KI <- KI - 1
  return(KI)
}

"kikap" <-
function (BE=.QKEP,NDX=1,NDY=3,NDZ=2) {
  N <- dim(BE)[1]
  NDD <- NDX + NDY
  ND <- NDD + NDZ
  MAX <- apply(BE, 2, max)
  KI <- array(0, dim=c(MAX[1:NDD], rep(2,NDZ)))
  KI[BE] <- 1
  return(KI)
}


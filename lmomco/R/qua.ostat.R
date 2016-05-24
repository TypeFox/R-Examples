"qua.ostat" <-
function(f,j,n,para=NULL) {
  betainv <- qbeta(f,j,n-j+1)
  return(par2qua(betainv,para))
}

# Assumes full column rank in the QR decomposition

backsolveSparse <- function(QR, b) {
  #R = qrR(QR)
  #x = solve(R, qr.qty(QR,b)[Seq(1,nrow(R))])
  #return(as.numeric(x))
  R = qr.R(QR)
  x = solve(R, qr.qty(QR,b)[Seq(1,nrow(R))])
  if (length(QR@q)==0) return(x)
  else return(x[Order(QR@q+1)])
}

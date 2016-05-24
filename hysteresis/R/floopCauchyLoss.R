floopCauchyLoss <- function (par,x,y) {
  t <- par[1:length(x)]
  cx <- par[length(x)+1]
  cy <- par[length(x)+2]
  b.x <- par[length(x)+3]
  b.y <- par[length(x)+4]
  logm <- par[length(x)+5]
  logn <- par[length(x)+6]
  retention <- par[length(x)+7]
  times <- cumsum(t)
  costp <- cos(times) 
  sintp <- sin(times) 
  direc <- sign(costp)
  direcsin <- sign(sintp)
  pred.x <- cx+b.x*costp
  pred.y <- cy+direcsin*retention*abs(sintp)^exp(logm)+direc*b.y*abs(costp)^exp(logn)
  logloss <- crossprod(y-pred.y[1:length(y)])+crossprod(x-pred.x)
  logloss
}


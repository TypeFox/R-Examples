coefinitial.ftn <-
function(dat, start){
  y=dat$dat[, 3]
  maxy=max(abs(y))
  sumdfs=sum(dat$dfs)
  A=3*maxy
  B=1
  D1=exp(-B*log(-(1+A/(y))))
  temp=cls(D1, as.matrix(dat$dat[, 4:(sumdfs+4)]) )
  beta=temp$betahat
  theta0=start
  coef=c(A, B, beta, theta0)
  return(coef)
}

logproby <-
function(y,XX,ZZ,aa,bb,s2,nu){
  df=nu
  n=length(y)
  U=aa*ZZ+bb*XX+diag(n)
  S=s2*U
  log.prob=dmt(y,mean=rep(0,n),df=df,S=S,log=TRUE)
  log.prob
}

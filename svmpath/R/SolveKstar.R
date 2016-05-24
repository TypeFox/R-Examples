SolveKstar<-function(Kstar,ridge=1e-10,...){
  onestar<-rep(1,ncol(Kstar))
  onestar[1]<-0
  if(ridge>0) Kstar=Kstar+diag(length(onestar))*ridge
  solve(Kstar,onestar)
}

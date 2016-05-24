## These are functions called by eqp.eq() to perform equipercentile equating

F.x<-function(x,S){
  max(0,min(1,mean(S<=x)))
}
 
P.x<-function(x,S){
  x.s<-floor(x+.5)
  max(0,min(100, 100*(F.x(x.s-1,S)+(x-(x.s-.5))*(F.x(x.s,S)-F.x(x.s-1,S)))))
}

x.s.up<-function(p.s,S){
  val<-S[max(1,floor((p.s/100)*length(S)))]
  while (F.x(val,S)<=(p.s/100)){
    val<-val+1
}
  val
}

x.s.lo<-function(p.s,S){
  val<-S[min(length(S),ceiling((p.s/100)*length(S)))]
  while (F.x(val,S)>=(p.s/100)){
    val<-val-1
}
  val
}

P.inv<-function(p.s,S,Kx=max(S)){
  if (p.s==100){x.up.ps<-Kx+.5}
  else{
    x.s.u<-x.s.up(p.s,S)
    x.up.ps<-(p.s/100-F.x(x.s.u-1,S))/(F.x(x.s.u,S)-F.x(x.s.u-1,S))+(x.s.u-.5)
  }
  if (p.s==0){x.lo.ps<--.5}
  else{
    x.s.l<-x.s.lo(p.s,S)
    x.lo.ps<-(p.s/100-F.x(x.s.l,S))/(F.x(x.s.l+1,S)-F.x(x.s.l,S))+(x.s.l+.5)
  }
  mean(c(x.lo.ps,x.up.ps))
}

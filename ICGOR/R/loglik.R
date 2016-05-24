loglik <-
function(x0,sdata,Xp,r,bl.Li,bl.Ri){
  P<-ncol(Xp)
  b0<-x0[1:P]
  g0<-x0[-c(1:P)]
  Lamd.Li<-t(bl.Li)%*%matrix(g0,ncol=1)
  Lamd.Ri<-t(bl.Ri)%*%matrix(g0,ncol=1)
  exb<-exp(Xp%*%b0)
  part1<-log((1-(1+r*Lamd.Ri*exb)^(-1/r))^sdata$d1)
  part2<-log(((1+r*Lamd.Li*exb)^(-1/r)-(1+r*Lamd.Ri*exb)^(-1/r))^sdata$d2)
  part3<-log((1+r*Lamd.Li*exb)^(-1/r*sdata$d3))
return(sum(part1+part2+part3))
}

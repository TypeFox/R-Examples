bb.gamma0 <-
function(bb,sdata,Xp,Eyil,Ewil,Eui,bl.Li,bl.Ri){
  N<-nrow(sdata)
  L<-nrow(bl.Li)
  r.mat<-function(x) matrix(x,ncol=N,nrow=L,byrow=TRUE)
  c.mat<-function(x) matrix(x,ncol=N,nrow=L)
  exb<-exp(Xp%*%bb)
  gg<-gg.beta0(bb,sdata,Xp,Eyil,Ewil,Eui,bl.Li,bl.Ri)
  prt1<-(c.mat(log(gg))+r.mat(log(exb)))*(Eyil+Ewil)
  prt2<-c.mat(gg)*(r.mat(exb*(1-sdata$d3))*bl.Ri+r.mat(exb*sdata$d3*Eui)*bl.Li)
  Q<-sum(prt1-prt2)
return(-Q)
}

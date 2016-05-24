bb.gamma <-
function(bb,r,sdata,Xp,Eyil,Ewil,Efi,Euifi,bl.Li,bl.Ri){
  N<-nrow(sdata)
  L<-nrow(bl.Li)
  r.mat<-function(x) matrix(x,ncol=N,nrow=L,byrow=TRUE)
  c.mat<-function(x) matrix(x,ncol=N,nrow=L)
  exb<-exp(Xp%*%bb)
  gg<-gg.beta(bb,r,sdata,Xp,Eyil,Ewil,Efi,Euifi,bl.Li,bl.Ri)
  prt1<-(Eyil+Ewil)*(c.mat(log(gg))+r.mat(log(exb)))
  prt2<-r*c.mat(gg)*(r.mat(exb*Efi*(1-sdata$d3))*bl.Ri+r.mat(exb*Euifi*sdata$d3)*bl.Li)
  Q<-sum(prt1-prt2)
return(-Q)
}

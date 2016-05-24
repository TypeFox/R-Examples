gg.beta <-
function(bb,sdata,Xp,Eyil,Ewil,Efi,r,bl.Li,bl.Ri){
  N<-nrow(sdata)
  L<-nrow(bl.Li)
  r.mat<-function(x) matrix(x,ncol=N,nrow=L,byrow=TRUE)
  c.mat<-function(x) matrix(x,ncol=N,nrow=L)
  numer<-apply(r.mat(sdata$d1)*Eyil+r.mat(sdata$d2)*Ewil,1,sum)
  exb<-exp(Xp%*%bb)
  denom<-apply(r.mat(r*exb*Efi)*(r.mat(1-sdata$d3)*bl.Ri+r.mat(sdata$d3)*bl.Li),1,sum)
  return(numer/denom)
}

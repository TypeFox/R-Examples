gg.beta0 <-
function(bb,sdata,Xp,Eyil,Ewil,Eui,bl.Li,bl.Ri){
  N<-nrow(sdata)
  L<-nrow(bl.Li)
  r.mat<-function(x) matrix(x,ncol=N,nrow=L,byrow=TRUE)
  c.mat<-function(x) matrix(x,ncol=N,nrow=L)
  exb<-exp(Xp%*%bb)
  numer<-apply(Eyil+Ewil,1,sum)
  denom<-apply(r.mat((1-sdata$d3)*exb)*bl.Ri+r.mat(sdata$d3*exb*Eui)*bl.Li,1,sum)
  gg<-numer/denom
return(gg)
}

bb.gamma <-
function(bb,sdata,Xp,Eyil,Ewil,Efi,r,bl.Li,bl.Ri){
  N<-nrow(sdata)
  L<-nrow(bl.Li)
  r.mat<-function(x) matrix(x,ncol=N,nrow=L,byrow=TRUE)
  c.mat<-function(x) matrix(x,ncol=N,nrow=L)
  gg<-gg.beta(bb,sdata,Xp,Eyil,Ewil,Efi,r,bl.Li,bl.Ri)
  part1<-(c.mat(log(gg))+r.mat(Xp%*%bb))*(r.mat(sdata$d1)*Eyil+r.mat(sdata$d2)*Ewil)
  part2<-r*c.mat(gg)*r.mat(exp(Xp%*%bb)*Efi)*(r.mat(1-sdata$d3)*bl.Ri+r.mat(sdata$d3)*bl.Li)
  Q<-sum(part1-part2)
return(-Q)
}

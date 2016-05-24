inters<-function(x, y) { 
  # Intersection of propositions
  # 
  I1<-rbind(x[,-1])
  I2<-rbind(y[,-1])
  N1<-array(c(I1),c(dim(I1),dim(I2)[1]))
  N2a<-array(c(I2),c(dim(I2),dim(I1)[1]))
  N2<-aperm(N2a,c(3,2,1))
  N12<-N1*N2 
}
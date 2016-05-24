"sym.P" <- function(x){
  alpha<-x[upper.tri(x)]
  x1<-t(x)
  beta<-x1[upper.tri(x1)]
  gamma<-(alpha+beta)/2
  x[upper.tri(x)]<-gamma
  x2<-t(x)
  x[lower.tri(x)]<-x2[lower.tri(x2)]
  x
}

near.bound<-function(X,Y,bX,bY){
  euc<-function(x1,x2,y1,y2)sqrt((x1-x2)^2+(y1-y2)^2)##euclidean distance
n<-length(X)
nb<-length(bX)
e<-matrix(ncol=n,nrow=nb)
  for(i in 1:n){
  e[,i]<-euc(X[i],bX,bY,Y[i])
  }
mn<-apply(e,2,function(x){which(x==min(x))})
xy.bound<-data.frame(near.X.coord=bX[mn],near.Y.coord=bY[mn])
xy.bound
}
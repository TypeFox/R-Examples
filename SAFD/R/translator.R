translator <-
function(X,nl=101,pic=0){
 #X...2-dim data.frame with columns "x" and "alpha" 
 #containing the vertexes of the polygonal fuzzy number
 #nl...number of levels that shall be calculated, uniformly in [0,1], at least 2
 ok<-checking2(X)
 if(nl<=1){
  print("Minimum number of output levels nl is 2")
  }
 if(nl>1&ok==1){
  levels<-seq(0,1,length=nl)
  A<-subset(X,X$alpha==1)
  cut1<-min(as.numeric(row.names(A)))
  cut2<-max(as.numeric(row.names(A)))
  Left<-X[1:cut1,]
  Right<-X[cut2:nrow(X),]
  L1<-approx(Left$alpha,Left$x, method="linear", rule = 2,n=nl)
  L2<-cbind(L1$y,L1$x)
  L3<-L2
  L4<-data.frame(x=L3[,1],alpha=L3[,2])

  R1<-approx(Right$alpha,Right$x, method="linear", rule = 2,n=nl)
  R2<-cbind(R1$y,R1$x)
  R3<-R2[nrow(R2):1,]
  R4<-data.frame(x=R3[,1],alpha=R3[,2])
  E<-rbind(L4,R4)
  if(pic==1){
  limx<-c(min(c(X$x,E$x))-0.25,max(c(X$x,E$x))+0.25)
 plot(X,type="l",xlim=limx, xlab=NA, ylab=expression(alpha), lwd=2, cex.main=1)
  titletxt <- substitute("Fuzzy number and their " * nl * "-Translator (in red)" , list(nl = as.character(nl)))
      title(main=titletxt,cex.main=1)
 lines(E,type="b",lwd=3, col="red", pch = 21)
   }
  invisible(E)
 }
}

ABWilcox <-
function(behavior,phaseX, v1,v2){
  t1<-table(phaseX)
  tmaxA<-t1[names(t1)==v1]
  startA<-match(v1,phaseX)
  endA<-tmaxA+startA-1
  A<-behavior[startA:endA]
  meanA=mean(A,na.rm=T)
  # x1=rep(1:tmaxA, each = v1)
  
  
  tmaxB<-t1[names(t1)==v2]
  startB<-match(v2,phaseX)
  endB<-tmaxB+startB-1
  
  B=(behavior[startB:endB])
  y=c(A,B)
  xa=rep(1:1,length(A))
  xb=rep(2:2,length(B))
  x=c(xa,xb)       
  meanB=mean(B,na.rm=T)
  u<-wilcox.test(y~x,correct=FALSE,exact=F)
  print(y)
  print(x)
  Means<-c(meanA,meanB) 
  
  print(u)
  
}

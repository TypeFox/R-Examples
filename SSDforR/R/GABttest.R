GABttest <-
function(behavior,phaseX,timeX,v1,v2){
  t1<-table(phaseX)
  tmaxA<-t1[names(t1)==v1]
  startA<-match(v1,phaseX)
  endA<-tmaxA+startA-1
  A1<-behavior[startA:endA]
  weekA<-timeX[startA:endA]
  A<-tapply(A1,weekA,mean,na.rm=T)
  meanA=mean(A,na.rm=T)
  # x1=rep(1:tmaxA, each = v1)
  
  
  tmaxB<-t1[names(t1)==v2]
  startB<-match(v2,phaseX)
  endB<-tmaxB+startB-1
  
  B1=(behavior[startB:endB])
 
  weekB<-timeX[startB:endB]
  B<-tapply(B1,weekB,mean,na.rm=T)
  meanB=mean(B,na.rm=T)
  
  t1<-t.test(A,B)
  var1<-var.test(A,B)  
  t2<-t.test(A,B,var.equal=TRUE)  
  print(t2)
  print(var1)
  print(t1)
  Means<-c(meanA,meanB) 
  
  
  graphics.off()
  barplot(Means,ylab="mean",names.arg=c(v1,v2))
}

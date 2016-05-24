GABrf2 <-
function(behavior,phaseX,timeX,v1){
  maxy=which.max(behavior)
  max<-behavior[maxy]+2
  miny=which.min(behavior)
  min<-behavior[miny]-1
  
  t1<-table(phaseX)
  tmaxA<-t1[names(t1)==v1]
  startA<-match(v1,phaseX)
  
  endA<-tmaxA+startA-1
  A1<-behavior[startA:endA]
  weekA<-timeX[startA:endA]
  
  A<-tapply(A1,weekA,mean,na.rm=T)
  
  firstA=A[1]
  bmean=mean(A)
  endA=length(A)
  
  p1=(.5*(A[1]-bmean)^2)+(.5*(A[endA]-bmean)^2)
  diff=seq(1:endA)
  diff1=seq(1:endA)
  diffsq=seq(1:endA)
  j=1
  while (j <= endA){
    
    diff[j]<-(A[j]-bmean)
    diffsq[j]<-((A[j]-bmean))^2
    j=j+1
  }
  
  for (i in 1:endA) {
    
    j=i-1
    while (j <= endA-1) {
      
      diff1[j]<-diff[j]*diff[i]
      
      j<-j+1
      
    }
  }
  diff1<-diff1[1:endA-1]
  
  p2<-sum(diff1)
  p3<-sum(diffsq)
  p4<-1+(5/(endA-1))
  rf2=((p1+p2)/p3)*p4
  
  tf2.1=((p1+p2)/p3)^2
  tf2.1=1-tf2.1
  tf2.2=tf2.1/(endA+1)
  tf2.3=tf2.2*(p4^2)
  tf2sd=sqrt(tf2.3)
  tf2=rf2/tf2sd
  x1=(c(seq(1:endA)))
  regA<-lm(A~x1)
  l1=c("tf2=",round(tf2,3))
  print(l1)
  
  l2=c("rf2=",round(rf2,3))
  print(l2)
  
  dfx=endA+5
  sig1<-pt(abs(tf2),df=dfx,lower.tail=FALSE)*2
  l3=c("sig of rf2=",round(sig1,3))
  print(l3)
  
  writeLines("----------regression------------")
  print(summary(regA))
  yA<-regA$coefficients[1]
  BetaA<-regA$coefficient[2]
  graphics.off()
  
  plot(x1,A,lwd=2,type="p",col="red", xlab="time", ylab="behavior", main=v1 )
  abline(c(yA,BetaA),col='Blue',lty="dashed")
  
}

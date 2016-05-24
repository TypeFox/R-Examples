diffchart <-
function(behavior,phaseX,v1){
  
  maxy=which.max(behavior)
  max<-behavior[maxy]+2
  miny=which.min(behavior)
  min<-behavior[miny]-1
  
  t1<-table(phaseX)
  tmaxA<-t1[names(t1)==v1]
  startA<-match(v1,phaseX)
  endA<-tmaxA+startA-1
  A<-behavior[startA:endA]
  bmean=mean(A)
  endA=length(A)
  j=2
  diff=seq(1:endA-1)
  for (i in 1:endA) {
    
    j=i-1
    while (j <= endA-1) {
      
      diff[j]<-A[i]-A[j]
      
      j<-j+1
      
    }
  }
  diff<-diff[1:endA-1] 
  maxD=which.max(diff)
  minD=which.min(diff)
  maxdiff=diff[maxD]
  mindiff=diff[minD]
  mindiff=mindiff-2
  maxdiff=maxdiff+2
 

  graphics.off()
  
  layout(rbind(1,2), heights=c(6,1))
  
  plot(diff,ylim=c(mindiff,max),type="l",col="red",bty='L',lty=2,main="Difference Chart")
  lines(A,ylim=c(min,max),type="l",col="blue")
  par(mar=c(1, 1, 1, 1))
  plot.new()
  legend("center", c("diff","Behavior"), col = c("red", "blue"), lty=c(2,1),lwd = 1,ncol=2,bty ="n")
        
  diff=c(diff,NA)
  behavior=c(A,NA)
  
  phase <- rep(v1, endA-1)
  phase <-c(phase,NA)
  transdat<-data.frame(diff,phase)
  a<-readline("Save results? (y/n) ")

  if (a=="y")
  {write.csv(transdat,file = tclvalue(tcl("tk_getSaveFile")),row.names=FALSE)}
  
  print(transdat)
  
}

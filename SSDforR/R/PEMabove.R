PEMabove <-
function(behavior,phaseX,v1,v2){
  
  t1<-table(phaseX)
  tmaxA<-t1[names(t1)==v1]
  startA<-match(v1,phaseX)
  endA<-tmaxA+startA-1
  A<-behavior[startA:endA]
  
  medianA=median(A,na.rm=T)
  
  
  
  tmaxB<-t1[names(t1)==v2]
  startB<-match(v2,phaseX)
  endB<-tmaxB+startB-1
  #tsxB<-behavior[startB:endB]
  B=(behavior[startB:endB])
  
  
  cdcl<-c(A,NA,B)
  y<-na.omit(cdcl)
  total=length(y)
  iv=(1:total)
  end<-which(is.na(cdcl))
  iv<-insert(iv,NA,end)
  
  
  
  #***********************below baseline
  nbelowline<-B>medianA
  p=(sum(nbelowline))/length(B)
  
  
  
  
  maxy=which.max(cdcl)
  
  max<-cdcl[maxy]+1
  numx<-sum(!is.na(cdcl))+3
  par(mfrow=c(3,3)) 
  
  maxy=which.max(behavior)
  
  max<-behavior[maxy]+1
  
  numx<-sum(!is.na(behavior))+3
  
  PND<-c("PEM Above = ",round(p,2))
  print(PND)
  writeLines("-------------------------------------------")
  writeLines(".90 or above = very effective" )
  writeLines(".70 to .89 = moderate effectiveness")
  writeLines(".50 to .69 = debatable effectiveness")
  writeLines(" below .50 = not effective")
  graphics.off()
  layout(rbind(1,2), heights=c(6,1))
  
  plot(iv,cdcl, ylim=c(0,max),lwd=2,type="o",col="red",bty="l", xlab="time", ylab="behavior", main="PEM" )
  
  abline(h=medianA,col="gray",lwd=3)
  
  ab<-NULL
  
  ab<<-recordPlot()
}

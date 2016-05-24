Gindex <-
function(behavior,phaseX,v1,v2){
  
  t1<-table(phaseX)
  tmaxA<-t1[names(t1)==v1]
  startA<-match(v1,phaseX)
  endA<-tmaxA+startA-1
  A<-behavior[startA:endA]
  
  meanA=mean(A,na.rm=T)
  medianA=median(A,na.rm=T)
  x1=(c(seq(1:tmaxA)))
  
  regA<-lm(A~x1)
  rA<-residuals(regA)
  yA<-regA$coefficients[1]
  BetaA<-regA$coefficient[2]
  
  
  
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
  x2=iv[end+1:total]
  regc<-lm(cdcl~iv)
  x1=iv[1:end-1]
  
  #************baseline regress line
  x1<-na.omit(x1)
  Ayhat<-yA+x1*BetaA
  Ayhat<-na.omit(Ayhat)
  #yhatA<-Byhat[startA:endA]
  #startB<-startB-1
  #endB<-endB-1
  #yhatB<-Byhat[startB:endB]
  
  #***********************below baseline
  nbelowline<-A<Ayhat
  nbelowmean<-A<meanA
  nbelowmedian<-A<medianA
  Abrline=(sum(nbelowline))/length(A)
  Abmean=(sum(nbelowmean))/length(A)
  Abmedian=(sum(nbelowmedian))/length(A)
  #*********************above baseline
  naboveline<-A>Ayhat
  nabovemean<-A>meanA
  nabovemedian<-A>medianA
  Aarline=(sum(naboveline))/length(A)
  
  Aamean=(sum(nabovemean))/length(A)
  Aamedian=(sum(nabovemedian))/length(A)
 
  
  #intervention line
  x2<-na.omit(x2)
  Byhat<-yA+x2*BetaA
  Byhat<-na.omit(Byhat)
  yhatA<-Byhat[startA:endA]
  startB<-startB-1
  endB<-endB-1
  yhatB<-Byhat[startB:endB]
  maxy=which.max(cdcl)
  
  max<-cdcl[maxy]+1
  numx<-sum(!is.na(cdcl))+3
  par(mfrow=c(3,3)) 
  
  maxy=which.max(behavior)
  
  max<-behavior[maxy]+1
  
  numx<-sum(!is.na(behavior))+3
  
  
  #***********************below intervention
  nbelowline<-B<Byhat
  nbelowmean<-B<meanA
  nbelowmedian<-B<medianA
  
  Bbrline=(sum(nbelowline))/length(B)
  Bbmean=(sum(nbelowmean))/length(B)
  Bbmedian=(sum(nbelowmedian))/length(B)
  #*********************above intervention
  naboveline<-B>Byhat
  nabovemean<-B>meanA
  nabovemedian<-B>medianA
  
  Barline=(sum(naboveline))/length(B)
 
  Bamean=(sum(nabovemean))/length(B)
  Bamedian=(sum(nabovemedian))/length(B)
  #*************************************
  gamean=Bamean-Aamean
  galine=Barline-Aarline
  
  gamedian=Bamedian-Aamedian
  
  gbmean=Bbmean-Abmean
  gbline=Bbrline-Abrline
  gbmedian=Bbmedian-Abmedian
  writeLines("")
  l1<-c("small effect size: <.3")
  l2<-c("medium effect size: .31 to .5 ")
  l3<-c("large effect size: >.51")
  writeLines(l1)
  writeLines(l2)
  writeLines(l3)
  writeLines("")
  writeLines("------------------------------------------------------g-index-------------------------------------------------------------")
  writeLines(" ")
  writeLines("-----------------------------------------------------Above Lines----------------------------------------------------------")
  l1<-c("G mean= ",round(gamean,3),"  G median= ",round(gamedian,3) ,"G Regression line= ",round(galine,3))
  print(l1)
  writeLines(" ")
  writeLines("-----------------------------------------------------Below Lines-----------------------------------------------------------")
  l2<-c("G mean= ",round(gbmean,3),"  G median= ", round(gbmedian,3) ,"G Regression line= ",round(gbline,3))
  print(l2)
  
  graphics.off()
  layout(rbind(1,2), heights=c(6,1))
  
  plot(iv,cdcl, ylim=c(0,max),lwd=2,type="o",col="red",bty="l", xlab="time", ylab="behavior", main="g-index" )
  
  abline(reg=regA,col='blue',lty="dotted",lwd=3)
  abline(h=meanA,col="gray",lwd=2)

  abline(h=medianA,col='black',lwd=2)
  par(mar=c(1, 1, 1, 1))
  plot.new()
  legend("center", c("median","mean","regression line"), col = c("black","gray","blue"), lty=c("solid","solid","dotted"),lwd = 3,ncol=3,bty ="n")
}

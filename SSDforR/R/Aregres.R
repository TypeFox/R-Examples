Aregres <-
function(behavior,phaseX,v1){
  t1<-table(phaseX)
  tmaxA<-t1[names(t1)==v1]
  startA<-match(v1,phaseX)
  endA<-tmaxA+startA-1
  A<-behavior[startA:endA]
  
  
  x1=(c(seq(1:tmaxA)))
  
  regA<-lm(A~x1)
  rA<-residuals(regA)
  yA<-regA$coefficients[1]
  BetaA<-regA$coefficient[2]
  

  graphics.off()
  plot(x1, A,lwd=2,type="o",col="red", xlab="time", ylab="behavior", bty='L',main=v1 )
  abline(c(yA,BetaA),col='Blue',lty="dashed")
  print(summary(regA))
  
  
}

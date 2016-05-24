ABanova <-
function(behavior,phaseX){
  meanA<-tapply(behavior,phaseX,mean,rm.na=T)
  aovx<-aov(behavior~as.factor(phaseX))
  print(summary(aovx))
  print(meanA)
  tukeyx<-TukeyHSD(aovx)
  graphics.off()
  plot(tukeyx)
  tukeyx
  
}

trimabove <-
function (behavior,phaseX,v1,v2) {
  medianab<-tapply(behavior,phaseX,mean,trim=.1)
  medianx=medianab[names(medianab)==v1]
  dzone<-behavior > medianx
  tm<-table(dzone,phaseX) 
  ctbl<-cbind(tm[,v1],tm[,v2])
  print(ctbl)
  print(prop.table(ctbl,1)*100)
  print(prop.table(ctbl,2)*100)
  c1<-chisq.test(ctbl,correct=FALSE)
  f1<-fisher.test(ctbl,alternative = "two.sided")
  print(c1)
  print(f1)
  graphics.off()
  layout(rbind(1,2), heights=c(6,1))
  plot(behavior,col="red",type="o")
  abline(h=medianx,col="green")
  par(mar=c(1, 1, 1, 1))
  plot.new()
  legend("center", c("behavior","trimmed mean"), col = c("red","green"), lwd = 1,ncol=2,bty ="n")
}

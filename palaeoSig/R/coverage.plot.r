coverage.plot<-function(mod,fos, rare=5, identify=FALSE){ 
  modfos<-Merge(mod, fos, split=TRUE)
  modmax<-sapply(modfos$mod,max)
  fosmax<-sapply(modfos$fos,max)
  n2<-Hill.N2(modfos$mod)
  n2[is.infinite(n2)]<-0
  cols<-rep(1,length(modmax))
  cols[n2<rare]<-4
  cols[n2==0]<-2

  plot(modmax,fosmax, xlab="Calibration set maximum abundance %", ylab="Fossil data maximum abundance %", col=cols, pch=ifelse(n2<=rare, 16,1))
  abline(0,1, col=2) #add a 1:1 line.
  if(identify)identify(modmax,fosmax, labels=names(modmax), cex=0.7)
  return(invisible(data.frame(modmax=modmax, fosmax=fosmax, n2=n2)))
}

#data(RLGH)
#coverage.plot(mod=SWAP$spec, fos=RLGH$spec, identify=FALSE)
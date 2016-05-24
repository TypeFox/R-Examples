plot.RNE<-function (x, which=1, ylim,...){
  z<-list()
  z$r<-x$random[,c(1,1+which)]
  z$n<-t(sapply(x$neigh,function(b)c(b$neigh,b$effn,b$hb.r2[which])))
  z$e<-sapply(x$neighbour,function(b)b$eb.r2[which]) 
  if(missing(ylim)) ylim<-range(c(z$r[,2],z$n[, 3]))
  plot(1-z$r[, 1],z$r[, 2],ylim=ylim,type="b",xlab="Fraction of sites deleted",ylab = expression("r"^2), ...)
  matpoints(z$n[, 2], cbind(z$n[,3],z$e), type = "b", pch = c(16, 3), lty = c(1, 2))
  text(z$n[, 2:3], label = paste(z$n[, 1], "km"), adj = -0.1)
  z
}


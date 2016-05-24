PLOT.cg <-
function(x){
    if(x$data[1,2]<=x$data[length(x$data[,2][x$data[,3]==unique(x$data[,3])[1]]),2])
       { m=1
    }else { print("cracks are not growing")}
  plot.data=function(y=x$data,z=x$a.F,u=x$Tc, ...){
      COL=c(1:7,colors()[82:150])
      plot.data=plot(subset(y[,c(1,2)],y[,3]==unique(y[,3])[1]),type="b",
                     xlim=c(min(y[,1]),max(y[,1])),ylim=c(min(y[,2]),max(y[,2])),
                     col=1,las=1,pch=20,cex=1.5,...)
      for(i in 2:length(unique(y[,3])))
      points(subset(y[,c(1,2)],y[,3]==unique(y[,3])[i]),type="b",col=COL[i],
      pch=20,cex=1.5)
      abline(h=z, col="gray30",lwd=3,lty=2)
      abline(v=u, col="gray30",lwd=3,lty=2)
  }
  plot.pred=function(y=x$data,z=x$a.F,u=x$Tc,v=x$crack.pred,w=x$F.emp, ...){
      COL=c(1:7,colors()[82:150])
      plot.pred=plot(subset(y[,c(1,2)],y[,3]==unique(y[,3])[1]),type="p",
                     xlim=c(min(y[,1]),max(w[,1])),ylim=c(min(y[,2]),max(y[,2])),
                     col=1,las=1,pch=20,cex=1.5,...)
      abline(h=z, col="gray30",lwd=3,lty=2)
      abline(v=u, col="gray30",lwd=3,lty=2)
      for(i in 1:length(unique(v[,3])))
      lines(subset(v[,c(1,2)],v[,3]==unique(v[,3])[i]), col=COL[i],lwd=2)
      for(i in 1:length(unique(y[,3])))
      points(subset(y[,c(1,2)],y[,3]==unique(y[,3])[i]),type="p",pch=20,cex=1.5)
      points(w[,1], rep(z,length(w[,1])), cex=2, col=2, lwd=2)
  }
  
  plot.F=function(y=x$F.est, u=x$Tc, w=x$F.emp, ...){
      plot.F=plot(y, xlim=c(min(y[,1]),max(w[,1])), col=1,las=1,pch=20,
                  cex=1.5,...)
      abline(v=u, col="gray30",lwd=3,lty=2)
      points(w, col=4, pch=20, cex=2)
  }
  plot.resid=function(y=x$crack.est, u=x$residuals, ...){
      plot.resid=plot(y[,2],u ,pch=20, ...)
      abline(h=0, col="gray30", lwd=2)
  }
  list(plot.data=plot.data, plot.pred=plot.pred, plot.F=plot.F,
       plot.resid=plot.resid)
}

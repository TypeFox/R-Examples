draw.rearranged.oriloc <- function(rearr.ori,breaks.gcfw=NA,breaks.gcrev=NA,breaks.atfw=NA,breaks.atrev=NA){

  rearat=cumsum(rearr.ori$atskew.rear)
  reargc=cumsum(rearr.ori$gcskew.rear)
  strand.rear=rearr.ori$strand.rear
  cds.rear=rep(1,length(strand.rear))
  cds.rear[strand.rear=="reverse"]=-1
  rearcds=cumsum(cds.rear)
  meancoord.rear=rearr.ori$meancoord.rear

  ymin <- min(rearat,reargc)
  ymax <- max(rearat,reargc)
  xmin <- min(meancoord.rear)
  xmax <- max(meancoord.rear)
  
  ticksrear <- pretty(rearcds)
  
  ticks.yrear <- (ymax-ymin)/(max(rearcds)-min(rearcds))*(ticksrear - min(rearcds)) + ymin
  cds.yrear   <- (ymax-ymin)/(max(rearcds)-min(rearcds))*(rearcds - min(rearcds)) + ymin

 plot(meancoord.rear,reargc, type="l", xlab="Map position (gene index)",
      ylab = "Cumulated normalized skew",xlim=c(xmin,xmax),ylim=c(ymin,ymax),cex.lab=1.35,col="blue",main="Rearranged nucleotide skews",lwd=2,cex.main=1.4)

  abline(v=sum(strand.rear=="forward"),col="black",lwd=2)
  
  axis(side = 4, at = ticks.yrear, labels = ticksrear, col = "black", col.axis ="black")

  tmp <- pretty(meancoord.rear)
  abline(v=tmp, col="grey", lty=3,lwd=1.5)
  tmp <- tmp[-length(tmp)] + diff(tmp)/2
  abline(v=tmp, col="grey", lty=3,lwd=1.5)
  
  lines(meancoord.rear,rearat, col="red",lwd=2)
  
  lines(meancoord.rear,cds.yrear, col="black",lwd=2)
  
  mtext("Cumul. A-T skew", col="red", adj=0,cex=1.1)
  mtext("Cumul. G-C skew", col="blue",cex=1.1)
  mtext("Cumul. CDS skew", col="black", adj=1,cex=1.1)

  if(sum(is.na(breaks.gcfw))==0){
    segments(x0=meancoord.rear[breaks.gcfw],y0=reargc[breaks.gcfw]-(ymax-ymin)/30,x1=meancoord.rear[breaks.gcfw],y1=reargc[breaks.gcfw]+(ymax-ymin)/30,col="blue",lwd=2,lty=1)
  }
  
  if(sum(is.na(breaks.gcrev))==0){
    segments(x0=meancoord.rear[breaks.gcrev],y0=reargc[breaks.gcrev]-(ymax-ymin)/30,x1=meancoord.rear[breaks.gcrev],y1=reargc[breaks.gcrev]+(ymax-ymin)/30,col="blue",lwd=2,lty=1)
  }
  
  if(sum(is.na(breaks.atfw))==0){
    segments(x0=meancoord.rear[breaks.atfw],y0=rearat[breaks.atfw]-(ymax-ymin)/30,x1=meancoord.rear[breaks.atfw],y1=rearat[breaks.atfw]+(ymax-ymin)/30,col="red",lwd=2,lty=1)
  }

  if(sum(is.na(breaks.atrev))==0){
    segments(x0=meancoord.rear[breaks.atrev],y0=rearat[breaks.atrev]-(ymax-ymin)/30,x1=meancoord.rear[breaks.atrev],y1=rearat[breaks.atfw]+(ymax-ymin)/30,col="red",lwd=2,lty=1)
  }

}

# Version: 30-11-2012, Daniel Fischer

estPlotSingle <- function(x,col,highlight,hlCol,pch,zoom,...){

  ifelse(is.matrix(x$probs),Nprobs <- dim(x$probs)[1], Nprobs <- length(x$probs))
  ifelse(is.matrix(x$probs),Nval <- dim(x$probs)[2], Nval <- 1)
  Nrows <- 1

  markThese <- rep(col,Nval)
  markThese[highlight] <- hlCol
  comb <- x$goi

 par(mfrow=c(Nprobs-1,Nprobs-1),
      pty="s",
      bty="c",
      oma=c(0,0,3,0),
      mar=c(0,0,0,0))

  for(xRun in 1:(Nprobs-1))
  { 
    yValues <- getSingleProb(x,comb[xRun])
    yMin <- ifelse(zoom,min(yValues),0)
    yMax <- ifelse(zoom,max(yValues),1)

    for(yRun in (xRun+1):Nprobs)
    {
      if(xRun==1 & yRun == 2) par(oma=c(0,3,3,0), mar=c(0,0,0,0))
      #if(xRun==(Nprobs-1) & yRun == Nprobs) par(bty="o")
      ifelse(yRun==(Nprobs), par(bty="o"),par(bty="c"))

      xValues <- getSingleProb(x,comb[yRun])
      xMin <- ifelse(zoom,min(xValues),0)
      xMax <- ifelse(zoom,max(xValues),1)

      plot(xValues,yValues,xlim=c(max(xMin-0.05,0),min(xMax+0.05,1)),ylim=c(max(yMin-0.05,0),min(yMax+0.05,1)),xaxt="n",yaxt="n",col=markThese,pch=pch,xlab="",ylab="",xaxt="n",yaxt="n")
      lines(c(0.5,0.5),c(-1,2),lty="dotted")
      lines(c(-1,2),c(0.5,0.5),lty="dotted")
      
      if(xRun==1) mtext(paste("P(",comb[yRun],")",sep=""), NORTH<-3, line=0.3, adj=0.5, cex=0.8, col="black", outer=FALSE)
      if(xRun==1 & yRun == 2) mtext(paste("P(",comb[xRun],")",sep=""), WEST<-2, line=0.3, adj=0.5, cex=0.8, col="black", outer=FALSE)
      if(xRun==1 & yRun == 2) axis(2,at=seq(0,1,0.2),labels=seq(0,1,0.2),tck=0.01,mgp=c(0,-1.5,0))
      axis(1,at=seq(0,1,0.2),labels=seq(0,1,0.2))
    }
    if(xRun!=(Nprobs-1))
    { 
      for(i in 1:xRun)
      {
	plot(-10,-10,xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n",pch=pch,xlab="",ylab="",xaxt="n",yaxt="n")
        if(i==1)
	{
	  mtext(paste("P(",comb[xRun+1],")",sep=""), WEST<-2, line=0.3, adj=0.5, cex=0.8, col="black", outer=FALSE)
	  axis(2,at=seq(0,1,0.2),labels=seq(0,1,0.2),tck=0.01,mgp=c(0,-1.5,0))
	}
      }
    }
  }
}
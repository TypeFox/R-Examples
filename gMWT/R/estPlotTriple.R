# Version: 30-11-2012, Daniel Fischer

estPlotTriple <- function(x,col="black",highlight=NULL,hlCol="red",pch=20,zoom=F){

  if(length(x$probs)==1) stop("You only have 1 probabilistic index, I cannot plot anything, sorry!!!\n")

  ifelse(is.matrix(x$probs),Nprobs <- dim(x$probs)[1], Nprobs <- length(x$probs))
  ifelse(is.matrix(x$probs),Nval <- dim(x$probs)[2], Nval <- 1)

  markThese <- rep(col,Nval)
  markThese[highlight] <- hlCol

  comb <- getComb(x$goi,type="pairs",order=T)
  tripleComb <- getTripleComb(x$goi,x$order)

  Ncols <- Nprobs*(Nprobs-1)/2

  par(mfrow=c(Nprobs-1,Nprobs-1),
      pty="s",
      bty="c",
      oma=c(0,0,3,0),
      mar=c(0,0,0,0))

  run <- 1
  for(xRun in 1:(Nprobs-1))
  { 
    # First fill up with empty fields
    if(xRun>1)
    {
      for(i in 1:(xRun-1))
      {
	plot(-10,-10,xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n",pch=pch,xlab="",ylab="",xaxt="n",yaxt="n")
	if(i==1)
	{
	    name1 <- paste("P(", tripleComb[run,1],"<",tripleComb[run,2],"<",tripleComb[run,3],")",sep="")
	    if(is.matrix(x$probs)){
	      mtext(paste(rownames(x$probs)[rownames(x$probs)==name1],sep=""), WEST<-2, line=0.3, adj=0.5, cex=0.8, col="black", outer=FALSE)
	    } else {
	      mtext(paste(names(x$probs)[names(x$probs)==name1],sep=""), WEST<-2, line=0.3, adj=0.5, cex=0.8, col="black", outer=FALSE)
	    }
	    axis(2,at=seq(0,1,0.2),labels=seq(0,1,0.2),tck=0.01,mgp=c(0,-1.5,0))
	}
      }
    }

    yValues <- getTripleProb(x,tripleComb[run,1],tripleComb[run,2],tripleComb[run,3])
    yMin <- ifelse(zoom,min(yValues),0)
    yMax <- ifelse(zoom,max(yValues),1)

    # Then fill the right pictures there!
    for(yRun in (xRun+1):Nprobs)
    {
      if(xRun==1 & yRun == 2) par(oma=c(0,3,3,0), mar=c(0,0,0,0))
#       if(xRun==(Nprobs-1) & yRun == Nprobs) par(bty="o")
      ifelse(yRun==(Nprobs), par(bty="o"),par(bty="c"))

      xValues <- getTripleProb(x,tripleComb[run,4],tripleComb[run,5],tripleComb[run,6])
      xMin <- ifelse(zoom,min(xValues),0)
      xMax <- ifelse(zoom,max(xValues),1)

      plot(xValues,yValues,xlim=c(max(xMin-0.05,0),min(xMax+0.05,1)),ylim=c(max(yMin-0.05,0),min(yMax+0.05,1)),xaxt="n",yaxt="n",col=markThese,pch=pch,xlab="",ylab="",xaxt="n",yaxt="n")
      lines(c(1/6,1/6),c(-1,2),lty="dotted")
      lines(c(-1,2),c(1/6,1/6),lty="dotted")
      if(!is.null(highlight)) {
	xChord <- xValues[highlight]
	yChord <- yValues[highlight]
	points(xChord,yChord,col=hlCol,pch=pch)
      }      

      name1 <- paste("P(", tripleComb[run,1],"<",tripleComb[run,2],"<",tripleComb[run,3],")",sep="")
      name2 <- paste("P(", tripleComb[run,4],"<",tripleComb[run,5],"<",tripleComb[run,6],")",sep="")

      if(xRun==1){
	if(is.matrix(x$probs)){
	   mtext(paste(rownames(x$probs)[rownames(x$probs)==name2],sep=""), NORTH<-3, line=0.3, adj=0.5, cex=0.8, col="black", outer=FALSE)
	} else {
	   mtext(paste(names(x$probs)[names(x$probs)==name2],sep=""), NORTH<-3, line=0.3, adj=0.5, cex=0.8, col="black", outer=FALSE)
	}
      }

      if(yRun==2){
	if(is.matrix(x$probs)){
	  mtext(paste(rownames(x$probs)[rownames(x$probs)==name1],sep=""), WEST<-2, line=0.3, adj=0.5, cex=0.8, col="black", outer=FALSE)
	} else {
	  mtext(paste(names(x$probs)[names(x$probs)==name1],sep=""), WEST<-2, line=0.3, adj=0.5, cex=0.8, col="black", outer=FALSE)
	}
      }
      if(xRun==1 & yRun == 2) axis(2,at=seq(0,1,0.2),labels=seq(0,1,0.2),tck=0.01,mgp=c(0,-1.5,0))
      axis(1,at=seq(0,1,0.2),labels=seq(0,1,0.2))
      run <- run + 1
    }
  }
}
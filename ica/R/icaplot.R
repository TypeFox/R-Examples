icaplot <-
  function(xseq=seq(-2,2,length.out=500),xlab="",ylab="",
           lty=1,lwd=1,col="black",...){
    
    if(length(lty)!=18L){ lty <- rep(lty[1],18) }
    if(length(lwd)!=18L){ lwd <- rep(lwd[1],18) }
    if(length(col)!=18L){ col <- rep(col[1],18) }
    xlim <- range(xseq)
    par(mfrow=c(6,3))
    for(i in 1:18){
      myfun <- as.character(letters[i])
      kurto <- icasamp(myfun,"kur")
      myden <- icasamp(myfun,"pdf",data=xseq)
      tit1p <- bquote("("*.(myfun)*")")
      tit2p <- bquote(k==.(round(kurto,2)))
      mytit <- bquote(.(tit1p)*"  "*.(tit2p))
      plot(xseq,myden,type="l",ylim=c(0,max(myden)+.1),
           xlab=xlab,ylab=ylab,main=mytit,lty=lty[i],
           lwd=lwd[i],col=col[i],...)
    }
    
  }
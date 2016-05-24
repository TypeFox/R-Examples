
plot.VPdtw <- function(x,type=c("All","Before","After","Shift"),xlim=NULL,...)
  {
    bgcol <- grey(0.9)
    type <- match.arg(type,c("All","Before","After","Shift"))

    pp <- switch(type,
                 All=par(mfrow=c(3,1)),
                 Before=par(mfrow=c(1,1)),
                 After=par(mfrow=c(1,1)),
                 Shift=par(mfrow=c(1,1)))
    
    if(is.null(xlim)) {
      xlim <- c(x$xVals[1],x$xVals[length(x$xVals)])
      ylim <- range(x$query,x$reference,na.rm=TRUE)
    } else {

      ind <- 1:length(x$reference)
      ylim <- range(x$reference[which(ind>=xlim[1] & ind <= xlim[2])],na.rm=TRUE)
      ind <- 1:length(x$query)
      wind <- which(ind>=xlim[1]-350 & ind <= xlim[2]+350)
      if(is.matrix(x$query)) {
        ylim <- range(c(ylim,x$query[wind,]),na.rm=TRUE)
      } else {
        ylim <- range(c(ylim,x$query[wind]),na.rm=TRUE)
      }
      
    }
    
    if(type=="All" | type=="Before") {
      
      if(is.matrix(x$query)) {
        plot(c(1,nrow(x$query)),c(1,ncol(x$query)),type="n",
             xlab="Index",ylab="Sample",
             main="Queries before Alignment",xlim=xlim)
        rect(-1000,-1000,nrow(x$query)+1000,ncol(x$query)+1000,col=bgcol)
        image(1:nrow(x$query),1:ncol(x$query),x$query,add=TRUE)
        box()
      }
      
      if(is.vector(x$query)) {
        
        plot(xlim,ylim,type="n",xlab="Index",ylab="Intensity",
             main="Query and Reference before Alignment")
        lines(x$xVals,x$reference,lwd=2,col=1)
        lines(x$query,col=2)
        
      }
    }

    if(type=="All" | type=="After") {
      if(is.matrix(x$query) & is.vector(x$penalty)) {
        plot(c(min(x$xVals),max(x$xVals)),c(1,ncol(x$warpedQuery)),type="n",
             xlab="Index",ylab="Sample",
             main="Queries after Alignment")
        rect(min(x$xVals)-1000,-1000,max(x$xVals)+1000,ncol(x$warpedQuery)+1000,col=bgcol)
        image(x$xVals,1:ncol(x$warpedQuery),x$warpedQuery,add=TRUE)
        box()
      }
      if(is.vector(x$query) & is.vector(x$penalty)) {

        plot(xlim,ylim,type="n",xlab="Index",ylab="Intensity",main="Query and Reference after Alignment")
        lines(x$xVals,x$reference,lwd=2,col=1)
        lines(x$xVals,x$warpedQuery,col=2)
        
      }
      if(is.matrix(x$query) & is.matrix(x$penalty)) {
        ## should happen
        return("Error, exiting...\n")
      }
      if(is.vector(x$query) & is.matrix(x$penalty)) {
        ##
        ncols <- ncol(x$warpedQuery)
        plot(xlim,ylim,type="n",xlab="Index",ylab="Intensity",main="Query and Reference after Alignment")
        lines(x$xVals,x$reference,lwd=2,col=1)
        matplot(x$xVals,x$warpedQuery,type="l",lty=1,add=TRUE,col=2:(ncols+1))
        legend("topright",legend=c("Reference",paste("penalty  #",1:ncols,sep="")),col=1:(ncols+1),lty=rep(1,ncols+1),lwd=c(2,rep(1,ncols)))
      }
    }

    if(type=="All" | type=="Shift") {
      if(is.matrix(x$shift)) {
        ## Image isn't very useful here, matplot is better, though crowded
        ##image(x$xVals,1:ncol(x$shift),x$shift,xlab="Index",ylab="Sample",main="Shifts after Alignment")
        ##box()
        ncols <- ncol(x$shift)
        ##xlim <- range(0,x$shift,na.rm=TRUE)
        
        matplot(x$xVals,x$shift,type="l",lty=1,xlab="Index",ylab="Shift",main="Shifts required for Alignment",col=2:(ncols+1),xlim=xlim)
        if(is.matrix(x$penalty)) legend("topleft",legend=paste("penalty  #",1:ncols,sep=""),col=2:(ncols+1),lty=rep(1,ncols))
        if(is.matrix(x$query)) legend("topleft",legend=paste("query  #",1:ncols,sep=""),col=2:(ncols+1),lty=rep(1,ncols))
        abline(h=0,lty=2,col=grey(0.75))
      }
      if(is.vector(x$shift)) {
        plot(x$xVals,x$shift,ylim=range(c(0,x$shift),na.rm=TRUE),xlab="Index",ylab="Shift",main="Shifts required for Alignment",type="n",xlim=xlim)
        abline(h=0,lty=2,col=grey(0.75))
        lines(x$xVals,x$shift)
               
        
      }
    }

    par(pp)
  }

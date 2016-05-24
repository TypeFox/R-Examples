`plot.freqMAP` <-
function(x,y=NULL,
         xlim=NULL,ylim=NULL,legend=NULL,
         show.p.value.legend=FALSE,type="freq",
         p.value.bar.alpha=c(.05,.01),
         p.value.bar.color=c("gray90","darkgray"),
         cex=1,cex.axis=1,cex.lab=1,cex.main=1,
         pch.x=2,lty.x=1,lwd.x=1,col.x="red",
         pch.y=1,lty.y=2,lwd.y=1,col.y="blue",
         cex.legend=1,
         layout.matrix=NULL,...){

  twomaps <- TRUE
  if(is.null(y)) twomaps <- FALSE

  if(!twomaps & type != "freq"){
    stop("only allowable type for single MAP plot is \"freq\"")
  }
  
  if(!(type %in% c("freq","or"))){
    stop("type must be either \"freq\" or \"or\"")
  }

  if(!is.null(legend) & (type=="or" | is.null(y))){
    cat("legend argument is ignored when type = \"or\" or y=NULL\n")
    legend=NULL
  }
    
  if(!is.null(legend)){
    if(length(legend)!=2){
      stop("legend must be either NULL or a string vector of length 2")
    }
  }

  #The x values to plot
  x1 <- x$cat.ma[,x$x.label]
  if(twomaps) x2 <- y$cat.ma[,y$x.label]

  if(show.p.value.legend & is.null(y)){
    cat("p-value legend not applicable in single non-comparison plots. Ignoring.\n")
    show.p.value.legend <- FALSE
  }
  

  
  if(is.null(xlim)){
    xlim = range(x1)
    if(twomaps) xlim = range(c(xlim,x2))
    #if either legend is active, then add some space on the left
    if(show.p.value.legend | !is.null(legend))
      xlim <- c(xlim[1]-0.2*(xlim[2]-xlim[1]),xlim[2])
  }else{
    if(length(xlim)!=2 ){
      stop("xlim must be NULL or a numeric vector of length 2")
    }
  }

  if(twomaps){
    if(any(x$cat.names!=y$cat.names)){
      stop("x and y must have the same category names, in the same order")
    }
    
    if(x$x.label!=y$x.label){
      stop("x.label must be identical in x and y")
    }

    if(x$hw != y$hw){
      stop("hw must be same in x and y")
    }

  }

  #or plots are made on a log scale, so make sure ylims are both positive
  if(type=="or"&any(ylim<=0))
    stop("for type=\"or\", ylim must be greater than 0")
  
  cat.names <- x$cat.names
  cat.short <- x$cat.short

  
  #This will be used to plot Bayesian p-values at the bottom of the figure
  if(twomaps) pc <- posterior.comparison.freqMAP(x,y)
  
  if(type=="freq"){ #Frequency plots

    if(!is.null(ylim)){
      if(ncol(ylim)!=2){
        stop("for type=\"freq\", ylim must be a numeric matrix with 2 columns and one row for each category")
      }
      if(nrow(ylim)!=length(cat.names)){
        stop("for type=\"freq\", ylim must be a numeric matrix with 2 columns and one row for each category")
      }
    }

    if(is.null(layout.matrix)) layout.matrix <- matrix(1:length(cat.names),ncol=1)
    
    layout(layout.matrix)
    
    i <- 0
    for(a in cat.names){
      i <- i + 1
      
      if(!is.null(ylim)){
        yl = ylim[i,]
      }else{
        yl <- range(c(range(x$cat.ma[,c(a,paste(a,c(".lpi",".upi"),sep=""))],
                            na.rm=TRUE)))

        if(twomaps)
          yl <- range(c(yl,
                        range(y$cat.ma[,c(a,paste(a,c(".lpi",".upi"),sep=""))],
                              na.rm=TRUE)))
      }
      
      
      plot(x1,x$cat.ma[,a],ylab="Freq.",
           main=paste("\"",cat.short[i],"\"",sep=""),
           xlab=x$x.label,ylim=yl,xlim=xlim,type="n",
           cex=cex,cex.axis=cex.axis,cex.lab=cex.lab,cex.main=cex.main)

      if(twomaps){
        draw.signif.post.boxes(post.dat=pc[,c(x$x.label,paste(a,".gr1.gt.gr2",sep=""))],
                               ylims=yl,alpha1=p.value.bar.alpha[1],alpha2=p.value.bar.alpha[2],
                               color1=p.value.bar.color[1],
			       color2=p.value.bar.color[2],
                               dens1=NULL,dens2=NULL)
      }
      
      points(x1,x$cat.ma[,a],pch=pch.x,col=col.x,cex=cex)
      lines(x1,x$cat.ma[,paste(a,".lpi",sep="")],col=col.x,lty=lty.x,lwd=lwd.x)
      lines(x1,x$cat.ma[,paste(a,".upi",sep="")],col=col.x,lty=lty.x,lwd=lwd.x)
      if(twomaps){
        points(x2,y$cat.ma[,a],col=col.y,cex=cex,pch=pch.y)
        lines(x2,y$cat.ma[,paste(a,".lpi",sep="")],col=col.y,lty=lty.y,lwd=lwd.y)
        lines(x2,y$cat.ma[,paste(a,".upi",sep="")],col=col.y,lty=lty.y,lwd=lwd.y)
      }
      
      abline(h=axTicks(side=2),lty=3)
      if(!is.null(legend)){
        legend(x="topleft",
               legend=legend,pch=c(pch.x,pch.y),
               col=c(col.x,col.y),cex=cex.legend,bty="n")
        
      }
      if(show.p.value.legend)
        add.p.value.legend(alpha1=p.value.bar.alpha[1],
                           alpha2=p.value.bar.alpha[2],
                           col1=p.value.bar.color[1],
			   col2=p.value.bar.color[2],
                           cex=cex.legend)      
    }
    
  }else if(type=="or"){ #odds ratio plots

    if(is.null(layout.matrix))
      layout.matrix <- matrix(1:(length(cat.names)*(length(cat.names)-1)/2),ncol=1)
    layout(layout.matrix)

    if(!is.null(ylim)){
      if(length(ylim)!=2){
        stop("for type=\"or\", ylim must be a vector of length 2")
      }
      yl=ylim
    }else{
      #a common set of y limits on all or plots
      yl <- range(exp(pc[,grep(".lor",names(pc))]),na.rm=TRUE)
    }
    
    
    for(i in 1:(length(cat.names)-1)){
      for(j in (i+1):length(cat.names)){
        
        #column names for stats on lor: mean, lower and upper
        #post. bounds, and Prob(lor>0)
        cn <- paste(cat.names[j],".",cat.names[i],".lor",
                    c(".mean",".lpi",".upi",".p.gt.0"),sep="")
        
        plot(0,1,main=paste("\"",cat.short[j],"\" vs. \"",cat.short[i],"\"",sep=""),
             ylim=yl,xlim=xlim,xlab=x$x.label,ylab="OR",type="n",log="y",
             cex.axis=cex.axis,cex.lab=cex.lab,cex.main=cex.main)
        
        draw.signif.post.boxes(post.dat=pc[,c(x$x.label,cn[4])],
                               ylims=yl,alpha1=p.value.bar.alpha[1],
                               alpha2=p.value.bar.alpha[2],
			       color1=p.value.bar.color[1],
                               color2=p.value.bar.color[2],
			       dens1=NULL,dens2=NULL,
                               ylog=TRUE)
        points(x1,exp(pc[,cn[1]]),pch=pch.x,col=col.x,cex=cex)
        lines(x1,exp(pc[,cn[2]]),lwd=lwd.x,lty=lty.x,col=col.x)
        lines(x1,exp(pc[,cn[3]]),lwd=lwd.x,lty=lty.x,col=col.x)
        abline(h=axTicks(side=2,log=TRUE),lty=3)
        abline(h=1,lwd=2)

        if(show.p.value.legend)
          add.p.value.legend(alpha1=p.value.bar.alpha[1],
			     alpha2=p.value.bar.alpha[2],
			     col1=p.value.bar.color[1],
			     col2=p.value.bar.color[2],
                             cex=cex.legend)

      }
    }

 
  }


  invisible()
}


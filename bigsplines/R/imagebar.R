imagebar <-
  function(x,y,z,xlim=NULL,ylim=NULL,zlim=NULL,
           zlab=NULL,zcex.axis=NULL,zcex.lab=NULL,
           col=NULL,ncolor=21,drawbar=TRUE,zline=2,
           pltimage=c(.2,.8,.2,.8),pltbar=c(.82,.85,.2,.8),...){
    
    ### define colors
    if(is.null(col[1])){
      col <- c("blueviolet","blue","cyan","green","yellow","orange","red")
    }
    col <- colorRampPalette(col)(ncolor)
    
    ### get z limits and label
    if(is.null(zlim)){zlim <- range(c(z))}
    if(is.null(zlab)){zlab <- "z"}
    if(is.null(zcex.lab)){zcex.lab <- 1}
    scales <- ncolor/(zlim[2]-zlim[1])
    breaks <- seq(zlim[1],zlim[2],length.out=(ncolor+1))
    
    # plot image and color bar
    oldplt <- par()$plt
    on.exit(par(plt=oldplt))
    if(drawbar){
      # plot color bar
      par(plt=pltbar)
      plot(1,1,t="n",ylim=zlim,xlim=c(0,1),xaxt="n",yaxt="n",xaxs="i",yaxs="i",xlab="",ylab="")
      axis(4,cex.axis=zcex.axis)
      mtext(zlab,side=4,line=zline,cex=zcex.lab*par()$cex)
      for (ii in 1:ncolor) {
        idx <- zlim[1] + (ii-1)/scales
        rect(0,idx,1,(idx+1/scales),col=col[ii],border=NA)
      }
      # plot image 
      par(plt=pltimage,new=TRUE)
      image(x,y,z,col=col,breaks=breaks,...)
    } else {
      # plot image
      par(plt=pltimage)
      image(x,y,z,col=col,breaks=breaks,...)
    }
    
  }
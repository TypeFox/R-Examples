eegspace <- 
  function(space,voltage,vlim=NULL,mycolors=NULL,ncolor=25,
           colorbar=TRUE,nctick=5,rtick=1,cex.axis=1,barloc=NULL,
           colorlab=NULL,cex.lab=1,plotaxes=FALSE,main="",
           xyzlab=NULL,cex.point=1,cex.main=1,nose=TRUE,ears=TRUE,
           head=TRUE,col.head="AntiqueWhite",mar=NULL,...){
    ###### Plots EEG Spatial Map
    ###### Nathaniel E. Helwig (helwig@umn.edu)
    ###### Last modified: February 16, 2015
    
    ### initial checks
    vlen <- length(voltage)
    space <- as.matrix(space)
    if(nrow(space)!=vlen){stop("Inputs 'space' and 'voltage' must have same number of rows.")}
    if(!any(ncol(space)==c(2L,3L))){stop("Input 'space' must be two- or three-column matrix of coordinates.")}
    if(!is.null(vlim[1])){
      if(vlim[1]>min(voltage)){warning("First entry of input 'vlim' is larger than min(voltage).")}
      if(vlim[2]<max(voltage)){warning("Second entry of input 'vlim' is smaller than max(voltage).")}
    } else {vlim <- range(voltage)}
    ncolor <- as.integer(ncolor)
    if(ncolor<1){stop("Input 'ncolor' must positive integer.")}
    nctick <- as.integer(nctick[1])
    if(nctick<0){stop("Input 'nctick' must nonnegative integer.")}
    rtick <- as.integer(rtick[1])
    if(rtick<1){stop("Input 'rtick' must positive integer.")}
    barloc <- barloc[1]
    if(is.null(barloc)){
      barloc <- ifelse(ncol(space)==2L,"right","backright")
    } else if(ncol(space)==2L){
      if(!any(barloc==c("right","left"))){"Input 'barloc' must be set to 'right' or 'left' for 2d input."}
    } else {
      if(!any(barloc==c("backright","backleft","frontright","frontleft"))){"Input 'barloc' must be set to 'backright', 'backleft', 'frontright', or 'frontleft', for 3d input."}
    }
    if(is.null(mar[1]) & ncol(space)==2L){
      if(barloc=="right"){mar <- c(4,2,4,5)} else{mar <- c(4,5,4,2)}
    }
    oldmar <- par()$mar
    on.exit(par(mar=oldmar))
    par(mar=mar)
    
    ### map voltages to colors
    if(is.null(mycolors[1])){
      mycolors <- colorRampPalette(c("blueviolet","blue","cyan","green","yellow","orange","red"))(ncolor)
    } else {
      mycolors <- colorRampPalette(mycolors)(ncolor)
    }
    voltcolors <- voltcol(voltage,vlim[1],vlim[2],mycolors,ncolor)
    
    ### create plot
    if(ncol(space)==2L){
      
      # check labels
      if(is.null(xyzlab)){
        xyzlab <- rep("",2)
      } else if(length(xyzlab)!=2L){
        stop("Input 'xyzlab' must be two element vector for 2D plots.")
      }
      
      # check colorbar location
      if(barloc=="right"){
        xlim <- c(-12.5,22.5)
        aidx <- 4
        adjx <- 12/35
        blim <- c(20.5,22.5)
      } else {
        xlim <- c(-22.5,12.5)
        aidx <- 2
        adjx <- 23/35
        blim <- c(-22.5,-20.2)
      }
      
      # plot head
      rad <- 12.5
      xx <- rad*cos(seq(0,2*pi,length.out=360))
      yy <- rad*sin(seq(0,2*pi,length.out=360))
      if(head){
        plot(xx,yy,asp=1,xlim=xlim,ylim=c(-12.5,12.5),type="l",
             xlab=xyzlab[1],ylab=xyzlab[2],axes=plotaxes,main="",...)
      } else {
        plot(1,1,asp=1,xlim=xlim,ylim=c(-12.5,12.5),type="n",
             xlab=xyzlab[1],ylab=xyzlab[2],axes=plotaxes,main="",...)
      }
      title(main,cex.main=cex.main,adj=adjx)
      if(nose){
        lines(c(xx[81],0),c(yy[81],rad*1.175))
        lines(c(-xx[81],0),c(yy[81],rad*1.175))
      }
      if(ears){
        xx <- 0.5*cos(seq(0,2*pi,l=360))
        yy <- 2.5*sin(seq(0,2*pi,l=360))
        lines(xx-13,yy)
        lines(xx+13,yy)
      }
      
      # plot color-mapped voltages
      points(space[,1],space[,2],pch=19,cex=cex.point,col=voltcolors)
      
      # plot colorbar (if needed)
      if(colorbar){
        
        # draw axis for colorbar
        if(is.null(colorlab)){colorlab <- expression("Voltage ("*mu*"V)")}
        ticklab <- round(seq(vlim[1],vlim[2],length.out=nctick),rtick)
        nctick <- length(ticklab)
        tickloc <- seq(-12.5,12.5,length.out=nctick)
        axis(aidx,cex.axis=cex.axis,at=tickloc,labels=as.character(ticklab))
        mtext(colorlab,side=aidx,line=3,cex=cex.lab*par()$cex)
        
        # draw colorbars
        scales <- ncolor/25
        for (ii in 1:ncolor) {
          y <- ((ii-1)/scales) - 12.5
          rect(blim[1],y,blim[2],y+1/scales,col=mycolors[ii],border=NA)
        }
        
      } # end if(colorbar)
      
    } else {
      
      # check labels
      if(is.null(xyzlab)){
        xyzlab <- rep("",3)
      } else if(length(xyzlab)!=3L){
        stop("Input 'xyzlab' must be three element vector for 3D plots.")
      }
      
      # load 3d dense cap and mesh
      eegdense <- eegmesh <- NULL
      data(eegdense,eegmesh,envir=environment())
      
      # map input coordinates to dense coordinates
      dmat <- matrix(0,977,vlen)
      for(j in 1:3){
        dmat <- dmat + (matrix(eegdense[,j],977,vlen)-t(matrix(space[,j],vlen,977)))^2
      }
      cidx <- apply(dmat,1,which.min)
      dencol <- voltcolors[cidx]
      
      # assign colors to mesh and plot
      itcol <- dencol[c(eegmesh$it)]
      eegmesh$material$color <- itcol
      plot3d(1,1,1,type="n",xlim=c(-12.5,12.5),ylim=c(-12.5,12.5),
             zlim=c(-17.5,12.5),axes=FALSE,xlab=xyzlab[1],
             ylab=xyzlab[2],zlab=xyzlab[3])
      shade3d(eegmesh)
      if(head){
        eeghead <- NULL
        data(eeghead,envir=environment())
        eeghead$material$color <- rep(col.head,length(eeghead$material$color))
        shade3d(eeghead)
      }
      
      # plot colorbar (if needed)
      if(colorbar){
        
        mybar <- colorbar3d(c(-17.5,12.5),mycolors,ncolor)
        
        # get barloc
        if(barloc=="backright"){
          mybar <- scale3d(mybar,1,-1,1)
          axloc <- "z+-"
        } else if(barloc=="backleft"){
          mybar <- scale3d(mybar,-1,-1,1)
          axloc <- "z--"
        } else if(barloc=="frontleft"){
          mybar <- scale3d(mybar,-1,1,1)
          axloc <- "z-+"
        } else {
          axloc <- "z++"
        }
        
        # get labels and plot
        shade3d(mybar)
        if(is.null(colorlab)){colorlab <- "Voltage (mcV)"}
        ticklab <- round(seq(vlim[1],vlim[2],length.out=nctick),rtick)
        nctick <- length(ticklab)
        tickloc <- seq(-17.5,12.5,length.out=nctick)
        axis3d(axloc,cex.axis=cex.axis,at=tickloc,labels=as.character(ticklab))
        mtext3d(colorlab,edge=axloc,line=3,cex=cex.lab*par()$cex,srt=6)
        
      } # end if(colorbar)
      
    } # end if(ncol(space)==2L)
  
}
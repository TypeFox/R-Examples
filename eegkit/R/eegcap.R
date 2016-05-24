eegcap <- 
  function(electrodes="10-10",type=c("3d","2d"),plotlabels=TRUE,
           plotaxes=FALSE,main="",xyzlab=NULL,cex.point=NULL,
           col.point=NULL,cex.label=NULL,col.label=NULL,nose=TRUE,
           ears=TRUE,head=TRUE,col.head="AntiqueWhite",index=FALSE,
           plt=c(0.03,0.97,0.03,0.97),...){
    ###### Plots EEG Cap with Selected Electrodes (2D or 3D)
    ###### Nathaniel E. Helwig (helwig@umn.edu)
    ###### Last modified: February 16, 2015
    
    ### initial checks
    eegcoord <- NULL
    data(eegcoord,envir=environment())
    enames <- rownames(eegcoord)
    if(electrodes[1]=="10-10"){
      eegidx <- 1:87
    } else if (electrodes[1]=="10-20"){
      eegidx <- match(c("A1","A2","FP1","FP2", "FPZ", "F7","F3","FZ",
                        "F4","F8", "NZ","T7","C3","CZ","C4","T8",
                        "P7","P3","PZ","P4","P8","O1","O2", "OZ"),enames)      
    } else {
      eegidx <- match(toupper(electrodes),enames)
    }
    type <- type[1]
    if(!any(type==c("2d","3d"))){stop("Incorrect 'type' input.")}
    
    ### plot 2d or 3d electrodes
    if(type=="2d"){
      
      # get old par and set new par
      oldplt <- par()$plt
      on.exit(par(plt=oldplt))
      par(plt=plt)
      
      # check labels
      if(is.null(xyzlab)){
        xyzlab <- rep("",2)
      } else if(length(xyzlab)!=2L){
        stop("Input 'xyzlab' must be two element vector for 2D plots.")
      }
      
      # plot head, nose, and ears
      rad <- 12.5
      xx <- rad*cos(seq(0,2*pi,length.out=360))
      yy <- rad*sin(seq(0,2*pi,length.out=360))
      if(head){
        plot(xx,yy,asp=1,xlim=c(-16,16),ylim=c(-16,16),type="l",
             xlab=xyzlab[1],ylab=xyzlab[2],axes=plotaxes,main=main,...)
      } else {
        plot(1,1,asp=1,xlim=c(-16,16),ylim=c(-16,16),type="n",
             xlab=xyzlab[1],ylab=xyzlab[2],axes=plotaxes,main=main,...)
      }
      if(nose){
        lines(c(xx[81],0),c(yy[81],rad*1.175))
        lines(c(-xx[81],0),c(yy[81],rad*1.175))
      }
      if(ears){
        xx <- 0.5*cos(seq(0,2*pi,l=360))
        yy <- 2.5*sin(seq(0,2*pi,l=360))
        lines(xx-13,yy-1.5)
        lines(xx+13,yy-1.5)
      }
      
      # plot electrodes
      if(plotlabels){
        if(is.null(cex.point)){cex.point <- 2.75}
        if(is.null(cex.label)){cex.label <- 0.5}
        if(is.null(col.label[1])){col.label <- "blue"}
        if(is.null(col.point[1])){col.point <- "green"}
        points(eegcoord[eegidx,4],eegcoord[eegidx,5],cex=cex.point,col=col.point,pch=19)
        points(eegcoord[eegidx,4],eegcoord[eegidx,5],cex=cex.point,pch=21)
        text(eegcoord[eegidx,4],eegcoord[eegidx,5],labels=enames[eegidx],cex=cex.label,col=col.label)
      } else {
        if(is.null(cex.point)){cex.point <- 1}
        if(is.null(col.point[1])){col.point <- "green"}
        points(eegcoord[eegidx,4],eegcoord[eegidx,5],cex=cex.point,col=col.point,pch=19)
        points(eegcoord[eegidx,4],eegcoord[eegidx,5],cex=cex.point,pch=21)
      }
      
    } else {
      
      # check labels
      if(is.null(xyzlab)){
        xyzlab <- rep("",3)
      } else if(length(xyzlab)!=3L){
        stop("Input 'xyzlab' must be three element vector for 3D plots.")
      }
      
      # plot electrodes
      if(is.null(cex.point)){cex.point <- 10}
      if(is.null(col.point[1])){col.point <- "green"}
      plot3d(eegcoord[eegidx,1],eegcoord[eegidx,2],eegcoord[eegidx,3],
             xlab=xyzlab[1],ylab=xyzlab[2],zlab=xyzlab[3],type="n",...)
      if(plotaxes){axes3d()}
      points3d(eegcoord[eegidx,1],eegcoord[eegidx,2],eegcoord[eegidx,3],
               size=cex.point,col=col.point,pch=19)
      if(plotlabels){
        if(is.null(cex.label)){cex.label <- 1.25}
        if(is.null(col.label[1])){col.label <- "blue"}
        text3d(eegcoord[eegidx,1]*1.05,eegcoord[eegidx,2]*1.05,eegcoord[eegidx,3]*1.05,
               texts=enames[eegidx],cex=cex.label,col=col.label)
      } else {
        plot3d(eegcoord[eegidx,1],eegcoord[eegidx,2],eegcoord[eegidx,3],
               xlab=xyzlab[1],ylab=xyzlab[2],zlab=xyzlab[3],size=cex.point,
               col=col.point,pch=19,...)
        if(plotaxes){axes3d()}
      }
      if(head){
        eeghead <- NULL
        data(eeghead,envir=environment())
        eeghead$material$color <- rep(col.head[1],length(eeghead$material$color))
        shade3d(eeghead)
      }
      par3d(scale=rep(1,3))
      
    } # end if(type=="2D")
  
    if(index){return(eegidx)}
    
}
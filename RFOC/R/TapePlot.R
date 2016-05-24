TapePlot<-function(TapeList=list(),add=FALSE, ann=TRUE,pcol= c(grey(0.85) , grey(0.95))  )
  {

    if(is.null(TapeList))
      {
        TapeList=TapeBase()
      }

    if(add==FALSE)
      {
        plot(range(c(  TapeList$Left1$x, TapeList$right2$x)), range(c(TapeList$top1$y, TapeList$bot1$y)),
             type='n', asp=1, axes=FALSE, ann=FALSE)
      }


    ## Underneath: filled Polygon patches 
    polygon(TapeList$POLYh1$x, TapeList$POLYh1$y, col=pcol[1] , border=NA)
    
    polygon(TapeList$POLYh2$x, TapeList$POLYh2$y, col=pcol[2] , border=NA)
    for(i in 1:length(TapeList$LONSp1))
      {
        
        segments(TapeList$LONSp1[[i]]$x, TapeList$LONSp1[[i]]$y, TapeList$LONSp2[[i]]$x, TapeList$LONSp2[[i]]$y)
      }

    for(i in 1:length(TapeList$LATSp1))
      {
        
        lines(TapeList$LATSp1[[i]]$x, TapeList$LATSp1[[i]]$y)
      }

    points(TapeList$PTSh1$x, TapeList$PTSh1$y,  pch=21, bg="white" )

     if(ann)
      {
        f1 = TapeList$d1$a1>0
        text(TapeList$PTSh1$x[f1], TapeList$PTSh1$y[f1], TapeList$d1$name[f1], pos=4)
        text(TapeList$PTSh1$x[!f1], TapeList$PTSh1$y[!f1], TapeList$d1$name[!f1], pos=2)
      }



     lines(TapeList$HOZh1$x, TapeList$HOZh1$y, lty=2, lwd=2)
 
 lines(TapeList$VERTh1$x, TapeList$VERTh1$y, lty=2, lwd=2)

    ###  add in a bold dotted line for LVD-1
    lines(TapeList$PATCH11$x, TapeList$PATCH11$y, lty=2, lwd=2)
###  add in a bold dotted line for LVD-2
    
    lines(TapeList$PATCH21$x, TapeList$PATCH21$y, lty=2, lwd=2)
###  connect up the Crack lines
    
    lines(TapeList$CRACKh1$x, TapeList$CRACKh1$y, lty=2, lwd=2)


  }

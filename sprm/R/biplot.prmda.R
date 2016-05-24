biplot.prmda <-
  function(x, comps=c(1,2), colors=list(scores1="orange",scores2="darkgreen", loadings="red",background="gray84"),textsize=6,arrowshapes=c(25,.03), labelpos=0.35, ...)    
{
    #  require(ggplot2)
    #  require(grid)
    V1 <- V2 <- V3 <- V4 <- V5 <- NULL
    if(!(class(x)=="prmda")){
      stop("The SPRM-DA biplot function can only be applied to sprmda class objects")
    }
    if (x$inputs$a<max(comps)){
      stop("Total number of componets in the model is smaller than specified components in comps.")
    }
    c1 <- comps[1]
    c2 <- comps[2]
    colors$scores <- ifelse(x$inputs$y0==1,colors$scores1, colors$scores2)
	
    plotscores <- as.data.frame(x$scores[,comps])
    colnames(plotscores) <- c(paste("Comp",c1,sep=""),paste("Comp",c2,sep=""))
    plotscalingfactors <- (apply(x$scores[,comps],2,max)-apply(x$scores[,comps],2,min))/(apply(x$loadings[,comps],2,max)-apply(x$loadings[,comps],2,min))*0.8
    plotloadings <- as.data.frame(x$loadings[,comps])
    colnames(plotloadings) <- c("V1", "V2")
    plotloadings[,1] <- plotloadings[,1]*plotscalingfactors[1]
    plotloadings[,2] <- plotloadings[,2]*plotscalingfactors[2]
    npl <- nrow(plotloadings)
    plotloadings$V3 <- rep(0,npl)
    plotloadings$V4 <- plotloadings$V1+sign(plotloadings$V1)* labelpos
    plotloadings$V5 <- plotloadings$V2+sign(plotloadings$V2)* labelpos
    plotty <- ggplot()
    plotty <- plotty + geom_text(data=plotscores, aes_string(x=colnames(plotscores)[1],y=colnames(plotscores)[2]), size=textsize,label=rownames(plotscores),color=colors$scores) 
    plotty <- plotty + labs(title=paste(names(x$YMeans)," Sparse PRM Components ", c1," and ",c2," Biplot"),x=colnames(plotscores)[1],y=colnames(plotscores)[2])
    plotty <- plotty + geom_segment(data=plotloadings,aes(x=V3,y=V3,xend=V1,yend=V2),color=colors$loadings,arrow=ggplot2::arrow(angle=arrowshapes[1],length=ggplot2::unit(arrowshapes[2],"npc")),na.rm=TRUE)
    plotty <- plotty + geom_text(data=plotloadings,aes(x=V4,y=V5),color=colors$loadings,label=rownames(plotloadings),size=textsize,na.rm=TRUE)
    plotty <- plotty + theme(panel.background=element_rect(fill=colors$background),plot.title=element_text(size=rel(1.5),face="bold"))
    print(plotty)
  }
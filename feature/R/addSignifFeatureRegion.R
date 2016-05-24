########## R function: addSignifFeatureRegion ##########

## For adding significant features region to a KDE plot.

## Last changed: 04 NOV 2005

addSignifFeatureRegion <- function(d,gridsize,SignifFeature,plot.inds,
    featureCol,dest,range.x.plot,add.bars=TRUE, trans.alpha)
{
  if (d==1)
  {
    SignifFeature.inds <- (1:gridsize)[SignifFeature==TRUE]
    SignifFeature.inds <- intersect(SignifFeature.inds,plot.inds[[1]])
    SGlen <- length(SignifFeature.inds)
    diff.vec <- diff(SignifFeature.inds)
    jump.inds <- (1:length(diff.vec))[diff.vec!=1]
    num.jumps <- length(jump.inds)
    
    if (num.jumps==0) lines(dest$x.grid[[1]][SignifFeature.inds], dest$est[SignifFeature.inds],col=featureCol,lwd=3)
    
    if (num.jumps>0)
    {
      curr.inds <- SignifFeature.inds[1:jump.inds[1]]
      lines(dest$x.grid[[1]][curr.inds],dest$est[curr.inds],
            col=featureCol,lwd=3)
      if (num.jumps>1) 
      { 
        for (j in 2:length(jump.inds))
        {
          curr.inds <- SignifFeature.inds[(jump.inds[j-1]+1):jump.inds[j]]
          lines(dest$x.grid[[1]][curr.inds],dest$est[curr.inds], col=featureCol,lwd=3)
        }
      }
      curr.inds <- SignifFeature.inds[(max(jump.inds)+1):SGlen]
      lines(dest$x.grid[[1]][curr.inds],dest$est[curr.inds], col=featureCol,lwd=3)
    }
  }

  if (d==2)
  {
    x.grid.1 <- dest$x.grid[[1]] ; x.grid.2 <- dest$x.grid[[2]]
    x.mesh <- expand.grid(x.grid.1,x.grid.2)
    inds.on <- (1:nrow(x.mesh))[as.vector(SignifFeature)]
    inds.on <- intersect(inds.on, (1:nrow(x.mesh))[x.mesh[,1]>=range.x.plot[[1]][1]])
    inds.on <- intersect(inds.on,(1:nrow(x.mesh))[x.mesh[,1]<=range.x.plot[[1]][2]])
    inds.on <- intersect(inds.on, (1:nrow(x.mesh))[x.mesh[,2]>=range.x.plot[[2]][1]])
    inds.on <- intersect(inds.on, (1:nrow(x.mesh))[x.mesh[,2]<=range.x.plot[[2]][2]])

    if (add.bars)
      points(x.mesh[inds.on,1],x.mesh[inds.on,2],pch=".",col=featureCol, cex=2)

    contour(x.grid.1[plot.inds[[1]]],x.grid.2[plot.inds[[2]]],
            SignifFeature[plot.inds[[1]],plot.inds[[2]]],add=TRUE,
            col=featureCol,levels=0.5,lwd=3,drawlabels=FALSE)  
  }

  if (d==3)
  {
    x.gd.1 <- dest$x.grid[[1]] ; x.gd.2 <- dest$x.grid[[2]]
    x.gd.3 <- dest$x.grid[[3]]

    if (!all(SignifFeature==FALSE))
      contour3d(SignifFeature,level=0.5,x=x.gd.1,color=featureCol,y=x.gd.2,z=x.gd.3,alpha=trans.alpha,add=TRUE)
  }

  ## for d==4, only significant curvature is calculated and plotted
  if (d==4)
  {
    x.gd.1 <- dest$x.grid[[1]] ; x.gd.2 <- dest$x.grid[[2]]
    x.gd.3 <- dest$x.grid[[3]] 

    if (!all(SignifFeature==FALSE))
    {  
      sig.levs <- vector()
      for (i in 1:gridsize[4])
        if (any(SignifFeature[,,,i])) sig.levs <- c(sig.levs,i)
      
      if (length(sig.levs)>0)
      {  
        ## somewhat arbitrary way of choosing which indices to plot...
        ##sig.levs <- as.integer(pretty(sig.levs, n=5)[-c(1,2)])
        if (length(sig.levs)>4)
          sig.levs <- sig.levs[-c(1,length(sig.levs))]
        num.levs <- length(sig.levs)
        
        if (length(trans.alpha)==2)
          alph.vec <- seq(trans.alpha[1], trans.alpha[2], length=num.levs)
        else if (length(trans.alpha)==1)
          alph.vec <- rep(trans.alpha, num.levs)
      
        j <- 1
        for (i in sig.levs)
        {
          contour3d(SignifFeature[,,,i],level=0.5,x=x.gd.1,color=featureCol,
                      y=x.gd.2,z=x.gd.3, alpha=alph.vec[j],add=TRUE)
          j <- j + 1
        }
      }
    }
  }
}

########## End of addSignifFeatureRegion ##########

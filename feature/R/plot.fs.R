plot.fs <-  function(x, ..., xlab, ylab, zlab, xlim, ylim, zlim, add=FALSE,
           addData=FALSE, scaleData=FALSE, addDataNum=1000,
           addKDE=TRUE, jitterRug=TRUE,  
           addSignifGradRegion=FALSE, addSignifGradData=FALSE,
           addSignifCurvRegion=FALSE, addSignifCurvData=FALSE,
           addAxes3d=TRUE,
           densCol, dataCol="black", gradCol="green4", curvCol="blue",
           axisCol="black", bgCol="white",
           dataAlpha=0.1, gradDataAlpha=0.3,
           gradRegionAlpha=0.2, curvDataAlpha=0.3, curvRegionAlpha=0.3) 
{
  fs <- x

  x <- as.matrix(fs$x)
  d <- ncol(x)
  n <- nrow(x)
  h <- fs$bw
  names.x <- fs$names
  
  if (d >1) gridsize <- dim(fs$fhat$est) 
  else gridsize <- length(fs$fhat$est) 
  
  ## Determine default axis labels.

  if (missing(xlab)) xlab <- NULL
  if (missing(ylab)) ylab <- NULL
  if (missing(zlab)) zlab <- NULL
  labs <- dfltLabs(d,names.x,xlab,ylab,zlab)
  xlab <- labs$xlab ; ylab <- labs$ylab ; zlab <- labs$zlab
 
  dest <- fs$fhat
  ESS <- n*dest$est*prod(h)*(sqrt(2*pi)^d)
  SigESS <- ESS >= 5
  
  ## random sample of data points used for display
  nsamp <- min(addDataNum, n)
  
  if (nsamp < n)
  {
    rand.inds <- 1:nsamp 
    x.rand <- as.matrix(x[rand.inds,])
  }
  else
    x.rand <- x
  
  if (missing(xlim))
    if (d==1)
      xlim <- c(min(x)-h[1],max(x)+h[1])
    else
      xlim <- c(min(x[,1])-h[1],max(x[,1])+h[1])

  if (missing(ylim))
    if (d==1)
      ylim <- c(0,1.5)*max(dest$est)
    else if (d>1)
      ylim <- c(min(x[,2])-h[2],max(x[,2])+h[2])
  
  if (missing(zlim) & d>2)
    zlim <- c(min(x[,3])-h[3],max(x[,3])+h[3])
  
  if (d==1)
    lims <- list(xlim)
  if (d==2)
    lims <- list(xlim, ylim)
  if (d==3)
    lims <- list(xlim, ylim, zlim)
  if (d==4)
    lims <- list(xlim, ylim, zlim, c(min(x[,4])-h[4],max(x[,4])+h[4]))

  plot.inds <- list()
  for (id in 1:d)
  {
    plot.inds.l <- (1:gridsize[id])[dest$x.grid[[id]]>=lims[[id]][1]]
    plot.inds.u <- (1:gridsize[id])[dest$x.grid[[id]]<=lims[[id]][2]]
    plot.inds[[id]] <- intersect(plot.inds.l,plot.inds.u)
  }
  
  if (missing(densCol))
    if (d==1)
      densCol <- "DarkOrange"
    else if (d==2)
      densCol <- rev(heat.colors(1000))
    else if (d==3)
      densCol <- rev(heat.colors(3))
 
  if (d==1)
  {
    par(bg=bgCol)
    if (addKDE)
    {
      plot(dest$x.grid[[1]][plot.inds[[1]]], dest$est[plot.inds[[1]]],
           type="n",bty="l" ,col=densCol, lwd=2, xlim=xlim, ylim=ylim,
           xlab=xlab,ylab="kernel density estimate")
    
      lines(dest$x.grid[[1]][plot.inds[[1]]],dest$est[plot.inds[[1]]],
            bty="l",col=densCol,lwd=2)
    }

    if (addData)
    {
      if (jitterRug) x.rug <- jitter(x.rand)
      else x.rug <- x.rand
      rug(x.rug)
    }  
  }
  else if (d==2)
  {
    par(bg=bgCol)
    x.grid.1 <- dest$x.grid[[1]] ; x.grid.2 <- dest$x.grid[[2]]

    if (addKDE)
      image(x.grid.1[plot.inds[[1]]],x.grid.2[plot.inds[[2]]],
            dest$est[plot.inds[[1]],plot.inds[[2]]],col=densCol,
            xlim=xlim, ylim=ylim, xlab=xlab,ylab=ylab,bty="n")
    if (!add & !addKDE)
      plot(x.grid.1, x.grid.2, xlim=xlim, ylim=ylim, xlab=xlab,ylab=ylab,type="n")
    
    box()

    if (addData)
      points(x.rand, col=dataCol)
  }
  else if (d==3)
  {
    if (!add)
    {
      ##clear3d()
      ##bg3d(bgCol)
      ##pop3d(type="lights")
      ##light3d(theta=0, phi=30)
      
      ##material3d(alpha=1)
      ##material3d(back="fill")

      plot3d(mean(xlim), mean(ylim), mean(zlim), xlab=xlab, ylab=ylab, zlab=zlab, xlim=xlim, ylim=ylim, zlim=zlim, axes=addAxes3d, box=addAxes3d, colors="transparent", alpha=0)
     
    }

    if (addKDE)
    {
        kde.temp <- kde(x, H=diag(h^2), binned=TRUE, gridsize=rep(31,3), compute.cont=TRUE, approx.cont=TRUE)
        plot(kde.temp, box=FALSE, axes=FALSE, add=TRUE)
    }
    if (addData)
      points3d(x.rand[,1],x.rand[,2],x.rand[,3],size=3,color=dataCol, alpha=dataAlpha)
     bg3d(bgCol)
  }

  SignifGradRegion.mat <- fs$grad
  SignifCurvRegion.mat <- fs$curv

  if (!is.null(SignifGradRegion.mat))
  {
    SignifGradData.mat <- SignifFeatureData(x.rand, d, dest,SignifGradRegion.mat)
    if (addSignifGradRegion)
      addSignifFeatureRegion(d,gridsize,SignifGradRegion.mat,plot.inds,gradCol, dest,lims, trans.alpha=gradRegionAlpha)
    if (addSignifGradData)
      addSignifFeatureData(x.rand,SignifGradData.mat,gradCol, trans.alpha=gradDataAlpha)
  }
  if (!is.null(SignifCurvRegion.mat))
  {
    SignifCurvData.mat <- SignifFeatureData(x.rand, d, dest,SignifCurvRegion.mat)
    if (addSignifCurvRegion)
      addSignifFeatureRegion(d,gridsize,SignifCurvRegion.mat,plot.inds,curvCol, dest,lims, trans.alpha=curvRegionAlpha)
    if (addSignifCurvData)
      addSignifFeatureData(x.rand,SignifCurvData.mat,curvCol, trans.alpha=curvDataAlpha)
  }
 
  invisible()
}
  

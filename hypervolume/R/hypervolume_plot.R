
plot.Hypervolume <- function(x, ...)
{
  templist = new("HypervolumeList")
  templist@HVList=list(x)
  plot.HypervolumeList(templist, ...)
}

extendrange <- function(x,factor=0.5)
{
  xmin <- min(x,na.rm=T)
  xmax <- max(x,na.rm=T)
  
  xminf <- xmin - (xmax - xmin)*factor
  xmaxf <- xmax + (xmax - xmin)*factor
  
  result <- c(xminf, xmaxf)
  
  return(result)
}

plot.HypervolumeList <- function(x, npmax_data = 1000, npmax_random = 2000, 
                                 colors=rainbow(length(x@HVList),alpha=0.8), names=NULL, 
                                 reshuffle=TRUE, showrandom=TRUE, showdensity=TRUE,showdata=TRUE,darkfactor=0.5,
                                 cex.random=0.5,cex.data=0.75,cex.axis=0.75,cex.names=1.0,cex.legend=0.75,
                                 legend=TRUE, varlims=NULL, showcontour=TRUE, contour.lwd=1, contour.filled=FALSE,contour.filled.alpha=0.5,contour.factor=0.05,
                                 showcentroid=TRUE, cex.centroid=3,
                                 pairplot=TRUE,whichaxes=NULL,...)
{
  sapply(x@HVList, function(z)
  {
    cat(sprintf("Showing %d random points of %d for %s\n",min(nrow(z@RandomUniformPointsThresholded), npmax_random), nrow(z@RandomUniformPointsThresholded), z@Name))
    if (showdata && length(z@Data) > 0)
    {
      npd <- ifelse(all(is.nan(z@Data)), 0, nrow(z@Data))
      cat(sprintf("Showing %d data points of %d for %s\n",min(npmax_data, npd), npd, z@Name))
    }    
    
  })
  
  alldims = sapply(x@HVList, function(z) { z@Dimensionality })
  allnames = sapply(x@HVList, function(z) { z@Name })
  stopifnot(all(alldims[1] == alldims))
  
  all <- NULL
  alldata <- NULL
  for (i in 1:length(x@HVList))
  {
    ivals = sample(nrow(x@HVList[[i]]@RandomUniformPointsThresholded), min(c(npmax_random, nrow(x@HVList[[i]]@RandomUniformPointsThresholded))))
    subsampledpoints = data.frame(x@HVList[[i]]@RandomUniformPointsThresholded[ivals,,drop=FALSE])
    densityvals = x@HVList[[i]]@ProbabilityDensityAtRandomUniformPoints[ivals]
    
    if (nrow(subsampledpoints) > 0)
    {  
      subsampledpoints = cbind(subsampledpoints, ID=rep(i, nrow(subsampledpoints)), Density=densityvals/max(densityvals,na.rm=T))
    
      all <- rbind(all, subsampledpoints)
    }
    
    thisdata=x@HVList[[i]]@Data
    alldata <- rbind(alldata, cbind(thisdata, ID=rep(i,nrow(thisdata))))
  }  
  all <- unique(all)
  alldata <- as.data.frame(alldata)
  alldata <- alldata[sample(nrow(alldata), min(c(npmax_data, nrow(alldata)))),]
  
  if (is.null(all))
  {
    stop('Nothing to plot.')
  }
  
  if (reshuffle==TRUE)
  {
    all <- all[sample(nrow(all)),] # reorder to shuffle colors
  }
  
  if (is.null(names))
  {
    names = names(all)[1:(ncol(all)-2)]
  }  
  
  if (!is.null(varlims) & !is.list(varlims))
  {
    varlimlist = vector('list',ncol(all)-2)
    for (i in 1:length(varlimlist))
    {
      varlimlist[[i]] <- varlims
    }
    varlims = varlimlist
  }
  
  colorlist <- colors[all$ID]
  alphavals <- all$Density
  if (showdensity)
  {
    colorlist <- rgb2rgba(colorlist, alphavals)
  }
  
  colorlistdata = colors[alldata$ID]
  colorlistdata <- rgb2rgbdark(colorlistdata, darkfactor)
  
  
  if (ncol(all) < 2)
  {
    stop('Plotting only available in n>=2 dimensions.')
  }
  
  if (pairplot)
  {
    op = par(no.readonly = T)
    
    par(mfrow=c(ncol(all)-2, ncol(all)-2))
    par(mar=c(0,0,0,0))
    
    for (i in 1:(ncol(all)-2))
    {
      for (j in 1:(ncol(all)-2))  
      {
        if (j > i)
        {
          # set up axes with right limits
          plot(all[,j], all[,i],type="n",axes=F,xlim=varlims[[j]], ylim=varlims[[i]])
  
          
          
          # draw random points
          if(showrandom==TRUE)
          {
            points(all[,j], all[,i], col=colorlist,cex=cex.random,pch=16)
          }
          
          # show data
          if (showdata & nrow(alldata) > 0)
          {
            points(alldata[,j], alldata[,i], col=colorlistdata,cex=cex.data,pch=16)
          }
          
          if (showcentroid == TRUE)
          {
            for (whichid in 1:length(unique(all$ID)))
            {
              allss <- subset(all, all$ID==whichid)
              centroid_x <- mean(allss[,j],na.rm=T) + rnorm(1)*diff(range(all[,j]))*0.01
              centroid_y <- mean(allss[,i],na.rm=T) + rnorm(1)*diff(range(all[,i]))*0.01
              
              # draw point
              points(centroid_x, centroid_y, col=colors[whichid],cex=cex.centroid,pch=16)
              # add a white boundary for clarity
              points(centroid_x, centroid_y, col='white',cex=cex.centroid,pch=1,lwd=1.5)
            }
          }
          
          
          # calculate contours
          if (showcontour==TRUE)
          {
            if (contour.filled==TRUE)
            {
              # draw shaded centers
              for (whichid in 1:length(unique(all$ID)))
              {
                allss <- subset(all, all$ID==whichid)
                
                if (nrow(allss) > 0)
                {     
                  contourx <- allss[,j]
                  contoury <- allss[,i]
                  
                  kde2dresults <- kde2d(contourx, contoury, n=50,lims=c(extendrange(contourx),extendrange(contoury)))
                  
                  .filled.contour(kde2dresults$x,kde2dresults$y, kde2dresults$z,
                                  col=c(NA,rgb2rgba(colors[whichid],contour.filled.alpha),NA),
                                  levels=c(0,min(kde2dresults$z)+diff(range(kde2dresults$z))*contour.factor,max(kde2dresults$z)))
                }
              }
            }
            
            # draw edges
            for (whichid in 1:length(unique(all$ID)))
            {
              allss <- subset(all, all$ID==whichid)
              
              if (nrow(allss) > 0)
              {         
                contourx <- allss[,j]
                contoury <- allss[,i]
                
                if (length(contourx) > 1 & length(contoury) > 1)
                {
                  kde2dresults <- kde2d(contourx, contoury, n=50,lims=c(extendrange(contourx),extendrange(contoury)))
                  
                  contour(kde2dresults,
                          col=colors[whichid],
                          levels=min(kde2dresults$z)+diff(range(kde2dresults$z))*contour.factor,
                          lwd=contour.lwd,
                          drawlabels=FALSE,add=TRUE)
                }
              }
            }
          }
          
          box()
        }
        else if (j == i)
        {
          plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),axes=F)
          text(0.5, 0.5, names[j],cex=cex.names)
        }
        else if (j==1 & i == (ncol(all) - 2))
        {
          plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),axes=F)
          
          if (legend == TRUE)
          {
            legend('topleft',legend=allnames,text.col=colors,bty='n',cex=cex.legend)
          }
        }
        else
        {
          plot(0,0,type="n",axes=F)    
        }
        
        if (j==i+1)
        {
          axis(side=1,cex.axis=cex.axis)
          axis(side=2,cex.axis=cex.axis)
        }
      }
    }  
    par(op)
  }
  else
  {
    if (is.null(whichaxes))
    {
      whichaxes=1:3  
    }
    if (is.null(names))
    {
      names <- names(data)
    }
    if(length(whichaxes)!=3) { stop('Must specify three axes') }
    
    if (all(is.numeric(whichaxes)))
    {
      axesnames <- names
    }

    rgl::plot3d(all[,whichaxes],col=colorlist,xlab=axesnames[1], ylab=axesnames[2], zlab=axesnames[3], xlim=varlims[[1]],ylim=varlims[[2]],zlim=varlims[[3]],size=cex.random,type='p',expand=1.05)
    
    if (legend==TRUE)
    {
      for (i in 1:length(allnames))
      {
        rgl::mtext3d(allnames[i],edge='x-+',line=1+i*cex.legend*1.25,color=colors[i],cex=cex.legend)  
      }
    }
    
    if (showdata)
    {
      if (!any(is.nan(as.matrix(alldata[,whichaxes]))))
      {
        rgl::points3d(x=alldata[,whichaxes[1]], y=alldata[,whichaxes[2]], z=alldata[,whichaxes[3]], col=colorlistdata,cex=cex.data,pch=16)
      }
    }
    
    if (showcentroid == TRUE)
    {
      for (whichid in 1:length(unique(all$ID)))
      {
        allss <- subset(all, all$ID==whichid)
        centroid_1 <- mean(allss[,whichaxes[1]],na.rm=T)
        centroid_2 <- mean(allss[,whichaxes[2]],na.rm=T)
        centroid_3 <- mean(allss[,whichaxes[3]],na.rm=T)
        
        # draw point
        rgl::points3d(x=centroid_1, y=centroid_2, z=centroid_3, col=colors[whichid],cex=cex.centroid,pch=16)
      }
    }
  }
}  

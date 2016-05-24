##-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-##
## Compute and visualize the high dimensional function
##-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-##
plot.msc.svm <- function(x, drawStdDev=FALSE, span=0.5, nsamples=50, plot=TRUE, colorMap=0,...) {
  plot.msc(x, drawStdDev, span, nsamples)
}

plot.msc.kd <- function(x, drawStdDev=FALSE, span=0.5, nsamples=50, plot=TRUE, colorMap=0,...) {
  plot.msc(x, drawStdDev, span, nsamples)
}

plot.msc <- function(x, drawStdDev=FALSE, span=0.5, nsamples=50, plot=TRUE, colorMap=0,...) {

  ## Compute the geometry
  geom <- geometry.mscPlot(x, span=span, nsamples=nsamples)
  
  ## Compute the vis scene
  scene <- scene.mscPlot(geom, colorMap=colorMap)
  
  ## Create the mscPlot object
  obj <- structure(list(geom=geom, scene=scene),class="mscPlot")

  ## Create the mscPlot
  if(plot)
    obj <- plot.mscPlot(obj, drawStdDev)

  ## Return the object
  invisible(obj)
}

##-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-##
## Plot the rgl scene
##-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-##
plot.mscPlot <- function(x, drawStdDev=FALSE, axesOn=TRUE, ...) {
  
  ## Initialize the 3D scene
  if(length(x$dev) == 0 || rgl.cur() == 0) {
    x$dev <- open3d(mouseMode=c("none","trackball", "zoom"), windowRect=c(0, 0, 512, 512))
    view3d(phi=0, theta=0, fov=15)
    light3d(phi=-45, theta=-45)
    light3d(phi=45, theta=45)
    bg3d(sphere=FALSE, color="white")
  }
  else
    rgl.clear()
  
  ## Draw the spheres
  spheres3d(x$geom$pMins3d, radius=0.01, col=rgb(x$scene$colormap(x$geom$pMins3d[,3]), maxColorValue=256))
  spheres3d(x$geom$pMaxs3d, radius=0.01, col=rgb(x$scene$colormap(x$geom$pMaxs3d[,3]), maxColorValue=256))

  ## Draw the cylinders
  cIds = c()
  for(i in 1:length(x$scene$cylinders))
    cIds[i] <- shade3d(addNormals(x$scene$cylinders[[i]]), col=x$scene$cylColors1[[i]])  
  x$cIds <- cIds

  ## Draw the optional standard deviation tubes
  if(drawStdDev) {
    for(i in 1:length(x$scene$tubes)) {
      shade3d(addNormals(x$scene$tubes[[i]]), col="white", alpha=0.25, back="fill", front="cull")
      shade3d(addNormals(x$scene$tubes[[i]]), col="white", alpha=0.25, back="cull", front="fill")
    }
  }

  ## Create the orientation axes
  if(axesOn) {
    text3d(1, 0, 0, texts="x1")
    text3d(0, 1, 0, texts="x2")
    text3d(0, 0, 1, texts="f(x)")
    lines3d(c(1, 0, 0),c(0, 0, 1), c(0,0, 0), color="black")
    segments3d(c(0, 0),c(0, 0), c(0, 1), color="red")
  }
  
  ## Create the mouse selection mode
  x=mouseSelect(1, x)

  ## Return the object
  invisible(x)
}


## @-@-@-@-@-@-@ Plot Curves @-@-@-@-@-@-@-@-@-@ ##
plotCurves <- function(obj) {
  
  ## Calculate the number of pages, columns and rows
  if(length(obj$plotList) == 0)
    obj$plotList <- c( 1:ncol(obj$geom$curvesXHD[[1]]))
  nc <- length(obj$plotList)
  numPages=nc%/%25
  if(nc%%25 != 0)
    numPages = numPages+1
  col=(nc/numPages)%/%5+1
  if(col > 1)
    row=5
  else
    row=nc%%5

  ## Create the number of pages
  if(length(dev.list()) == 0) {
    for(d in 1:numPages)
      dev.new()
  }
  for(c in 1:nc) {
    if(c%%25 == 1) {
      dev.set(which=c%/%25+2)
      par(bg = "white", mfcol=c(row, col), mar=c(5, 4, 1, 1))
    }
    cindex <- obj$plotList[c]
    
    title <- paste("Plot #", cindex, ": ", obj$geom$colNames[cindex], " vs ", "Y")
    yLabs <- c(sprintf("%.2f", obj$geom$ylim[1]), sprintf("%.2f", (obj$geom$ylim[2]-obj$geom$ylim[1])/2),
               sprintf("%.2f", obj$geom$ylim[2]))
    plot(y=NA, x=NA, xlim=obj$geom$xlim[cindex,], ylim=range(0:1), 
         main=title, xlab=obj$geom$colNames[cindex], ylab="Y", yaxt="n", xaxt="n")
    axis(1, labels=c(0.0,0.5,1.0), at=c(0.0,0.5,1.0))
    axis(2, labels=yLabs, at=c(0.0,0.5,1.0))      
    
    for(i in 1:length(obj$scene$cylinders)){
      if(obj$scene$state[[i]] == 2){
        lightColor=hex(mixcolor(0.25, hex2RGB(obj$scene$cylColors2[[i]]), RGB(1,1,1)))
        lines(x=obj$geom$curvesXHD[[i]][,cindex], y=obj$geom$crvs3d[[i]][,3], 
              col=obj$scene$cylColors2[[i]], lwd=3)
        lines(x=obj$geom$curvesXHD[[i]][,cindex] +
            obj$geom$residuals[[i]][,cindex],
              y=obj$geom$crvs3d[[i]][,3], col=lightColor, lwd=1) 
        lines(x=obj$geom$curvesXHD[[i]][,cindex] -
            obj$geom$residuals[[i]][,cindex],
              y=obj$geom$crvs3d[[i]][,3], col=lightColor, lwd=1) 
      }
    }
  }
}

## - - - - - - - MOUSE SELECT - - - - - ##
mouseSelect <- function(button, obj) {
  
  

    ## (-)-(-)-(-)-(-) MOUSE BEGIN (-)-(-)-(-)-(-) ##
    mouseBegin <- function(x, y) {
      
      ##-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-##
      ## Define the area of the window that is selected,
      ## return a function that tests if the cylinders is selected
      createSelectFunction <- function(xa, ya, w=3, h=3) {
        proj <- rgl.projection()
        vp <- proj$view
        llx <- (xa - vp[1]) / vp[3]
        lly <- 1-(ya - vp[2]) / vp[4]
        urx <- (xa -vp[1] + w) / vp[3]
        ury <- 1-(ya - vp[2] +h) / vp[4]
        
        if ( llx > urx ){
          temp <- llx
          llx <- urx
          urx <- temp
        }
        if ( lly > ury ){
          temp <- lly
          lly <- ury
          ury <- temp
        }
        
        function(x,y=NULL,z=NULL) {
          pixel <- rgl.user2window(x,y,z,projection=proj)
          x <- pixel[,1]
          y <- pixel[,2]
          z <- pixel[,3]
          (llx <= x) & (x <= urx) & (lly <= y) & (y <= ury) & 
          (0 <= z) & (z <= 1)
        }
      }
      
      ## Create the picker function
      picker <- createSelectFunction(x, y)
      selected <- c(0,1)
      
      ## Cycle through all curves
      for(i in 1:length(obj$scene$cylinders)){
        
        ## Test to see if this curve was selected
        selector <- picker(t(obj$scene$cylinders[[i]]$vb[1:3,]))

        ## Check each point. If picked, save smallest Z
        for(j in 1:length(selector)){
          if(selector[j] == TRUE){
            thisZ <- obj$scene$cylinders[[i]]$vb[3,j]
            if(thisZ < selected[2])
              selected <- c(i, thisZ)
          }
        }
      }
      
      ## If we have selected more than is supported, throw an error
      if(selected[1] != 0){
        noSelect <- FALSE
        if(length(obj$selColors)==1){
          if(obj$scene$state[[selected[1]]] == 1)
            {
              noSelect <- TRUE
              print("ERROR: Too Many Selected Curves!! Unselect a curve.")
            }
        }

      ## If in selected state, pop the selected object
      if(!noSelect){
        rgl.pop(id=obj$cIds[[selected[1]]])
      }

        ## Re-draw the tubes if selected
        addColor <- FALSE
        takeColor <- FALSE
        
        if(obj$scene$state[[selected[1]]] == 1 && !noSelect){
          obj$scene$state[[selected[1]]]    <<- 2
          obj$scene$cylColors2[selected[1]] <<- obj$scene$selectColors[1]
          obj$cIds[[selected[1]]]           <<- shade3d(addNormals(obj$scene$cylinders[[selected[1]]]),
                                                        col=obj$scene$cylColors2[[selected[1]]])
          takeColor <- TRUE
        }
        else if(obj$scene$state[[selected[1]]] == 2)
          {
            obj$scene$state[[selected[1]]] <<- 1
            obj$cIds[[selected[1]]]        <<- shade3d(addNormals(obj$scene$cylinders[[selected[1]]]),
                                                            col=obj$scene$cylColors1[[selected[1]]])
            addColor <- TRUE
          }

        ## Add or take colors from the color vector  
        if(addColor)
          obj$scene$selectColors <<- c(obj$scene$cylColors2[[selected[1]]], obj$scene$selectColors)
        else if(takeColor)
          obj$scene$selectColors <<- obj$scene$selectColors[2:length(obj$scene$selectColors)]
        
        ## Plot the curves
        plotCurves(obj)
    
      }
      invisible(obj)
    }
    
    ## Set the rgl mouse call backs
    rgl.setMouseCallbacks(button, mouseBegin)
    invisible(obj)
  }
  

##-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-##
## Generate the vis pieces
##-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-##
scene.mscPlot <- function(geom, cSides=8, colorMap=0) {
    ## Create a list of selected colors
    selectColors <- c()  
    brewerSets=c("Set1", "Dark2", "Set3")
    for(b in 1:length(brewerSets)){ 
      cols <- brewer.pal(brewer.pal.info[brewerSets[b],]$maxcolors, brewerSets[b])
      selectColors <- c(selectColors, cols)
    }
    selectColors <- c(selectColors, "#000000")
    
    ## The blue, green, red colormap
    if(colorMap == 0)
      colormap = colorRamp( c('#0066CC', '#CCCC00', '#D22905') , interpolate="linear", bias=0.5)
    else if(colorMap == 1)
      colormap = colorRamp(c('#0520B0', '#f7f7f7', '#CA0020'), interpolate="linear", bias=0.5)
    else if(colorMap == 2)
      colormap = colorRamp(c('#7b3294', '#f7f7f7', '#008837'), interpolate="linear", bias=0.5)
    
    ## Create the cylinders
    cylinders <- c()
    cylColors1 <- c()
    cylColors2 <- c()
    state <- c()
    tubes <- c()  
    for(i in 1:geom$npartition){
      color <- rgb(colormap(geom$crvs3d[[i]][,3]), maxColorValue=255)
      fullColors <- c()
      for(p in 1:length(geom$crvs3d[[i]][,3]))
        fullColors[((p-1)*(cSides*4)+1) : (p*(cSides*4))] <- color[p]
      
      ## Create the tubes
      cylColors1[[i]] <- fullColors
      cylColors2[[i]] <- "black"
      cylinders[[i]] <- cylinder3d(geom$crvs3d[[i]], twist = 0, sides = cSides,
                                   e2=rbind(c(1,0,0),c(1,0,0)), radius = 0.008)
      tubes[[i]]     <- cylinder3d(geom$crvs3d[[i]], twist=0, sides=cSides,
                                   e2=rbind(c(1,0,0),c(1,0,0)), radius = geom$totResiduals[[i]])
      state[[i]] <- 1
    }    
    scene <- structure(list(geom=geom,cylinders=cylinders, cylColors1=cylColors1, cylColors2=cylColors2,
                            state=state, tubes=tubes, colormap=colormap, selectColors=selectColors))
    invisible(scene)
  }

##-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-##
## Generate the geometry for the scene
##-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-##
geometry.mscPlot <- function(ms, span = 0.35, nsamples = 50)
{
  msLevel <- ms$level[[ms$predictLevel]]
 
  ## Unique min and max indicies into ms$x/ms$y (domain/function value)
  uMins <- unique(msLevel$mins)
  uMaxs <- unique(msLevel$maxs)
  nmins <- length(uMins)
  nmaxs <- length(uMaxs)
  
  ## Do PCA on extrema
  extrema <- ms$x[c(uMins, uMaxs), ]
  pcaE <- prcomp(extrema)
  
  ## Location of mins and maxs
  pMins <- pcaE$x[1:nmins, 1:2]
  pMaxs <- pcaE$x[(nmins+1):(nmins+nmaxs), 1:2]

  # Scale the function values to be within the bounding box (0-1)
  miny <- min(ms$y)
  maxy <- max(ms$y)
 # ms$y <- (ms$y - miny) / (maxy-miny)
  
  ## Function values of mins and maxs
  fMins <- ms$y[uMins]
  fMaxs <- ms$y[uMaxs]
  
  ## The number of partition
  npartition <- length(msLevel$partitionSize)
  
  ## For each crystal store curve points (curvesX) and associated function value (curvesY)
  ## Curves are of the from x_i = c_i(y), e.g. a n-dimensional curve
  ## curve projected into 2D PCA sapce
  crvs3d <- c()
  ## curves in original space
  curvesXHD <- c()
  
  ## residuals for each curve, for each dimension 
  residuals <- c()
  for(i in 1:npartition){
    crvs3d[[i]] <- matrix(nrow=nsamples, ncol=3)

    residuals[[i]] <- matrix(nrow=nsamples, ncol=ncol(ms$x))
    eMax <- msLevel$maxs[i]
    eMin <- msLevel$mins[i]
    miny <- ms$y[eMin]
    maxy <- ms$y[eMax]
    h <- (maxy-miny) / (nsamples-1)
    sharedMax <- unique( which(msLevel$maxs == eMax) )
    sharedMin <- unique( which(msLevel$mins == eMin) )
    y <- ms$y[msLevel$partition == i]
    x <- ms$x[msLevel$partition == i, ]
    for(k in sharedMin){
      if(k != i){
        y2 <- ms$y[msLevel$partition == k]
        y2 <-  miny - (y2 - miny)
        y <- c(y, y2)
        x <- rbind(x, ms$x[msLevel$partition == k, ]) 
      }
    }
    for(k in sharedMax){
      if(k != i){
        y2 <- ms$y[msLevel$partition == k]
        y2 <-  maxy + (maxy-y2)
        y <- c(y, y2)
        x <- rbind(x, ms$x[msLevel$partition == k, ] )
      }
    }
    
    ye <- miny + (0:(nsamples-1))*h
    ye[nsamples] <- maxy
    xp <- matrix(nrow=nsamples, ncol=ncol(ms$x))
    for(k in 1:ncol(xp)){
      l <- loess(x[,k] ~ y, span=span, normalize=FALSE, degree=1)
      rl <- loess(I(abs(l$residuals)) ~ y, span=span/2.5, normalize=FALSE, degree=1)
      residuals[[i]][,k] <- predict(rl, ye)
      xp[, k] <- predict(l, ye)
    }
    curvesXHD[[i]] <- xp
    pcaC <- prcomp(xp)
    
    iMin <- which(uMins == eMin) 
    iMax <- which(uMaxs == eMax)
    xp  <- pcaC$x[, 1:2]
    if(nmins == 1)
      anchor <- pMins      
    else
      anchor <- pMins[iMin, ]       
    
    xp[, 1] <- xp[, 1] - xp[1, 1] + anchor[1]
    xp[, 2] <- xp[, 2] - xp[1, 2] + anchor[2]
    if(nmaxs == 1){
      stretch <- pMaxs - xp[nsamples, ]       
    }
    else
      stretch <- pMaxs[iMax, ] - xp[nsamples, ]       
    
    for(k in 1:nsamples)
      xp[k, ] <- xp[k, ] + (k-1)/(nsamples-1) * stretch
    
    crvs3d[[i]][, 1:2] <- xp 
    crvs3d[[i]][,3] <- ye
  }


  ## Extend the points & curves to 3D
  if(nmins ==1){
    pMins3d <- t(as.matrix(c(pMins, fMins)))
  }
  else{
    pMins3d <- cbind(pMins, fMins)
  }
  if(nmaxs == 1){
    pMaxs3d <- t(as.matrix(c(pMaxs, fMaxs)))
  }
  else{
    pMaxs3d <- cbind(pMaxs, fMaxs)
  }
  pTmins <- c(min(c(pMins3d[,1], pMaxs3d[,1])), min(c(pMins3d[,2], pMaxs3d[,1])), min(ms$y))
  pTmaxs <- c(max(c(pMins3d[,1], pMaxs3d[,1])), max(c(pMins3d[,2], pMaxs3d[,2])), max(ms$y))

  for(i in 1:npartition){    
    crvs3d[[i]][,1] <- (crvs3d[[i]][,1] - pTmins[1])/(pTmaxs[1] - pTmins[1])
    crvs3d[[i]][,2] <- (crvs3d[[i]][,2] - pTmins[2])/(pTmaxs[2] - pTmins[2])
    crvs3d[[i]][,3] <- (crvs3d[[i]][,3] - pTmins[3])/(pTmaxs[3] - pTmins[3])
  }
  
  for(i in 1:nmins)
    pMins3d[i,] <- (pMins3d[i,] - pTmins)/ (pTmaxs - pTmins)
  
  for(i in 1:nmaxs)
    pMaxs3d[i,] <- (pMaxs3d[i,] - pTmins)/ (pTmaxs - pTmins)
  

  ## Get the complete residual and scale appropriatly
  totResiduals <- c()
  
  totScale <- 1 / mean(dist(extrema))
  for(i in 1:npartition)
    totResiduals[[i]] <- sqrt(rowSums(residuals[[i]]^2))*totScale


  # Set the limits of the x
  xlim <- matrix(nrow=ncol(ms$x), ncol=2)
  for( i in 1:ncol(ms$x))
    xlim[i,] <- range(ms$x[, i])
  
  ## Return lots of info
  mscGeom <- structure(list(npartition=npartition,
                            curvesXHD = curvesXHD, 
                            crvs3d=crvs3d,
                            residuals = residuals,
                            totResiduals=totResiduals,
                            pMins3d=pMins3d,
                            pMaxs3d=pMaxs3d, 
                            ylim= c(pTmins[3], pTmaxs[3]),
                            xlim = xlim,
                            colNames = colnames(ms$x)))
}


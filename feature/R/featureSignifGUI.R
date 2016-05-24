featureSignifGUI <- function(x, scaleData=FALSE)
{
  fscreate.tcl <- function()
  {
    xlim <- as.numeric(c(tclvalue(xlim1.tcl), tclvalue(xlim2.tcl)))
    bw <- as.numeric(tclvalue(bw1.tcl))
    if (d>=2)
    {
      ylim <- as.numeric(c(tclvalue(ylim1.tcl), tclvalue(ylim2.tcl)))
      bw <- as.numeric(c(tclvalue(bw1.tcl), tclvalue(bw2.tcl)))
    }
    if (d>=3)
    {
      zlim <- as.numeric(c(tclvalue(zlim1.tcl), tclvalue(zlim2.tcl)))
      bw <- as.numeric(c(tclvalue(bw1.tcl), tclvalue(bw2.tcl), tclvalue(bw3.tcl)))
    }
    gs <- rep(as.numeric(tclvalue(gridsize.tcl)), d)
    fs <- featureSignif(x, bw=bw, addSignifGrad=TRUE, addSignifCurv=TRUE, gridsize=gs)
    assign("fs", fs, envir=fs.env)

    if (d==1) plot(fs, addKDE=TRUE, add=FALSE, xlim=xlim, xlab=tclvalue(xlab.tcl))
    if (d==2) plot(fs, addKDE=TRUE, add=FALSE, xlim=xlim, ylim=ylim, xlab=tclvalue(xlab.tcl), ylab=tclvalue(ylab.tcl))
    if (d==3) plot(fs, addKDE=TRUE, add=FALSE, xlim=xlim, ylim=ylim, zlim=zlim, xlab=tclvalue(xlab.tcl), ylab=tclvalue(ylab.tcl), zlab=tclvalue(zlab.tcl), addAxes3d=as.logical(as.numeric(tclvalue(addAxes3d.tcl))))
    tkmessageBox(title="featureSignifGUI", message="Feature significance computations complete.", type="ok")
    return()
  }
  
  fssizer.tcl <- function()
  {
    xlim <- as.numeric(c(tclvalue(xlim1.tcl), tclvalue(xlim2.tcl)))
    xlab <- tclvalue(xlab.tcl)
    gs <- rep(as.numeric(tclvalue(gridsize.tcl)), d)
    fs.env$fs.SiZer <- SiZer(as.vector(x), bw=bw.range[[1]], gridsize=gs, xlim=xlim, xlab=xlab, plotSiZer=TRUE)
    tkmessageBox(title="featureSignifGUI", message="SiZer map computed.", type="ok")
    return()
  }

  fssicon.tcl <- function()
  {
    xlim <- as.numeric(c(tclvalue(xlim1.tcl), tclvalue(xlim2.tcl)))
    xlab <- tclvalue(xlab.tcl)
    gs <- rep(as.numeric(tclvalue(gridsize.tcl)), d)
    fs.env$fs.SiCon <- SiCon(as.vector(x), bw=bw.range[[1]], gridsize=gs, xlim=xlim, xlab=xlab, plotSiCon=TRUE)
    tkmessageBox(title="featureSignifGUI", message="SiCon map computed.", type="ok")
    return()
  }
  
  fsdata.tcl <- function()
  {
    addDataNum <- as.numeric(tclvalue(addDataNum.tcl))
    if (d>=3) dataAlpha <- as.numeric(tclvalue(dataAlpha.tcl))
    if (!is.null(fs.env$fs))
      if (d<3)
        plot(fs.env$fs, addKDE=FALSE, addData=TRUE, add=TRUE, addDataNum=addDataNum)
      else
        plot(fs.env$fs, addKDE=FALSE, addData=TRUE, add=TRUE, dataAlpha=dataAlpha, addDataNum=addDataNum)
    return()
  }

  fsgraddata.tcl <- function()
  {
    addDataNum <- as.numeric(tclvalue(addDataNum.tcl))
    if (d>=3) gradDataAlpha <- as.numeric(tclvalue(gradDataAlpha.tcl))
    if (!is.null(fs.env$fs))
      if (d<3)
        plot(fs.env$fs, addKDE=FALSE, addSignifGradData=TRUE, add=TRUE, addDataNum=addDataNum)
      else
        plot(fs.env$fs, addKDE=FALSE, addSignifGradData=TRUE, add=TRUE, gradDataAlpha=gradDataAlpha, addDataNum=addDataNum)
    return()
  }

  fsgradregion.tcl <- function(panel)
  {
    if (d>=3) gradRegionAlpha <- as.numeric(tclvalue(gradRegionAlpha.tcl))
    if (!is.null(fs.env$fs))
      if (d<3)
        plot(fs.env$fs, addKDE=FALSE, addSignifGradRegion=TRUE, add=TRUE)
      else
        plot(fs.env$fs, addKDE=FALSE, addSignifGradRegion=TRUE, add=TRUE, gradRegionAlpha=gradRegionAlpha, addAxes=as.logical(as.numeric(tclvalue(addAxes3d.tcl))))
    return()
  }

  fscurvdata.tcl <- function(panel)
  {
    addDataNum <- as.numeric(tclvalue(addDataNum.tcl))
    if (d>=3) curvDataAlpha <- as.numeric(tclvalue(curvDataAlpha.tcl))
    if (!is.null(fs.env$fs))
      if (d <3)
        plot(fs.env$fs, addKDE=FALSE, addSignifCurvData=TRUE, add=TRUE, addDataNum=addDataNum)
      else
        plot(fs.env$fs, addKDE=FALSE, addSignifCurvData=TRUE, add=TRUE, curvDataAlpha=curvDataAlpha, addDataNum=addDataNum)
    return()
  }

  fscurvregion.tcl <- function(panel)
  {
    if (d>=3) curvRegionAlpha <- as.numeric(tclvalue(curvRegionAlpha.tcl))
    if (!is.null(fs.env$fs))
      if (d <3)
         plot(fs.env$fs, addKDE=FALSE, addSignifCurvRegion=TRUE, add=TRUE)
      else
        plot(fs.env$fs, addKDE=FALSE, addSignifCurvRegion=TRUE, add=TRUE, curvRegionAlpha=curvRegionAlpha, addAxes=as.logical(as.numeric(tclvalue(addAxes3d.tcl))))
    return(panel)
  }

  fsclearall1.tcl <- function()
  {
    xlim <- as.numeric(c(tclvalue(xlim1.tcl), tclvalue(xlim2.tcl)))
    if (d>=2) ylim <- as.numeric(c(tclvalue(ylim1.tcl), tclvalue(ylim2.tcl)))
    if (d>=3) zlim <- as.numeric(c(tclvalue(zlim1.tcl), tclvalue(zlim2.tcl)))
           
    if (!is.null(fs.env$fs))
    {
      if (d==1) plot(fs.env$fs, addKDE=TRUE, add=FALSE, xlim=xlim, xlab=tclvalue(xlab.tcl))
      if (d==2) plot(fs.env$fs, addKDE=TRUE, add=FALSE, xlim=xlim, ylim=ylim, xlab=tclvalue(xlab.tcl), ylab=tclvalue(ylab.tcl))
      if (d==3) plot(fs.env$fs, addKDE=TRUE, add=FALSE, xlim=xlim, ylim=ylim, zlim=zlim, xlab=tclvalue(xlab.tcl), ylab=tclvalue(ylab.tcl), zlab=tclvalue(zlab.tcl),  addAxes3d=as.logical(as.numeric(tclvalue(addAxes3d.tcl))))
    }
    return()
  }
  fsclearall2.tcl <- function(panel)
  {
    xlim <- as.numeric(c(tclvalue(xlim1.tcl), tclvalue(xlim2.tcl)))
    if (d>=2) ylim <- as.numeric(c(tclvalue(ylim1.tcl), tclvalue(ylim2.tcl)))
    if (d>=3) zlim <- as.numeric(c(tclvalue(zlim1.tcl), tclvalue(zlim2.tcl)))
           
    if (!is.null(fs.env$fs))
    {
      if (d==1) plot(fs.env$fs, addKDE=FALSE, add=FALSE, xlim=xlim, xlab=tclvalue(xlab.tcl))
      if (d==2) plot(fs.env$fs, addKDE=FALSE, add=FALSE, xlim=xlim, ylim=ylim, xlab=tclvalue(xlab.tcl), ylab=tclvalue(ylab.tcl))
      if (d==3) plot(fs.env$fs, addKDE=FALSE, add=FALSE, xlim=xlim, ylim=ylim, zlim=zlim, xlab=tclvalue(xlab.tcl), ylab=tclvalue(ylab.tcl), zlab=tclvalue(zlab.tcl),  addAxes3d=as.logical(as.numeric(tclvalue(addAxes3d.tcl))))
    }
    return()
  }
  
  ############################################################################
  ## GUI
  ############################################################################

  fs.env <- new.env()
    
  button.col <- "lightblue"
  bg.col <- "white"
  ##fg.col <- "black"
  heading.col <- "blue"
  heading.font <- tkfont.create(family="helvetica", weight="bold")
  space.font <- tkfont.create(size=10)
  textwidth <- 10
  sf <- 4
  #space.label <- paste(rep(" ", scalebar.width), collapse="")
  tcl("tk_setPalette", "background", bg.col, "selectColor", "grey75") 

  ## set defaults
  
  if (is.vector(x))
  { d <- 1;  names.x <- deparse(substitute(x)); x <- as.matrix(x); if (scaleData)  x <- (x-min(x))/(max(x) - min(x))}
  else
  { d <- ncol(x); names.x <- colnames(x)
    if (is.null(names.x))
    {
       names.xx <- deparse(substitute(x))
       names.xx <- strsplit(names.xx, "\\[")[[1]][1]
       names.x <- paste(names.xx, "[,", 1:d, "]", sep="")
    }
    if (scaleData)
      for (i in 1:d)
        x[,i] <- (x[,i]-min(x[,i]))/(max(x[,i]) - min(x[,i]))
  }
  if (d>4)
    stop("Feature significance only available for 1- to 4-d data")
  
  tau <- 5

  if (d==1) gridsize <- 401
  if (d==2) gridsize <- 151
  if (d==3) gridsize <- 31
  if (d==4) gridsize <- 21
  bw.range <- dfltBWrange(x, tau)
  bw <- matrix(unlist(bw.range), nrow=2, byrow=FALSE)
  h.low <- bw[1,]
  h.upp <- bw[2,]
  hmix.prop <- 1/4
  h <- h.low^(hmix.prop)*h.upp^(1-hmix.prop)

  if (d==1)
  {  
    xlim <- c(min(x)-h[1],max(x)+h[1])
    dfltCounts.out <- dfltCounts(x,gridsize, apply(bw, 2, max))
    gcounts <- dfltCounts.out$counts
    range.x <- dfltCounts.out$range.x  
    dest <- drvkde(gcounts, rep(0,d), bandwidth=h, binned=TRUE, range.x=range.x, se=FALSE)
    ylim <- c(0,1.5)*max(dest$est)
  }
  else
  {
    xlim <- c(min(x[,1])-h[1],max(x[,1])+h[1])
    ylim <- c(min(x[,2])-h[2],max(x[,2])+h[2])
    if (d>2)
      zlim <- c(min(x[,3])-h[3],max(x[,3])+h[3])
  }

  ## panel dimensions
  if (d <3)
  {  
    fspanel.width <- 740
    fspanel.height <- 410
    fspanel.button.width <- (fspanel.width-40)/3
  }
  else
  {
    fspanel.width <- 950
    fspanel.height <- 410
    fspanel.button.width <- (fspanel.width-50)/4
  }
  
  tt <- tktoplevel(width=fspanel.width, height=fspanel.height)
  tktitle(tt) <- paste("featureSignif GUI")
  tkgrid(tklabel(tt,text="featureSignif GUI: feature significance for multivariate kernel density estimation", font=heading.font, foreground=heading.col), row=1, columnspan=7+2*(d>=3))
  tkgrid(tklabel(tt,text=" "), row=2, column=1)
  tkgrid(tklabel(tt,text=" "), row=2, column=3)
  tkgrid(tklabel(tt,text=" "), row=2, column=5)
  tkgrid(tklabel(tt,text=" "), row=2, column=7)
  if (d==3) tkgrid(tklabel(tt,text=" "), row=2, column=9)
  tkgrid(tklabel(tt,text=" "), row=11, columnspan=7+2*(d>=3))
  
  ## Bandwidths
  
  col1.frame <- tkframe(tt)
  
  bw1.tcl <- tclVar(h[1])
  bw2.tcl <- tclVar(h[2])
  bw3.tcl <- tclVar(h[3])

  if (d>=1)
  {
    col1b.frame <- tkframe(col1.frame)
    bw1.s <- tkscale(col1b.frame, from=h.low[1], to=h.upp[1], showvalue=TRUE, resolution=(h.upp[1]-h.low[1])/100, length=fspanel.button.width, variable=bw1.tcl, orient="horizontal", label=paste("Bandwidth for", names.x[1]))
    tkconfigure(bw1.s, variable=bw1.tcl) 
    tkgrid(bw1.s, sticky="we")
    tkgrid(col1b.frame)
    tkgrid(tklabel(col1.frame, text=" "), sticky="w")
    tkgrid(col1.frame, row=3, column=2, sticky="nw")
  }
  if (d>=2)
  {
    col1b.frame <- tkframe(col1.frame)
    bw2.s <- tkscale(col1.frame, from=h.low[2], to=h.upp[2], showvalue=TRUE,  resolution=(h.upp[2]-h.low[2])/100, length=fspanel.button.width, variable=bw2.tcl, orient="horizontal", label=paste("Bandwidth for", names.x[2]))
    tkconfigure(bw2.s, variable=bw2.tcl) 
    tkgrid(bw2.s, sticky="we")
    tkgrid(col1b.frame)
    tkgrid(tklabel(col1.frame, text=" "), sticky="w")
    tkgrid(col1.frame, row=5, column=2, sticky="nw")
  }
  if (d>=3)
  {
    bw3.s <- tkscale(col1.frame, from=h.low[3], to=h.upp[3], showvalue=TRUE,  resolution=(h.upp[3]-h.low[3])/100, length=fspanel.button.width, variable=bw3.tcl, orient="horizontal", label=paste("Bandwidth for", names.x[3]))
    tkconfigure(bw3.s, variable=bw3.tcl) 
    tkgrid(bw3.s, sticky="we")
    tkgrid(col1b.frame)
    tkgrid(col1.frame, row=7, column=2, sticky="we")
  }
  tkgrid(tklabel(col1.frame, text=" "), sticky="w")
  
  ## Grid size

  col1b.frame <- tkframe(col1.frame)
  gridsize.tcl <- tclVar(gridsize)
  grid.e <- tkentry(col1b.frame, textvariable=gridsize.tcl, width=textwidth/2)
  tkgrid(tklabel(col1b.frame, text="Grid size "), grid.e, sticky="nw")
  tkgrid(col1b.frame, sticky="w")
  tkgrid(col1.frame, sticky="nw")


  ## Axis limits
  col2.frame <- tkframe(tt)

  if (d>=1)
  { 
    col2b.frame <- tkframe(col2.frame)
    xlab.tcl <- tclVar(names.x[1])
    xlab.e <- tkentry(col2b.frame, textvariable=xlab.tcl)
    tkgrid(tklabel(col2b.frame, text="Axis label", justify="left"), xlab.e, sticky="nw")
    tkgrid(col2b.frame, sticky="w")

    col2b.frame <- tkframe(col2.frame)
    xlim1.tcl <- tclVar(signif(xlim[1], digits=sf))
    xlim1.e <- tkentry(col2b.frame, textvariable=xlim1.tcl, width=textwidth)
    tkgrid(tklabel(col2b.frame, text="Axis lower limit ", justify="left"), xlim1.e, sticky="nw")
    tkgrid(col2b.frame, sticky="w")
    
    col2b.frame <- tkframe(col2.frame)
    xlim2.tcl <- tclVar(signif(xlim[2], digits=sf))
    xlim2.e <- tkentry(col2b.frame, textvariable=xlim2.tcl, width=textwidth)
    tkgrid(tklabel(col2b.frame, text="Axis upper limit", justify="left"), xlim2.e, sticky="nw")
    tkgrid(col2b.frame, sticky="w")
    
    tkgrid(tklabel(col2.frame, text=" "), sticky="w")
    tkgrid(col2.frame, row=3, column=4, sticky="ne")
    
  }
  if (d>=2)
  {  
    col2b.frame <- tkframe(col2.frame)
    ylab.tcl <- tclVar(names.x[2])
    ylab.e <- tkentry(col2b.frame, textvariable=ylab.tcl)
    tkgrid(tklabel(col2b.frame, text="Axis label "), ylab.e)
    tkgrid(col2b.frame, sticky="w")

    col2b.frame <- tkframe(col2.frame)
    ylim1.tcl <- tclVar(signif(ylim[1], digits=sf))
    ylim1.e <- tkentry(col2b.frame, textvariable=ylim1.tcl, width=textwidth)
    tkgrid(tklabel(col2b.frame, text="Axis lower limit "), ylim1.e)
    tkgrid(col2b.frame, sticky="w")
    
    col2b.frame <- tkframe(col2.frame)
    ylim2.tcl <- tclVar(signif(ylim[2], digits=sf))
    ylim2.e <- tkentry(col2b.frame, textvariable=ylim2.tcl, width=textwidth)
    tkgrid(tklabel(col2b.frame, text="Axis upper limit"), ylim2.e)
    tkgrid(col2b.frame, sticky="w")
   
    tkgrid(tklabel(col2.frame, text=" "), sticky="w")
    tkgrid(col2.frame, row=5, column=4, sticky="ne")
  }
  if (d>=3)
  {  
    col2b.frame <- tkframe(col2.frame)
    zlab.tcl <- tclVar(names.x[3])
    zlab.e <- tkentry(col2b.frame, textvariable=zlab.tcl)
    tkgrid(tklabel(col2b.frame, text="Axis label "), zlab.e)
    tkgrid(col2b.frame, sticky="w")
    
    col2b.frame <- tkframe(col2.frame)
    zlim1.tcl <- tclVar(signif(zlim[1], digits=sf))
    zlim1.e <- tkentry(col2b.frame, textvariable=zlim1.tcl, width=textwidth)
    tkgrid(tklabel(col2b.frame, text="Axis lower limit "), zlim1.e)
    tkgrid(col2b.frame, sticky="w")
    
    col2b.frame <- tkframe(col2.frame)
    zlim2.tcl <- tclVar(signif(zlim[2], digits=sf))
    zlim2.e <- tkentry(col2b.frame, textvariable=zlim2.tcl, width=textwidth)
    tkgrid(tklabel(col2b.frame, text="Axis upper limit"), zlim2.e)
    tkgrid(col2b.frame, sticky="w")
    tkgrid(col2.frame, row=7, column=4, sticky="w")
    tkgrid(tklabel(col2.frame, text=" "), sticky="w")
  }
  
  col2b.frame <- tkframe(col2.frame)
  addDataNum.tcl <- tclVar(1000)
  addDataNum.e <- tkentry(col2b.frame, textvariable=addDataNum.tcl, width=textwidth/2)
  tkgrid(tklabel(col2b.frame, text="No. of display data points "), addDataNum.e)
  tkgrid(col2b.frame, sticky="w")
  tkgrid(col2.frame,  sticky="nw")

  ## Plot options
  
  col3.frame <- tkframe(tt)
  col3b.frame <- tkframe(col3.frame)
  fsdata.b <- tkbutton(col3b.frame, text=" ", command=fsdata.tcl, background=button.col)
  fsgradregion.b <- tkbutton(col3b.frame, text=" ", command=fsgradregion.tcl, background=button.col)
  fsgraddata.b <- tkbutton(col3b.frame, text=" ", command=fsgraddata.tcl, background=button.col)
  fscurvregion.b <- tkbutton(col3b.frame, text=" ", command=fscurvregion.tcl, background=button.col)
  fscurvdata.b <- tkbutton(col3b.frame, text=" ", command=fscurvdata.tcl, background=button.col)
  tkgrid(fsdata.b, tklabel(col3b.frame, text="Add data\npoints", justify="left"), sticky="nw")
  tkgrid(fsgradregion.b, tklabel(col3b.frame, text="Add significant\ngradient region", justify="left"), sticky="nw")
  tkgrid(fsgraddata.b, tklabel(col3b.frame, text="Add significant\ngradient data", justify="left"), sticky="nw")
  tkgrid(fscurvregion.b, tklabel(col3b.frame, text="Add significant\ncurvature region", justify="left"), sticky="nw")
  tkgrid(fscurvdata.b, tklabel(col3b.frame, text="Add significant\ncurvature data", justify="left"), sticky="nw")
  tkgrid(col3b.frame)
  if (d==3) tkgrid(col3.frame, row=2*d+1, column=8, sticky="nw")
  else tkgrid(col3.frame, row=2*d+1, column=6, sticky="nw")

  ## Transparency options

  if (d==3)
  {
    col4.frame <- tkframe(tt)

    dataAlpha.tcl <- tclVar(0.1)
    dataAlpha.s <- tkscale(col4.frame, from=0, to=1, showvalue=TRUE,  resolution=0.1, length=fspanel.button.width, variable=dataAlpha.tcl, orient="horizontal", label="Transparency of data points")
    tkconfigure(dataAlpha.s, variable=dataAlpha.tcl) 
    tkgrid(dataAlpha.s, sticky="we")

    gradRegionAlpha.tcl <- tclVar(0.2)
    gradRegionAlpha.s <- tkscale(col4.frame, from=0, to=1, showvalue=TRUE,  resolution=0.1, length=fspanel.button.width, variable=dataAlpha.tcl, orient="horizontal", label="Transparency of signif. grad. regions")
    tkconfigure(gradRegionAlpha.s, variable=gradRegionAlpha.tcl) 
    tkgrid(gradRegionAlpha.s, sticky="we")

    gradDataAlpha.tcl <- tclVar(0.2)
    gradDataAlpha.s <- tkscale(col4.frame, from=0, to=1, showvalue=TRUE,  resolution=0.1,length=fspanel.button.width, variable=dataAlpha.tcl, orient="horizontal", label="Transparency of signif. grad. data")
    tkconfigure(gradDataAlpha.s, variable=gradDataAlpha.tcl) 
    tkgrid(gradDataAlpha.s, sticky="we")

    curvRegionAlpha.tcl <- tclVar(0.3)
    curvRegionAlpha.s <- tkscale(col4.frame, from=0, to=1, showvalue=TRUE,  resolution=0.1, length=fspanel.button.width, variable=dataAlpha.tcl, orient="horizontal", label="Transparency of signif. curv. regions")
    tkconfigure(curvRegionAlpha.s, variable=curvRegionAlpha.tcl) 
    tkgrid(curvRegionAlpha.s, sticky="we")

    curvDataAlpha.tcl <- tclVar(0.3)
    curvDataAlpha.s <- tkscale(col4.frame, from=0, to=1, showvalue=TRUE,  resolution=0.1, length=fspanel.button.width, variable=dataAlpha.tcl, orient="horizontal", label="Transparency of signif. curv. data")
    tkconfigure(curvDataAlpha.s, variable=curvDataAlpha.tcl) 
    tkgrid(curvDataAlpha.s, sticky="we")

    col4b.frame <- tkframe(col4.frame)
    addAxes3d.tcl <- tclVar(TRUE)
    addAxes3d.cb <- tkcheckbutton(col4b.frame, variable=addAxes3d.tcl)
    tkgrid(addAxes3d.cb, tklabel(col4b.frame, text="Add 3D axes "))
    tkgrid(col4b.frame, sticky="nw")
    tkgrid(col4.frame, row=7, column=6, sticky="nw")
  }

  ## Computation buttons
  
  col1.frame <- tkframe(tt)
  fscreate.b <- tkbutton(col1.frame, text=" ", command=fscreate.tcl, background=button.col)
  tkgrid(fscreate.b, tklabel(col1.frame, text="Compute significant\nfeatures", justify="left"), sticky="nw")
  tkgrid(col1.frame, row=9, column=2, sticky="nw")

  col2.frame <- tkframe(tt)
  fsclearall1.b <- tkbutton(col2.frame, text=" ", command=fsclearall1.tcl, background=button.col)
  tkgrid(fsclearall1.b, tklabel(col2.frame, text="Reset plot (except KDE)", justify="left"), sticky="nw")
  tkgrid(col2.frame, row=9, column=4, sticky="nw")

  col3.frame <- tkframe(tt)
  if (d==1)
  {
    fssizer.b <- tkbutton(col3.frame, text=" ", command=fssizer.tcl, background=button.col)
    fssicon.b <- tkbutton(col3.frame, text=" ", command=fssicon.tcl, background=button.col)
    tkgrid(fssizer.b, tklabel(col3.frame, text="Compute SiZer map", justify="left"), sticky="nw")
    tkgrid(fssicon.b, tklabel(col3.frame, text="Compute SiCon map", justify="left"), sticky="nw")
  }
  else
  {
    fsclearall2.b <- tkbutton(col3.frame, text=" ", command=fsclearall2.tcl, background=button.col)
    tkgrid(fsclearall2.b, tklabel(col3.frame, text="Reset plot", justify="left"), sticky="nw")
  }
  tkgrid(col3.frame, row=9, column=6, sticky="nw")
}



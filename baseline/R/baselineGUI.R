### $Id: baselineGUI.R 193 2012-06-24 21:13:42Z kristl $

## Baseline parameters, expandable list
baselineAlgorithmsGUI <- list()
baselineAlgorithmsGUI$irls					<- as.data.frame(matrix(c(0,10,0.1,5, 5,15,0.1,8, 0,0.5,0.01,0.05, 50,200,25,100), 4,4, byrow=TRUE))
dimnames(baselineAlgorithmsGUI$irls)		<- list(par=c("lambda1", "lambda2", "wi", "maxit"),val=c("min","max","step","default"))
baselineAlgorithmsGUI$irls$current		 	<- c(5,8,0.05,100)
baselineAlgorithmsGUI$irls$name				<- c("Primary smoothing", "Main smoothing", "Weighting", "Maximum iterations")
baselineAlgorithmsGUI$modpolyfit			<- as.data.frame(matrix(c(0,10,1,4, 0,15,1,4, 25,200,25,100), 3,4, byrow=TRUE))
dimnames(baselineAlgorithmsGUI$modpolyfit)	<- list(par=c("degree", "tol", "rep"),val=c("min","max","step","default"))
baselineAlgorithmsGUI$modpolyfit$current	<- c(4,4,100)
baselineAlgorithmsGUI$modpolyfit$name		<- c("Polynomial degree", "Update tolerance 10^-", "Max #iterations")
baselineAlgorithmsGUI$als					<- as.data.frame(matrix(c(0,15,1,6, 0,0.5,0.001,0.05), 2,4, byrow=TRUE))
dimnames(baselineAlgorithmsGUI$als)			<- list(par=c("lambda", "p"),val=c("min","max","step","default"))
baselineAlgorithmsGUI$als$current			<- c(6,0.05)
baselineAlgorithmsGUI$als$name				<- c("Smoothing parameter", "Residual weighting")
baselineAlgorithmsGUI$rollingBall 			<- as.data.frame(matrix(c(0,500,10,300, 0,500,10,200), 2,4, byrow=TRUE))
dimnames(baselineAlgorithmsGUI$rollingBall)	<- list(par=c("wm", "ws"),val=c("min","max","step","default"))
baselineAlgorithmsGUI$rollingBall$current	<- c(300,200)
baselineAlgorithmsGUI$rollingBall$name		<- c("Min/max window width", "Smoothing window width")
baselineAlgorithmsGUI$medianWindow			<- as.data.frame(matrix(c(0,500,10,300, 0,500,10,200), 2,4, byrow=TRUE))
dimnames(baselineAlgorithmsGUI$medianWindow)<- list(par=c("hwm", "hws"),val=c("min","max","step","default"))
baselineAlgorithmsGUI$medianWindow$current	<- c(300,200)
baselineAlgorithmsGUI$medianWindow$name		<- c("Median window half width", "Smoothing window half width")
baselineAlgorithmsGUI$fillPeaks		 		<- as.data.frame(matrix(c(0,15,1,6, 10,500,10,50, 1,100,1,10, 100,5000,100,2000), 4,4, byrow=TRUE))
dimnames(baselineAlgorithmsGUI$fillPeaks)	<- list(par=c("lambda", "hwi", "it", "int"),val=c("min","max","step","default"))
baselineAlgorithmsGUI$fillPeaks$current		<- c(6,50,10,2000)
baselineAlgorithmsGUI$fillPeaks$name		<- c("Primary smoothing", "Half width of local windows", "Maximum number of iterations", "Number of buckets")
baselineAlgorithmsGUI$peakDetection		 	<- as.data.frame(matrix(c(50,500,50,300, 10,200,10,50), 2,4, byrow=TRUE))
dimnames(baselineAlgorithmsGUI$peakDetection)<- list(par=c("left.right", "lwin.rwin"),val=c("min","max","step","default"))
baselineAlgorithmsGUI$peakDetection$current <- c(300,50)
baselineAlgorithmsGUI$peakDetection$name	<- c("Peak window width", "Smoothing window width")
baselineAlgorithmsGUI$rfbaseline		 	<- as.data.frame(matrix(c(100,5000,10,1000, 1,5,0.5,3.5), 2,4, byrow=TRUE))
dimnames(baselineAlgorithmsGUI$rfbaseline)	<- list(par=c("NoXP", "b"),val=c("min","max","step","default"))
baselineAlgorithmsGUI$rfbaseline$current	<- c(1000,3.5)
baselineAlgorithmsGUI$rfbaseline$name		<- c("Number of regression points", "Relative weighting")
baselineAlgorithmsGUI$shirley			 	<- as.data.frame(matrix(c(10,1000,10,50, 1e-7,1e-5,1e-7,1e-6), 2,4, byrow=TRUE))
dimnames(baselineAlgorithmsGUI$shirley)		<- list(par=c("maxit", "err"),val=c("min","max","step","default"))
baselineAlgorithmsGUI$shirley$current		<- c(50,1e-6)
baselineAlgorithmsGUI$shirley$name			<- c("Max #iterations", "Error")

baselineGUI <- function(spectra, method='irls', labels, rev.x=FALSE){
  if(requireNamespace("gWidgets", quietly = TRUE)){
    if(exists("baselineAlgorithmsGUI",envir=.GlobalEnv)){
      bAGUI <- get("baselineAlgorithmsGUI",envir=.GlobalEnv)
    } else {
      bAGUI <- baselineAlgorithmsGUI
    }
    
    
    ##
    ## Ititialise variables that are shared between the elements of the GUI:
    ##
    
    
    ## Spectrum parameters
    Y <- spectra; specNo <- 1
    n <- length(Y[1,]); n1 <- length(Y[1,])
    x <- 1:n1
    if(missing(labels)) # X-axis labels
      labels <- 1:n
    
    ## Plotting parameters
    xz <- 1; yz <- 1; xc <- 0; yc <-0
    setZoom <- function(){
      xz <<- 1; yz <<- 1; xc <<- 0; yc <<-0
    }
    gridOn <- FALSE           # Is the grid on?
    visibleZoom <- FALSE      # Has the zoom tools been activated?
    visibleCustomize <- FALSE # Has the parameter customization been activated?
    visibleExport <- FALSE	  # Has the export of all spectra been started?
    
    ##
    ## Define functions that are used by the GUI elements
    ##
    
    
    ## --=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--
    ## Baseline computation
    baseline.compute <- function(){
      clearPlot()
      ## Compute baseline based on current method and settings
      spec <- list()
      command <- "spec <- baseline(Y[specNo,,drop=FALSE]"
      for(i in 1:dim(bAGUI[[method]])[1]){
        command <- paste(command, ", ", rownames(bAGUI[[method]])[i], "=", bAGUI[[method]]$current[i], sep="")
      }
      command <- paste(command, ", method='", method, "')", sep="")
      eval(parse(text=command))
      ## Kludge to aviod warnings from R CMD check:
      # assign("baseline.result", spec, .GlobalEnv)
      putBaselineEnv("baseline.result", spec)
      putBaselineEnv("baseline.current", list(method=method, parNames=rownames(bAGUI[[method]]), parValues=bAGUI[[method]]$current))
      updatePlot()
    }                                   # end of baseline.compute
    
    ## --=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--
    ## Clear plot
    clearPlot <- function(){
      par(new = TRUE, mfrow = c(1,1))
      plot(0, 0, xlim=c(-1,1), ylim=c(-1,1), xlab="", ylab="", main="",
           axes=FALSE, col='white')
      par(new=FALSE)
      C <- as.list(par("usr")); names(C) <- c("xmin", "xmax", "ymin", "ymax")
      ##print(unix.time(
      rect(C$xmin, C$ymin, C$xmax, C$ymax, angle = 0, density = 20,
           col = 'white', lwd = 2)
      ##))
      text(0, 0, labels = "Calculating baseline...", col = 'blue', cex = 2)
    }
    
    ## --=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--
    ## Update plot
    updatePlot <- function() {
      ## FIXME: get is a kludge to avoid warnings from R CMD check:
      plot(getBaselineEnv("baseline.result"), grid = gridOn, labels = labels, rev.x = rev.x,
           zoom = list(xz = xz, yz = yz, xc = xc, yc = yc))
    }                                   # end of updatePlot
    
    ## --=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--
    ## Zoom control window
    zoomControl <- function(){
      ## Initialize zoom sliders
      visibleZoom <<- TRUE
      zoomX <- gWidgets::gslider(from=1,to=100,by=.1, value=1, handler = function(h,...){ xz <<- gWidgets::svalue(zoomX); updatePlot()})
      zoomY <- gWidgets::gslider(from=1,to=100,by=.5, value=1, handler = function(h,...){ yz <<- gWidgets::svalue(zoomY); updatePlot()})
      centerX <- gWidgets::gslider(from=-100,to=100,by=.1, value=0, handler = function(h,...){ xc <<- gWidgets::svalue(centerX); updatePlot()})
      centerY <- gWidgets::gslider(from=-100,to=100,by=.1, value=0, handler = function(h,...){ yc <<- gWidgets::svalue(centerY); updatePlot()})
      resetZoom <- gWidgets::gbutton(text = "Reset zoom and center", handler = function(h,...){ gWidgets::svalue(zoomX)<-1; gWidgets::svalue(zoomY)<-1; gWidgets::svalue(centerX)<-0; gWidgets::svalue(centerY)<-0; updatePlot()})
      gridCheck <- gWidgets::gcheckbox('Grid', handler = function(h,...){ gridOn <<- gWidgets::svalue(gridCheck); updatePlot()})
      
      ## Make zoom window
      zoomWindow <- gWidgets::gwindow("Plot properties", width=300)
      superGroup <- gWidgets::ggroup(horizontal=FALSE,container=zoomWindow)
      
      ## Add zoom sliders
      #	gWidgets::add(superGroup,gWidgets::glabel("X"),expand=TRUE)
      subgroupXz <- gWidgets::gframe("X zoom",horizontal=FALSE)
      gWidgets::add(subgroupXz,zoomX,expand=TRUE)
      subgroupXc <- gWidgets::gframe("X center",horizontal=FALSE)
      gWidgets::add(subgroupXc,centerX,expand=TRUE)
      gWidgets::add(superGroup,subgroupXz,expand=TRUE)
      gWidgets::add(superGroup,subgroupXc,expand=TRUE)
      gWidgets::addSpace(superGroup,20,horizontal=FALSE)
      #	gWidgets::add(superGroup,gWidgets::glabel("Y"),expand=TRUE)
      subgroupYz <- gWidgets::gframe("Y zoom",horizontal=FALSE)
      gWidgets::add(subgroupYz,zoomY,expand=TRUE)
      subgroupYc <- gWidgets::gframe("Y center",horizontal=FALSE)
      gWidgets::add(subgroupYc,centerY,expand=TRUE)
      gWidgets::add(superGroup,subgroupYz,expand=TRUE)
      gWidgets::add(superGroup,subgroupYc,expand=TRUE)
      subgroup3 <- gWidgets::ggroup(horizontal=TRUE,expand=TRUE)
      gWidgets::add(subgroup3,resetZoom,expand=TRUE)
      gWidgets::add(subgroup3,gridCheck,expand=FALSE)
      gWidgets::add(superGroup,subgroup3,expand=TRUE)
      gWidgets::addhandlerdestroy(zoomWindow, handler=function(h,...){visibleZoom <<- FALSE})
    }                                   # end of zoomControl
    
    ## --=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--
    ## Apply to all spectra
    exportControl <- function(){
      visibleExport <<- TRUE
      spec <- list()
      command <- "spec <- baseline(Y"
      for(i in 1:dim(bAGUI[[method]])[1]){
        command <- paste(command, ", ", rownames(bAGUI[[method]])[i], "=", bAGUI[[method]]$current[i], sep="")
      }
      command <- paste(command, ", method='", method, "')", sep="")
      eval(parse(text=command))
      
      exportName <- gWidgets::gedit(text="corrected.spectra", width=20)
      doExport   <- gWidgets::gbutton(text = "Apply and export", handler = function(h,...){the.name <- gWidgets::svalue(exportName); cat("\nCorrecting ..."); putBaselineEnv('the.export', spec);
                                                                                           eval(parse(text = paste(the.name, ' <- getBaselineEnv("the.export")', sep="")),envir = .GlobalEnv)
                                                                                           #assign(the.name, the.export,envir = .GlobalEnv);
                                                                                           gWidgets::dispose(exportWindow);cat("\nSaved as: ",the.name, sep="")})
      exportWindow <- gWidgets::gwindow("Apply correction to all spectra", width=300)
      superGroup   <- gWidgets::ggroup(horizontal=FALSE,container=exportWindow)
      subgroup     <- gWidgets::gframe("Object name",horizontal=FALSE)
      gWidgets::add(subgroup,exportName,expand=TRUE)
      gWidgets::add(subgroup,doExport,expand=TRUE)
      gWidgets::add(superGroup,subgroup,expand=TRUE)
      gWidgets::addhandlerdestroy(exportWindow, handler=function(h,...){visibleExport <<- FALSE})
    }
    
    ## --=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--
    ## Slider configuration window
    parameterControl <- function(){
      ## Initialize spectrum type chooser and parameter fields
      visibleCustomize <<- TRUE
      
      saveParameters <- function(){
        theSet <- bAGUI[method][[1]]
        for(i in 1:length(parameterList)){
          for(j in 1:4){
            theSet[i,j] <- gWidgets::svalue(parameterList[[i]][[j]])
          }
        }
        bAGUI[method][[1]][,1:3] <<- theSet[,1:3]
        bAGUI[method][[1]]$current <<- theSet[,4]
        gWidgets::delete(outerParam,remParam); createMethodSliders(); gWidgets::dispose(parameterWindow)# ; baseline.compute()
      }                               # end of saveParameters
      
      addParameterGroup <- function(nameStr, lineNo){
        parameterList[[lineNo]]  <<- c(gWidgets::gedit(text = "", width=1,coerce.with=as.numeric),gWidgets::gedit(width=5,coerce.with=as.numeric),gWidgets::gedit(width=10,coerce.with=as.numeric),gWidgets::gedit(width=15,coerce.with=as.numeric))
        parameterGroup[lineNo+1,1] <<- gWidgets::glabel(text=nameStr)
        size(parameterList[[lineNo]][[1]]) <- c(60,20)
        size(parameterList[[lineNo]][[2]]) <- c(60,20)
        size(parameterList[[lineNo]][[3]]) <- c(60,20)
        size(parameterList[[lineNo]][[4]]) <- c(60,20)
        parameterGroup[lineNo+1,2] <<- parameterList[[lineNo]][[1]]
        parameterGroup[lineNo+1,3] <<- parameterList[[lineNo]][[2]]
        parameterGroup[lineNo+1,4] <<- parameterList[[lineNo]][[3]]
        parameterGroup[lineNo+1,5] <<- parameterList[[lineNo]][[4]]
      }                               # end of gWidgets::addParameterGroup
      
      setParameters <- function(){
        theSet <- bAGUI[method][[1]]
        for(i in 1:length(parameterList)){
          for(j in 1:4){
            gWidgets::svalue(parameterList[[i]][[j]]) <<- theSet[i,j]
          }
        }
      }
      
      saveButton <- gWidgets::gbutton(text = "Apply parameters", handler = function(h,...) saveParameters())
      
      ## Make parameter window
      parameterWindow <- gWidgets::gwindow("Slider configuration", width=300)
      superGroup <- gWidgets::ggroup(horizontal=FALSE,container=parameterWindow)
      #        choiceGroup <- gWidgets::gframe("Type of spectra", container=superGroup, horizontal=FALSE)
      #        gWidgets::add(choiceGroup,typeChooser,expand=FALSE)
      gWidgets::addSpace(superGroup,10,horizontal=FALSE)
      
      ## Add edit fields
      parameterList <- list()
      parameterGroup <- gWidgets::glayout(homogeneous = FALSE, spacing = 5, container=superGroup)
      parameterGroup[1,1] <- gWidgets::glabel("")
      parameterGroup[1,2] <- gWidgets::glabel("From:")
      parameterGroup[1,3] <- gWidgets::glabel("To:")
      parameterGroup[1,4] <- gWidgets::glabel("Spacing:")
      parameterGroup[1,5] <- gWidgets::glabel("Start:")
      nameStrs <- rownames(bAGUI[method][[1]])
      for(i in 1:length(nameStrs))
        addParameterGroup(nameStrs[i],i)
      setParameters()
      gWidgets::visible(parameterGroup) <- TRUE
      gWidgets::addSpring(superGroup)
      
      ## Add buttons
      subgroupButtons <- gWidgets::ggroup(horizontal=TRUE)
      gWidgets::add(subgroupButtons,saveButton,expand=FALSE)
      gWidgets::add(superGroup,subgroupButtons,expand=FALSE)
      gWidgets::addhandlerdestroy(parameterWindow, handler=function(h,...){visibleCustomize <<- FALSE})
    }                                   # end of parameterControl
    
    ## --=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--
    ## Set up sliders according to chosen baseline correction method
    sVals <- list()
    sliders <- resets <- tmps <- list()
    iI <- 0
    
    createMethodSliders <- function(){
      remParam <<- gWidgets::ggroup(horizontal=FALSE, container=outerParam)
      
      #############
      ## General ##
      #############
      
      ## Collect values for sliders
      sVals <<- bAGUI[[method]]
      sliders <<- resets <<- tmps <<- list()
      iI <<- dim(sVals)[1]
      for(i in 1:iI){
        sliders[[i]] <<- gWidgets::gslider(from=sVals[i,1], to=sVals[i,2], by=sVals[i,3], value=bAGUI[[method]]$current[i],
                                           handler = function(h,...){ for(a in 1:iI){bAGUI[[method]]$current[a] <<- gWidgets::svalue(sliders[[a]]); bc <- getBaselineEnv("baseline.current");bc$parValues[a] <- gWidgets::svalue(sliders[[a]]);putBaselineEnv("baseline.current",bc)}})
        resets[[i]]  <<- gWidgets::gbutton(text = "Reset", handler = function(h,...){for(a in 1:iI) gWidgets::svalue(sliders[[a]])<<-sVals[a,4]})
        tmps[[i]]    <<- gWidgets::gframe(paste(sVals[i,6], " (", rownames(sVals)[i], ")", sep=""), container=remParam, horizontal=TRUE)
        gWidgets::add(tmps[[i]], sliders[[i]], expand=TRUE); gWidgets::add(tmps[[i]], resets[[i]], expand=FALSE)
      }
    }                                   # end of createMethodSliders
    
    
    
    ##
    ## Create main GUI window
    ##
    
    
    ## Initialize spectrum chooser slider and method chooser
    if(dim(Y)[1] > 1)
      spectrumNo <- gWidgets::gslider(from=1, to=dim(Y)[1], by=1, value=1, handler = function(h,...) specNo<<-gWidgets::svalue(spectrumNo))
    if(exists("baselineAlgorithms",envir=.GlobalEnv)){
      bA <- get("baselineAlgorithms",envir=.GlobalEnv)
    } else {
      bA <- baselineAlgorithms
    }
    GUI.names <- sort(names(bAGUI))
    names <- character(length(GUI.names))
    for(i in 1:length(GUI.names)){ # Let bAGUI control, and bA have descriptions -------------
                                   names[i] <- paste("'", ifelse(is.null(bA[[GUI.names[i]]]@description),"",bA[[GUI.names[i]]]@description), " (", GUI.names[i], ")'", sep="")
    }
    methodChooser <- gWidgets::gdroplist(names,
                               selected=which(GUI.names==method), handler = function(h,...){method <<- GUI.names[gWidgets::svalue(methodChooser,index=TRUE)]; gWidgets::delete(outerParam,remParam); createMethodSliders(); setZoom(); baseline.compute()})
    
    ## Initialize window and main containers
    window <- gWidgets::gwindow("Baseline correction", width=300)
    superGroup2 <- gWidgets::ggroup(horizontal=FALSE,container=window)
    if(dim(Y)[1] > 1){
      tmp <- gWidgets::gframe("Spectrum number", container=superGroup2, horizontal=TRUE)
      gWidgets::add(tmp,spectrumNo, expand=TRUE)
    }
    gWidgets::add(superGroup2,methodChooser, expand=FALSE)
    gWidgets::addSpring(superGroup2)
    Settings <- gWidgets::ggroup(container=superGroup2, horizontal=FALSE)
    outerParam <- gWidgets::ggroup(horizontal=FALSE, container=Settings)
    remParam <- gWidgets::ggroup()
    
    createMethodSliders()               # Add algorithm slides
    
    ## Add buttons
    gWidgets::addSpring(superGroup2)
    buttonGroup <- gWidgets::ggroup(horizontal=TRUE, container=superGroup2)
    plotButton <- gWidgets::gbutton(text = "Update plot", handler=function(h,...) baseline.compute())
    gWidgets::add(buttonGroup,plotButton,expand=TRUE)
    newZoom <- gWidgets::gbutton(text = "Zoom", handler = function(h,...){ if(visibleZoom==FALSE) zoomControl() })
    gWidgets::add(buttonGroup,newZoom,expand=FALSE)
    newParameter <- gWidgets::gbutton(text = "Customize", handler = function(h,...){ if(visibleCustomize==FALSE) parameterControl()})
    gWidgets::add(buttonGroup,newParameter,expand=FALSE)
    newExport <- gWidgets::gbutton(text = "Apply to all", handler = function(h,...){ if(visibleExport==FALSE) exportControl()})
    gWidgets::add(buttonGroup,newExport,expand=FALSE)
    
    ## Display initial plot
    plot(0,0, xlim=c(-1,1), ylim=c(-1,1), xlab="", ylab="", main="", axes=FALSE, col='white')
    baseline.compute()
  } else {                                       # end of baselineGUI
    warning('Package gWidgets not installed')
    return(list())
  }
}
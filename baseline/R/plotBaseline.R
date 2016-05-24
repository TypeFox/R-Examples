## $Id: plotBaseline.R 188 2011-08-11 07:50:28Z kristl $
### Plot function for estimated baseline (result of call to baseline).

###
### The workhorse plot functions
###


## (The y and ... arguments are needed to register this as a plot method)
## FIXME: How can we avoid having y as the second argument?  Or can we use it
## for something useful?
plotBaseline <- function(x, y, specNo = 1, grid = FALSE, labels = 1:n,
                         rev.x = FALSE,
                         zoom = list(xz = 1, yz = 1, xc = 0, yc = 0),
                         ...) {
  ## Warning if y is specified
  if(!missing(y)) warning("Argument 'y' is ignored")
  
  spectrum  <- getSpectra(x)[specNo,]
  baseline <- getBaseline(x)[specNo,]
  corrected <- getCorrected(x)[specNo,]
  n <- length(spectrum)
  if (is.numeric(labels)) {
    xvals <- labels
    xaxt <- par("xaxt")
  } else {
    xvals <- 1:n
    xaxt <- "n"
  }
  
  par(bg="white", cex.main=0.7, pty="m", mar=c(3,2,2,1), mfrow=c(2,1))
  
  ## Vertical placing of 'Original spectrum'
  yu <- max(spectrum); yl <- min(spectrum)
  p  <- yu-yl
  ylim <- c(-p*0.025/zoom$yz + p*zoom$yc/100 + yl,
            p*1.025/zoom$yz + p*zoom$yc/100 + yl)
  
  ## Vertical placing of 'Baseline corrected spectrum'
  yuc <- max(corrected); ylc <- min(corrected)
  pc  <- yuc-ylc
  ylimc <- c(-pc*0.025/zoom$yz + pc*zoom$yc/100 + ylc,
             pc*1.025/zoom$yz + pc*zoom$yc/100 + ylc)
  
  ## Horizontal placing of both spectra
  xr <- range(xvals)
  xc <- mean(xr)
  xhwd <- diff(xr) / 2
  xlim <- (xr - xc)/zoom$xz + xc + xhwd * zoom$xc/100
  if(rev.x) xlim <- rev(xlim)
  
  ## Original spectrum and baseline
  plot(xvals, spectrum, type = "l", xlim = xlim, ylim = ylim, xaxt = xaxt,
       xlab = "", ylab = "", main = "Original spectrum")
  lines(xvals, baseline, col = 2)     # Baseline
  if (!is.numeric(labels)) {
    ticks <- axTicks(1)
    ## Remove tick marks outside the plot range:
    ticks <- ticks[ticks >= 1 & ticks <= length(labels)]
    axis(1, ticks, labels[ticks], ...)
  }
  if(grid) grid()
  
  ## Corrected spectrum and baseline
  plot(xvals, corrected, type = "l", xlim = xlim, ylim = ylimc, xaxt = xaxt,
       xlab = "", ylab = "", main = "Baseline corrected spectrum")
  lines(xvals, numeric(n), col = 2)   # Baseline
  if (!is.numeric(labels)) axis(1, ticks, labels[ticks], ...)
  if(grid) grid()
}


###
### The plot method
###

setMethod("plot", "baseline",
          function(x, y, specNo = 1, grid = FALSE, labels = 1:n,
                   rev.x = FALSE, zoom = NULL,
                   ...) {
            ## Warning if y is specified
            if(!missing(y)) warning("Argument 'y' is ignored")
            n <- dim(getSpectra(x))[2]
            if (!isTRUE(zoom)) {                # Just do a single plot
              if (missing(zoom)) zoom <- list(xz = 1, yz = 1, xc = 0, yc = 0)
              plotBaseline(x, specNo = specNo, grid = grid,
                           labels = labels, rev.x = rev.x,
                           zoom = zoom)
            } else {                        # Set up a zoom window
              if(requireNamespace("gWidgets", quietly = TRUE)){
                ## Parametres used by updatePlot and GUI:
                xz <- 1; yz <- 1; xc <- 0; yc <- 0
                if (length(dev.list()) == 0) dev.new() # Make sure a device is open
                mydevice <- dev.cur()
                updatePlot <- function() {
                  currdevice <- dev.cur()
                  if (currdevice != mydevice) dev.set(mydevice)
                  plotBaseline(x, specNo = specNo, grid = grid,
                               labels = labels, rev.x = rev.x,
                               zoom = list(xz = xz, yz = yz, xc = xc, yc = yc))
                  if (currdevice != mydevice) dev.set(currdevice)
                }
                ## Initialize zoom sliders
                zoomX <- gWidgets::gslider(from=1,to=100,by=.1, value=1, handler = function(h,...){ xz <<- gWidgets::svalue(zoomX); updatePlot()})
                zoomY <- gWidgets::gslider(from=1,to=100,by=.5, value=1, handler = function(h,...){ yz <<- gWidgets::svalue(zoomY); updatePlot()})
                centerX <- gWidgets::gslider(from=-100,to=100,by=.1, value=0, handler = function(h,...){ xc <<- gWidgets::svalue(centerX); updatePlot()})
                centerY <- gWidgets::gslider(from=-100,to=100,by=.1, value=0, handler = function(h,...){ yc <<- gWidgets::svalue(centerY); updatePlot()})
                resetZoom <- gWidgets::gbutton(text = "Reset zoom and center", handler = function(h,...){ gWidgets::svalue(zoomX)<-1; gWidgets::svalue(zoomY)<-1; gWidgets::svalue(centerX)<-0; gWidgets::svalue(centerY)<-0; updatePlot()})
                gridCheck <- gWidgets::gcheckbox('Grid', handler = function(h,...){ grid <<- gWidgets::svalue(gridCheck); updatePlot()})
                
                ## Make zoom window
                zoomWindow <- gWidgets::gwindow(paste("Zoom: Device", mydevice), width=300)
                superGroup <- gWidgets::ggroup(horizontal=FALSE,container=zoomWindow)
                
                ## Add zoom sliders
                #	gWidgets::add(superGroup,glabel("X"),expand=TRUE)
                subgroupXz <- gWidgets::gframe("X zoom",horizontal=FALSE)
                gWidgets::add(subgroupXz,zoomX,expand=TRUE)
                subgroupXc <- gWidgets::gframe("X center",horizontal=FALSE)
                gWidgets::add(subgroupXc,centerX,expand=TRUE)
                gWidgets::add(superGroup,subgroupXz,expand=TRUE)
                gWidgets::add(superGroup,subgroupXc,expand=TRUE)
                addSpace(superGroup,20,horizontal=FALSE)
                #	gWidgets::add(superGroup,glabel("Y"),expand=TRUE)
                subgroupYz <- gWidgets::gframe("Y zoom",horizontal=FALSE)
                gWidgets::add(subgroupYz,zoomY,expand=TRUE)
                subgroupYc <- gWidgets::gframe("Y center",horizontal=FALSE)
                gWidgets::add(subgroupYc,centerY,expand=TRUE)
                gWidgets::add(superGroup,subgroupYz,expand=TRUE)
                gWidgets::add(superGroup,subgroupYc,expand=TRUE)
                subgroup3 <- ggroup(horizontal=TRUE,expand=TRUE)
                gWidgets::add(subgroup3,resetZoom,expand=TRUE)
                gWidgets::add(subgroup3,gridCheck,expand=FALSE)
                gWidgets::add(superGroup,subgroup3,expand=TRUE)
                updatePlot()
              } else {
                warning('Package gWidgets not installed')
                return(list())
              }
            }
          }
)

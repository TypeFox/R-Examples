# optionally, graphs can be plotted with zplot and lines with zlines
# allowing to go back and forth with bf()

# to use it: copy it in your folder and source it:
# source("zoom.r")

# WARNING: when you zoom on a multiplot window, you can only zoom on the last plot and the other plots will be rescaled
# accordingly, can be very usefull or anoying but I don't plan to change this soon

#=============================
# Navigation "session"
#=============================
# multipancPoint
# @rdname replot

multipancPoint<-function(ancien,fact=1,new,point=NULL){
	# cat("ancien",ancien,"fact",fact,"new:\n")
	# print(new)
  	if(is.null(point)) point<-mean(ancien)
	fact<-1/fact
	newRange<- (1-fact)*point+fact*ancien

	return(newRange);
}
centerOnPoint<-function(ancien,fact=1,new,point=NULL){
  	currentCenter<-mean(ancien)
	if(is.null(point)) point<-currentCenter
	fact<-1/fact
	newRange<- fact*ancien-currentCenter+point

	return(newRange);
}

keepanc<-function(ancien,fact,new,point=NULL){
	# cat("ancien",ancien,"fact",fact,"new:\n")
	# print(new)
	return(ancien);
}
usenew<-function(ancien,fact,new,point=NULL){
	# cat("ancien",ancien,"fact",fact,"new:\n")
	# print(new)
	return(new);
}
# to avoid repeating oneself, specially on something quickly changing
# in case of need for change here, also check "locator" in zoomplot.zoom()
is.plot.window<-function(alst,fn){
  ## for former versions of R
  # tmp <- all.equal(.Primitive("plot.window"), fn)
  # tmp <- (length(grep("plot.window",deparse(fn)))>0)
  version<-sessionInfo()$R.version
  Maj<-as.numeric(version$major)
  Min<-as.numeric(version$minor)
  if(Maj>=3 & Min>= 0.1){
    tmp <- all.equal("C_plot_window",alst[[1]]$name)
    attributes(tmp)$lims<-c(2,3)
  }else if(Maj==3 & Min <=0.1){
    tmp <- (length(grep("plot.window",deparse(fn)))>0)
    attributes(tmp)$lims<-c(1,2)
  }else if(Maj<3){
    tmp <- all.equal(.Primitive("plot.window"), fn)
    attributes(tmp)$lims<-c(1,2)
  }
  return(tmp)
}
is.locator<-function(alst,fn){
  ## for former versions of R
  ## beginning 3.0.1
  version<-sessionInfo()$R.version
  Maj<-as.numeric(version$major)
  Min<-as.numeric(version$minor)
  if(Maj>=3 & Min>= 0.1){
    tmp <- all.equal("C_locator",alst[[1]]$name)
  }else if(Maj==3 & Min <=0.1){
    tmp <- (length(grep("locator",deparse(fn)))>0)
  }else if(Maj<3){
    tmp <- all.equal(.Primitive("locator"), fn)
  }
  return(tmp)
}

# get the limits of the plot in argument
getalst<-function(tmp=recordPlot()[[1]]){
	for (i in seq(along = tmp)) {
		fn <- tmp[[i]][[1]]
		alst <- as.list(tmp[[i]][[2]])
		tmp2<- is.plot.window(alst,fn)
		if (is.logical(tmp2) && tmp2) {
			alstwin<-alst
		}
	}
	return(alstwin);
}

#' Central low level function of the zoom package.
#'
#' This function allow to replot the current or a saved plot with specific
#' boundaries, magnification factor and possibly arround a user defined x/y.
#'
#' This function is not necessarily easy to use by hand. It is designed to work
#' well when called from higher level functions. End user should always use
#' zm().
#'
#' @param xlim A vector with min and max x
#' @param ylim A vector with min and max y
#' @param fact A scalar giving the magnification factor (>1 brings you closer)
#' @param rp A previously recorded plot with recordPlot(). With all the
#' corresponding warnings in ?recordPlot.
#' @param x x of a fix point when rescaling, by default the center.
#' @param y y of a fix point when rescaling, by default the center.
#' @param xlimfn a function using x, y and/or fact to generate new x lim if NULL and xlim/ylim not given will use multipancPoint
#' @param ylimfn a function using x, y and/or fact to generate new y lim, if NULL will use xlimfn
#' @param \dots Additional parameters not implemented, just in case.
#' @return Not guaranted for now.
#' @note This function is the heart of the zoom package and the one that can be.
#' affected by R version changes.
#' It is inspired by the zoomplot function in TeachingDemos package
#' @author Corentin M. Barbu
#' @seealso zm, in.zoom
#' @keywords zoom plot
#' @examples
#'
#' plot(rnorm(1000),rnorm(1000))
#' zoomplot.zoom(fact=2,x=0,y=0)
#'
#' @export zoomplot.zoom
zoomplot.zoom <- function (xlim=NULL, ylim = NULL,fact=NULL,rp=NULL,x=NULL,y=NULL,xlimfn=NULL,ylimfn=NULL,...)
{
  # cat("using zoomplot.zoom")
  # rp is a recorded plot
  # fact is a factor of magnification/outzoom
  # fact has priority on xlim/ylim
  if(is.null(xlimfn)){
    if(! is.null(fact)&& is.numeric(fact) && length(fact)==1){
      xlimfn <-multipancPoint;
      ylimfn <-multipancPoint;
    }else{
      if(!is.null(xlim)&& is.numeric(xlim) && length(xlim)==2){
	xlimfn <- usenew;
      }else{
	xlimfn <-keepanc;
      }
      if(!is.null(ylim) && is.numeric(ylim) && length(ylim)==2){
	ylimfn <-usenew;
      }else{
	ylimfn <-keepanc;
      }
    }
  }
  if(is.null(ylimfn)){
    ylimfn<-xlimfn
  }

  if(is.null(rp)){
    tmp <- recordPlot()[[1]]
  }else{
    tmp<-rp[[1]]
  }

  plotOk<-NULL
  for (i in seq(along = tmp)) {
    # cat("i:",i,"\n")
    fn <- tmp[[i]][[1]]
    alst <- as.list(tmp[[i]][[2]])
    tmp1 <- is.locator(alst,fn)
    if (is.logical(tmp1) && tmp1) {
      next # will not like do.call
    }
    tmp2<- is.plot.window(alst,fn)
    if (is.logical(tmp2) && tmp2) {
      # print(alst)
      # cat("alst orig:",alst[[1]],alst[[2]],"\n")
      locx<-attributes(tmp2)$lims[1]
      locy<-attributes(tmp2)$lims[2]
      alst[[locx]] <- xlimfn(alst[[locx]],fact,xlim,x)
      alst[[locy]] <- ylimfn(alst[[locy]],fact,ylim,y)
    }
    plotOk<-try(do.call(fn, alst))
  }
  if(isError(plotOk)){
    box() # finish the graph even if errors in between
  }
}

is.out.of.plot.click<-function(loc){
	plotLim<-par("usr")

	if(loc$x<min(plotLim[1:2])
	|| loc$x>max(plotLim[1:2])
	|| loc$y>max(plotLim[3:4])
	|| loc$y<min(plotLim[3:4])){
		return(TRUE)
	}else{
		return(FALSE)
	}
}
# to avoid repeating this painful piece
other.option.session.message<-function(){
  devType<-names(dev.cur())
  terminate.key<-switch(EXPR=devType,
	 X11cairo = "Any other click",
	 X11 = "Any other click",
	 quartz = "Hit Esc",
	 RStudioGD = "Click on Finish",
  	 windows = "Middle click","Terminate as in locator()")

  cat(terminate.key,"for other options.\n")
}

#' @title Direct access to zoom functionalities.
#'
#' @description Direct selection of a zoom method of the "session" type.
#' Possibly of use in scripts?
#'
#' @aliases in.zoom out.zoom set.zoom sq.zoom
#' @param \dots Extra arguments to zoomplot.zoom.
#' @note Each function starts a different interactive sequence
#'  	\itemize{
#' 	 	\item{inout.zoom(): }{left click within bounds zooms in, outside bounds zoom out}
#' 	 	\item{move.to.click.zoom(): }{center plot arround left click}
#' 	 	\item{in.zoom(): }{each left click zooms in}
#' 	 	\item{out.zoom(): }{each left click zooms out}
#' 	 	\item{set.zoom(): }{ask for a magnification factor}
#' 		\item{sq.zoom(): }{allow to click on the two corners of the desired region to zoom on}
#'	}
#' @return NULL
#' @author Corentin M. Barbu
#' @seealso zm(), session.zoom().
#' @export in.zoom
in.zoom<-function(...){
  # Ideally later should center arround the point selected
  cat("Left click to zoom in\n")
  other.option.session.message()

  center<-locator(1)
  if(length(center$x)==1){
    zoomplot.zoom(fact=2,x=center$x,y=center$y,...);
    in.zoom()
  }
  return()
}
#' @rdname in.zoom
#' @export move.to.click.zoom
# center the plot at the point of click
move.to.click.zoom<-function(...){
  cat("Click on the point you want to center the plot arround.\n")
  other.option.session.message()

  center<-locator(1)
  if(length(center$x)==1){
    zoomplot.zoom(fact=1,x=center$x,y=center$y,xlimfn=centerOnPoint,...);
    move.to.click.zoom(...)
  }
  return()
}

#' allow interactive in/out zoom in "session" mode
#' @rdname in.zoom
#' @export inout.zoom
inout.zoom<-function(...){
  # Ideally later should center arround the point selected
  cat("Left click in plot to zoom in\n")
  cat("Left click out of plot to zoom out\n")
  other.option.session.message()

  center<-locator(1)
  if(length(center$x)==1){
	  if(is.out.of.plot.click(center)){
		  zoomplot.zoom(fact=0.5,...);
	  }else{
		  zoomplot.zoom(fact=2,x=center$x,y=center$y,...);
	  }
    inout.zoom(...)
  }
  return()
}
#' @rdname in.zoom
#' @export out.zoom
out.zoom<-function(...){
  # Ideally later should center arround the point selected
  cat("Left click to zoom out\n")
  other.option.session.message()

  center<-locator(1)
  print(center)
  if(length(center$x)==1){
    zoomplot.zoom(fact=0.5,x=center$x,y=center$y,...);
    out.zoom()
  }
  return()
}
#' @rdname in.zoom
#' @export out.zoom
set.zoom<-function(...){
	cat("Enter magnification factor: \n")
	f<-scan(n=1)
	zoomplot.zoom(fact=f,...);
	return()
}

# Double-click detection
# Test is the output of locator corresponds to two identical points
is.double.click<-function(loc){
	if(length(loc$x)!=2){
		cat("1\n")
		return(FALSE)
	}else{
		if(loc$x[1]==loc$x[2] &&
		loc$y[1]==loc$y[2]){
			return(TRUE)
		}else{
		cat("2\n")
		print(loc)
			return(FALSE)
		}
	}
}

#' @rdname in.zoom
#' @export sq.zoom
sq.zoom<-function(...){
	# use locator to zoom with the mouse (two left clicks)
	# specially, ... can be used to pass a recorded plot rp
	cat("Click left over opposite corners of zoom area.\n");
	cat("Double left click for zoom out\n")
	other.option.session.message()
	square<-locator(2)
	print(square)
	if(length(square)==2){
	  if(is.double.click(square)){
		zoomplot.zoom(fact=0.5,x=square$x,y=square$y,...);
	  	out.zoom()
	  }else{
		  xmin<-min(square$x)
		  xmax<-max(square$x)
		  ymin<-min(square$y)
		  ymax<-max(square$y)
		  zoomplot.zoom(xlim=c(xmin,xmax),ylim=c(ymin,ymax),...)
	  }
	  sq.zoom(...)
	}
}
orig.zoom<-function(orig){
	replayPlot(orig)
	return()
}

#' @import tools
print.zoom<-function(orig=NULL,dev=NULL,fileName=NULL,...){
  if(is.null(fileName)){
    fileName<-file.choose()
  }
  if(is.null(dev)){
    if (isError(try(dev<-eval(parse(text=tools::file_ext(fileName))),
                    silent=TRUE))) {
      cat("Error: extension not recognized, try png or pdf\n")
      return(1)
    }
  }

  if(identical(dev,pdf)){
    dev.print(device=dev,fileName,...)
  }else{
    devSize<-dev.size(units="px")
    dev.print(device=dev,fileName,width=devSize[1],height=devSize[2],...)
  }
  cat("Successfully printed in",fileName,"\n")
}

png.zoom<-function(orig=NULL,fileName=NULL,...){
  devSize<-dev.size(units="px")
  print.zoom(dev=png,fileName,width=devSize[1],height=devSize[2],...)
}

pdf.zoom<-function(orig=NULL,fileName=NULL){
  print.zoom(dev=pdf)
}

#' @rdname navigation.zoom
#' @export session.zoom
session.zoom<-function(...){
	orig <- recordPlot()
	go_on<-TRUE
	inout.zoom(...)
	while(go_on){
		cat("Do you want to?\n")
		cat(" zoom in/out: 1\n")
		cat("move arround: 2\n")
		cat(" zoom square: 3\n")
		cat("        save: 6\n")
		cat("back to init: 9\n")
		cat("        Exit: Enter\n")
		sel<-scan(n=1)
		if(length(sel)==0){
			go_on<-FALSE;
		}else{
			if(exists("exec.zoom")){
				rm(exec.zoom)
			}
			exec.zoom<-switch(sel,
				inout.zoom,move.to.click.zoom,
				sq.zoom,NULL,
				NULL,print.zoom,NULL,NULL,
				orig.zoom)
			if(!is.null(exec.zoom)){
				exec.zoom(orig=orig,...);
			}else{
				cat("Say it again?\n")
			}
		}
	}
	return(recordPlot());
}

# Transform the button event code passed by the event handler into a
# meaningful string.
#
# The event handler code returns a hard to understand code. I made a lot of
# trials and errors and came up with the corresponding function to triage
# these codes. This may evolve as I better understand this system.
#
# Used for all analysis of the events in navigation.zoom The returned values
# should be self-explaining.
#
# @param buttons button event code as passed to event functions by
# getGraphicsEvent()
# @return The returned value is a string that can be:\itemize{
# \item{left }{left button of the mouse}
# \item{right }{right button of the mouse}
# \item{middle }{middle button or scrolling weel *pressed*}
# \item{scrollDown }{scrolling weel turned down}
# \item{scrollUp }{scrolling weel turned up}
# }
# @author Corentin M. Barbu
# @seealso setCallBack
# @keywords event
# @examples
#
# labelButton(c(0,1)) # should return "right"
#
# @export labelButton
labelButton<-function(buttons){
  # label the buttons with easy to remember names
  label<-""
  # cat("buttons:")
  # print(buttons)
  if(length(buttons)>=2){ # rightbutton or scrolling
    if(buttons[2]==2){ # scroll down
      label<-"scrollDown"
    }else if(buttons[2]==1){ # right button
      label<-"right"
    }
  }else if(length(buttons)==1){
    if(buttons==1){ # middle button
      label<-"middle"
    }else if(buttons==0){
      label<-"left"
    }else if(buttons==2){# scroll up
      label<-"scrollUp"
    }
  }else{
    label<-NULL
  }
  # cat("label:",label,"\n")
  return(label)
}

setCallBack<-function(..., xlim = NULL, ylim = NULL, xaxs = "r", yaxs = "r"){
  startx <- NULL
  starty <- NULL
  usr <- NULL
  rp <- recordPlot()

  #---------------------
  # Navigation functions
  #---------------------
    dragmousemove <- function(buttons, x, y) {
    devset()
    # cat("In dragmousemove\n")
    deltax <- diff(grconvertX(c(startx, x), "ndc", "user"))
    deltay <- diff(grconvertY(c(starty, y), "ndc", "user"))
    xlim<<-usr[1:2]-deltax
    ylim <<-usr[3:4]-deltay
    zoomplot.zoom(xlim=xlim,ylim=ylim,...)
    NULL
  }
  zoomDyn <- function(buttons, x, y) {
    devset()
    usr <<- par("usr")
    mevent<-labelButton(buttons)
    if(!is.null(mevent)){
      # cat("mevent:",mevent,"\n")
      if(mevent=="scrollDown" || mevent=="middle"){
	fact<-0.7
      }else if(mevent=="scrollUp" || mevent == "right"){
	fact<-1.5
      }else {
	deltay <- diff(grconvertY(c(starty, y), "ndc", "user"))
	fact<-max(min(1+deltay/(usr[2]-usr[1]),10),0.1)
      }
      xuser<-grconvertX(x, "ndc", "user")
      yuser<-grconvertY(y, "ndc", "user")

      zoomplot.zoom(fact=fact,x=xuser,y=yuser)
      # cat("fact:",fact,"\n")
      # cat(usr[1:2],"->",xlim,"\n")
      # cat(usr[3:4],"->",ylim,"\n")
    }
    NULL
  }
  mouseDownNavig <- function(buttons, x, y) {
    startx <<- x
    starty <<- y
    devset()
    usr <<- par("usr")
    # cat("buttonPress:",buttons,"\n")
    mevent<-labelButton(buttons)
    if(mevent=="scrollDown"){
      eventEnv$onMouseMove <- zoomDyn
      zoomDyn(buttons,x,y)
    }else if(mevent=="scrollUp"){
      eventEnv$onMouseMove <- zoomDyn
      zoomDyn(buttons,x,y)
    }else if(mevent=="middle"){
      eventEnv$onMouseMove <- zoomDyn
      zoomDyn(buttons,x,y)
      # cat("Turn on zoomDyn\n")
    }else if(mevent=="left"){
      # cat("Turn on dragmousemove\n")
      eventEnv$onMouseMove <- dragmousemove
    }else if(mevent=="right"){
      # cat("Closing...")
      # return(invisible(1))
    }
    NULL
  }

  mouseup <- function(buttons, x, y) {
    eventEnv$onMouseMove <- NULL
    NULL
  }

  #---------------------
  # Keyboard mode commuter
  #---------------------

  keydown <- function(key) {
    switch(key,
      ## quit (\033 == ESC)
      "\033" =, "ctrl-C" =, "q" = { return(invisible(1)) },
      ## print to file
      "p" = {
        message("Entering printing mode:")
        eventEnv$prompt <<- "Printing mode"
        setGraphicsEventEnv(env=eventEnv)
        print.zoom()
      },
      ## show limits
      "s" = { message("xlim=c(",
                      paste0(round(.xlim(), digits=3), collapse=", "), "), ",
                      "ylim=c(",
                      paste0(round(.ylim(), digits=3), collapse=", "), ")") },
      ## restaure initial size
      "r" = { orig.zoom(rp) },
      ## zoom in (ctrl-* == [CTRL]+[+])
      "ctrl-*" =, "i" =,"+" = { zoomplot.zoom(fact=1.1) },
      ## zoom out (ctrl-_ == [CTRL]+[-])
      "ctrl-_" =, "o" =,"-" = { zoomplot.zoom(fact=0.9) },
      ## zoom in (x-axis only)
      "L" = { zoomplot.zoom(xlim=.zoomXlim(1.1)) },
      ## zoom out (x-axis only)
      "H" = { zoomplot.zoom(xlim=.zoomXlim(0.9)) },
      ## zoom in (y-axis only)
      "K" = { zoomplot.zoom(ylim=.zoomYlim(1.1)) },
      ## zoom out (y-axis only)
      "J" = { zoomplot.zoom(ylim=.zoomYlim(0.9)) },
      ## move left
      "Left" =, "h" = { zoomplot.zoom(xlim=.moveXlim(-0.1)) },
      ## move right
      "Right" =, "l" = { zoomplot.zoom(xlim=.moveXlim(+0.1)) },
      ## move down
      "Down" =, "j" = { zoomplot.zoom(ylim=.moveYlim(-0.1)) },
      ## move up
      "Up" =, "k" = { zoomplot.zoom(ylim=.moveYlim(+0.1)) },
      ## default (nothing)
      {}
    )
    return(NULL)
  }

  #---------------------
  # Set event handler
  #---------------------
  setGraphicsEventHandlers(prompt = "p to print, q to quit",
			   onMouseDown = mouseDownNavig,
			   onMouseUp = mouseup,
			   onMouseMove = NULL,
			   onKeybd = keydown)

  eventEnv <- getGraphicsEventEnv()

  devset <- function(){
    if (dev.cur() != eventEnv$which) dev.set(eventEnv$which)
  }
}


#' Opening of an interactive zoom/navigate session.
#'
#' To launch an interactive session you should use zm() but if you are sure of
#' your device you can launch directly one of these functions.
#'
#' session.zoom launch an interactive console menu to navigate a plot.
#'
#' navigation.zoom allows to interactively navigate a plot with the mouse.
#'
#' @aliases navigation.zoom session.zoom
#' @param \dots Everything that can be accepted by sq.zoom.
#' @return Returns the final plot, as saved by recordPlot().
#' @author Corentin M. Barbu, Sebastian Gibb
#' @seealso zm().
#' @keywords session navigation
#' @examples
#' \dontrun{
#' plot(rnorm(100),rnorm(100))
#' session.zoom()
#' }
#'
#' @export navigation.zoom
navigation.zoom<-function(...){
  message("Mouse:")

  if(names(dev.cur())=="windows"){
    zoom.in.out.mes<-"Right to zoom in, Middle or Hold Left + click right to zoom out"
  }else{
    zoom.in.out.mes<-"Scroll to zoom in and out"
  }
  message(zoom.in.out.mes,"\nHold left mouse button to move")

  message("\nKeyboard:")
  ## keyboard usage
  keys <- c("Left/Right (h/l)", "Up/Down (k/j)",
            "+ or i / - or o", "L/H", "K/J",
            "p", "r", "s")
  usage <- c("move left/right", "move up/down",
             "zoom in/out", "zoom in/out (x-axis only)",
             "zoom in/out (y-axis only)",
             "print to file", "reset limits", "show limits")
  message(paste(format(keys, justify="left"),
                format(usage, justify="right"),
                sep=": ", collapse="\n"))

  g<-0
  while(length(g)!=1 || g!=1){
    g<-getGraphicsEvent(consolePrompt="q on the graphic window to quit")
  }
  out<-recordPlot()
  return(out)
}
# try to replot the graph into a Xlib device to allow events handling


# Zoom package utility functions
#
# Different functions, should not be used by end users, definitions may change
# in future versions.
#
# This functions are more or less tightly connected to R internals functions
# that are subject to a lot of changes and should not be relied on by a end
# user.  In particular is.plot.windown() is responsible for the delicate task
# of identifying in a recording of a plot to identify the part actually
# responsible for plotting the graphic window. This single thing has changed
# three times between 2010 and 2013.
#
# setCallBack() sets up the call handler. It should probably never be called
# alone.
#
# keepanc, usenew and multipancPoint are three functions possibly called by
# zoomplot.zoom() in an interchangeable way, none of them make use of all the
# arguments.  keepanc() and usenew() return ancien or new according to their
# names. multipancPoint() allows to transform "ancien" into new coordinates,
# after magnification by "fact" (ignores "new") but keeping point invariant.
#
# @aliases replot getalst is.plot.window is.locator keepanc usenew
# setCallBack orig.zoom multipancPoint
# @param rp Saved plot as generated by recordPlot().
# @param tmp Usable part (equivalent to rp[[1]]) of a saved plot as generated
# by recordPlot().
# @param alst Arguments part of an item in an recorded plot object.
# @param fn Function part of an item in an recorded plot object.
# @param ancien ancien set of coordinates
# @param fact factor of magnification
# @param new new set of coordinates
# @param point coordinates of the pointer
# @param xlim x limits of the plot
# @param ylim y limits of the plot
# @param xaxs style of the x axis see xaxs in par() help
# @param yaxs style of the y axis see yaxs in par() help
# @param orig like rp, orig.zoom is at this point just a convenience wrapper
# for replayPlot()
# @param ...  Possible addicional arguments to zoomplot.zoom
# @return Can't be relied upon at this point.
# @author Corentin M. Barbu, Sebastian Gibb
# @seealso zm, zoomplot.zoom
replot <- function(rp=NULL) {
  if (!interactive()) {
    stop("R has to run in interactive mode.")
  }

  if (is.null(rp)) {
    rp <- recordPlot()
  }

  initDev <- dev.cur()

  ## choose correct device for each supported platform
  if (.Platform$OS.type == "windows") {
    cl <- call("windows")
  } else {
    cl <- call("X11", type="Xlib")
  }

  initWarn<-options("warn")
  options(warn=2) # avoid overlooking a warning that may be deal killer
  if (!isError(try(eval(cl), silent=TRUE)) &&
      !isError(try(replayPlot(rp), silent=TRUE))) {
    message("Replot successful, use ", deparse(cl), " to avoid this step.")
    return(!isError(try(setCallBack(), silent=TRUE)))
  } else {
    if(dev.cur() != initDev) {
      dev.off()
    }
    message("Fall back to classical interface.")
    message("Use ", deparse(cl), " to enable full navigation.")
    return(FALSE)
  }
  options(warn=initWarn)
}

# zm() Main function, choosing between navigation or old "session" interactions

#' Launch interaction on a plot
#'
#' Allow to zoom/navigate in any open plot. The controls should be intuitive:
#' \itemize{
#'   \item{zoom in:}{ scroll up, or right click if no scrolling weel.}
#'   \item{zoom out:}{ scroll down, or Hold left + right click if no working weel.}
#'   \item{move:}{ left click and move }
#' }
#'
#' By default, zm() try to open a mouse interactive session. If the current
#' device is not interactive, will try to replot the current plot in a
#' \code{X11(type="Xlib")} device. If it fails it will open a console menu
#' based interactive session.
#'
#' Zoom handle multiple plots on a device together. You need to navigate the
#' last one plotted and all the other plots will be navigated according to the
#' last one: that can be pretty amazing too if you want to explore multiple
#' layers at the same time.
#'
#' @param type the type of interaction with the plot. Possible types are:
#' \itemize{
#'    \item{session}{ for console menu}
#'    \item{navigation}{ for mouse interaction}
#' }
#'
#' Or any short names for these. By default will try to
#' launch a "navigation" session.
#' @param rp plot to navigate, saved using \code{rp<-recordPlot()}. By default
#' (NULL) will use the current device.
#' @return The recording of the final plot. Can be reploted using replayPlot().
#' The most useful may be to get the xlim and ylim of the final plot. That can
#' be simply got using: \code{par("usr")} after \code{zm()} ends.
#' @note This function relies on pretty low level functions in R that change
#' quite often with new versions. New version of R can break this package but I
#' got used to it and fix it quickly.
#'
#' In case you close the device before striking q, just hit Ctrl-C on the
#' command line.
#' @author Corentin M. Barbu
#' @keywords zoom zm navigate navigation plot
#' @examples
#'
#' \dontrun{
#' # basic example
#' plot(rnorm(1000),rnorm(1000)) # could be any plot
#' zm() # navigate the plot
#'
#' # use the same xlim/ylim as ended up in the zoom session
#' xylim<-par("usr") # xmin,xmax,ymin,ymax of the final version of the plot
#' dev.off()
#' plot(rnorm(1000),rnorm(1000),xlim=xylim[1:2],ylim=xylim[3:4])
#'
#' # navigate two layers of data at the same time
#' par(mfrow=c(1,2))
#' plot(1,type="n",xlim=c(-3,3),ylim=c(-3,3),main="First Track")
#' polygon(c(-1,1,1,-1)*2,c(-1,-1,1,1)*2,col="blue")
#' lines(rnorm(100),rnorm(100))
#' plot(1,type="n",xlim=c(-3,3),ylim=c(-3,3),main="Second Track")
#' polygon(c(-1,1,1,-1)*2,c(-1,-1,1,1)*2,col="green")
#' lines(rnorm(100),rnorm(100))
#' zm() # it flickers quite a bit as it needs to replot everything every time...
#'
#' # one might want to use the older interface
#' # if attached to cairo under linux or MacOS
#' # it is also sometimes helpful to just define a square you want to zoom on
#' zm(type="s")
#' }
#'
#' @export zm
zm <- function(type="navigation", rp=NULL) {

  if (missing(type)) {
    if (isError(try(setCallBack(), silent=TRUE))) {
      if (replot(rp=rp)) {
        type <- "navigation"
      } else {
        type <- "session"
      }
    }
  }

  if (type == "navigation") {
    return(try(navigation.zoom()))
  } else {
    return(try(session.zoom()))
  }
}


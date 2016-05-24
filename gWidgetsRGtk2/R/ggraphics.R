## cairo graphics device
## would like to get size from par("fin"), but this isn't so easy as it
## seems to pop up a new plot container

## pas through ... some means to control:
## rubber banding (do.rubber.banding=FALSE)
## menu popup (no_popup=TRUE)

setClass("gGraphicsRGtk",
         contains="gComponentRGtk",
         prototype=prototype(new("gComponentRGtk"))
         )

setMethod(".ggraphics",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   width=dpi*6, height=dpi*6,
                   dpi=75, ps=12,
                   container=NULL,...) {

            force(toolkit)
            
            da <- gtkDrawingAreaNew()
            asCairoDevice(da, pointsize=ps)
            ## set size
            if(!is.null(width) & !is.null(height))
              da$setSizeRequest(width, height)              
            ## allow events on this widget
            da$AddEvents(GdkEventMask["all-events-mask"])

            
            obj <- as.gWidgetsRGtk2(da)
#            obj = new("gGraphicsRGtk",block=da, widget=da, toolkit=toolkit)


            ## Woah Nelly! since 2.0.1 the device needs to be realized before we can make it
            ## so we put this in: 
            ## when a device is clicked.
            
            addhandler(obj,signal="map-event",handler = function(h, w, e, ...) {
              da <- h$action
              ## in cairoDevice (>= 2.2.0) the device is stored in da$GetData(".devnum")
              if(is.null(da$GetData(".devnum"))) {
                asCairoDevice(da, pointsize=ps) # turn into cairo device
                tag(obj,"device") <- da$GetData(".devnum")
              }
              return(TRUE)             # don't propogate
             }, action=da)

            ## handlers to raise device when clicked upon. This seems a natural way to interact with
            ## the device
            .getDevNo <- function(da) da$getData(".devnum")
            .setDevNo <- function(da, ...) {
              dev.set(.getDevNo(da))
              ## indicate?
              
              FALSE}
            ## raise when click into window
            gSignalConnect(da, "button-press-event", f=.setDevNo)
            ## raise when motion over device -- CONFUSING, leave out
#            da$addEvents(GdkEventMask['enter-notify-mask'])
#            gSignalConnect(da, "enter-notify-event", f=.setDevNo)
            ## close device when destroyed
            gSignalConnect(da, "destroy-event", f=function(da, ...) {
              dev.off(.getDevNo(da))
              return(FALSE)
            })


            ## Add rubber banding
            ## This code is borrowed from the excellent playwith package by Felix Andrews
            theArgs <- list(...)
            doRubberBanding <- getWithDefault(theArgs$do.rubber.banding, TRUE)
            
            ## add environment and values to da
            e <- environment()
            e$dragging <- FALSE
            e$x0 <- e$y0 <- e$x <- e$y <- 0
            da$setData("env", e)

            ## need to bind drag actions: click, motion, release

            if(doRubberBanding)
              gSignalConnect(da, "button-press-event", function(w, e) {
                if(isRightMouseClick(e))
                  return(FALSE)
                da <- w
                daClearRectangle(da)

                wh <- daGetWidthHeight(da)
                da.w <- wh[1]
                da.h <- wh[2]
                
                ## assign("da", da, envir=.GlobalEnv)
                ## buf <- gdkPixbufGetFromDrawable(src=da$window, src.x=0, src.y=0,
                ##                                 dest.x=0, dest.y=0, width=da.w, height=da.h)
                
                ## w$setData("buf", buf)
                env <- w$getData("env")
                env$x0 <- env$x <- e$x
                env$y0 <- env$y <- e$y
                env$dragging <- TRUE
                return(FALSE)
              })

            if(doRubberBanding)
              gSignalConnect(da, "motion-notify-event", function(w, e) {
                env <- w$getData("env")
                ## are we dragging?
                if(env$dragging) {
                  daClearRectangle(w)
                  
                  env$x <- e$x
                  env$y <- e$y
                  
                  ## did we move enough? 10 pixels say
                  
                  if(max(abs(env$x - env$x0), abs(env$y - env$y0)) > 10)
                  daDrawRectangle(w, env$x0, env$x, env$y0, env$y)
                  
                  
                }
                return(FALSE)
              })
            
            if(doRubberBanding)            
              gSignalConnect(da, "button-release-event", function(w, e) {
                if(isRightMouseClick(e))
                  return(FALSE)
                env <- w$getData("env")
                ## remove draggin
                env$dragging <- FALSE
                daClearRectangle(w)       # tidy up
                return(FALSE)
            })
            
            ## Right mouse menu -- some means to prevent
            if(is.null(theArgs$no_popup)) {
              l <- list()
              l$copyAction <- gaction("Copy", "Copy current graph to clipboard", icon="copy",
                                      handler=function(h, ...) copyToClipboard(obj))
              l$printAction <- gaction("Save", "Save current graph", icon="save",
                                       handler=function(h,...) {
                                         fname <- gfile(gettext("Filename to save to"), type="save")
                                         if(nchar(fname)) {
                                           if(!file.exists(fname) || gconfirm(gettext("Overwrite file?")))
                                             svalue(obj) <- fname
                                         }
                                       })
              
              .add3rdmousepopupmenu(obj, toolkit, l)
            }
              
              
            ## Add to container if requested
            if (!is.null(container)) {
              if(is.logical(container) && container == TRUE)
                container = gwindow()
              add(container, obj, ...)
            }

            gSignalConnect(da, "realize", function(...) {
              gdkWindowProcessAllUpdates()
              while (gtkEventsPending()) gtkMainIterationDo(blocking=FALSE)
            })
            
            return(obj)
          })

as.gWidgetsRGtk2.GtkDrawingArea <- function(widget,...) {
  obj <- new("gGraphicsRGtk",block=widget, widget=widget,
             toolkit=guiToolkit("RGtk2"))
  return(obj)
}


as.gGd = function(obj) {
  if(inherits(obj,"GtkDrawingArea")) {
    newobj = list(ref = obj, device = obj$GetData("device"))
    class(newobj) <- c("gGd", "gComponent")
    return(newobj)
  } else {
    cat(gettext("conversion failed\n"))
    return(obj)
  }
}


### methods

## adding to a group is funny, we intercept here
setMethod(".add",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gGroupRGtk", value="gGraphicsRGtk"),
          function(obj, toolkit, value, ...) {
            getWidget(obj)$PackStart(value@block, TRUE, TRUE, 0) # expand to fill if TRUE
          })



## raise this device
setReplaceMethod(".visible",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gGraphicsRGtk"),
                 function(obj, toolkit, ..., value) {
                   if(is.logical(value) == TRUE) {
                     da <- obj@widget
                     devnum <- da$GetData(".devnum")
                     if(!is.null(devnum))
                       dev.set(devnum)
                   }
                   return(obj)
                 })

## save Current Page
## This uses GTK -- not R to save.
## need to have window fully shown
setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gGraphicsRGtk"),
                 function(obj, toolkit, index=NULL,  ..., value) {
                   if(length(value) > 1) {
                     file = value$file
                     extension  = value$extension
                   } else {
                     file = value;
                     extension <- strsplit(file, ".", fixed = TRUE)[[1L]]
                     if(n <- length(extension)) {
                       extension <- extension[n]
                     } else {
                       cat(gettext("No file extension given"))
                       return()
                     }
                   }
                   ## check that filename is okay
                   if(!is.null(file) && !is.null(extension)) {
                     tmp = unlist(strsplit(file,"\\."))
                     if(tmp[length(tmp)] != extension) {
                       filename = Paste(file,".",extension)
                     } else {
                       filename = file
                     }
                   } else {
                     return()
                   }

                   da <- getWidget(obj)
                   wh <- daGetWidthHeight(da)
                   da.w <- wh[1]
                   da.h <- wh[2]
                   pixbuf <- gdkPixbufGetFromDrawable(src=da$window, src.x=0, src.y=0,
                                                   dest.x=0, dest.y=0, width=da.w, height=da.h)


                   out <- try(pixbuf$Save(filename = filename,type=extension), silent=TRUE)
                   if(inherits(out, "try-error")) {
                     galert(sprintf("Error in saving: %s", out), parent=obj)
                   }
                   
                   ##   switch(extension,
                   ##          "ps" = dev.copy2eps.hack(file=filename,
                   ##            onefile=onefile, horizontal=horizontal,
                   ##            width=width, height = height),
                   ##          "eps" = dev.print.hack(postscript,file=filename,
                   ##            onefile=onefile, horizontal=horizontal,
                   ##            width=width, height = height),
                   ##          "pdf" = dev.print.hack(pdf,file=filename,
                   ##            onefile=onefile, horizontal=horizontal,
                   ##            width=width, height = height),
                   ##          "jpg" = dev.print.hack(jpeg,file=filename,
                   ##            onefile=onefile, horizontal=horizontal,
                   ##            width=width, height = height),
                   ##          "jpeg" = dev.print.hack(jpeg,file=filename,width=width,height=height),
                   ##          "png" = dev.print.hack(png,file=filename,width=width,height=height),
                   ##          cat("***\n Don't know this extension:", type,"\n\n")
                   ##          )
                   
                   return(obj)
                 })


### handlers
## add this expose event for graph
setMethod(".addhandlerexpose",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gGraphicsRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            addhandler(obj,"expose-event",handler,action,...)
          })

## applies a handler to the mouse click. The handler gets extra
## argument h$x, h$y passed into it. These are in [0,1] coordinates
setMethod(".addhandlerclicked",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gGraphicsRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            ## handler has $obj for obj clicked on, $x, $y, $action


            f = function(h,w,e,...) {
              if(!isFirstMouseClick(e))
                return(FALSE)

              ## changes to allocation storage with newer RGtk2
              xclick = e$GetX()
              yclick = e$GetY()
              da <- getWidget(obj)
              wh <- daGetWidthHeight(da)
              width <- wh[1]
              height <- wh[2]

              
              x = xclick/width
              y = (height - yclick)/height

              ## put into usr coordinates
              h$x <- grconvertX(x, from="ndc", to="user")
              h$y <- grconvertY(y, from="ndc", to="user")
              
              handler(h, w, e, ...)
              return(FALSE)
            }
            
            id = addhandler(obj,signal = "button-press-event",handler=f, action=action)
            invisible(id)
          })


##' Changed handler is called after rubber band selection is updated
##'
##' Just click and drag out a rubber band
##' The "h" list has components
##' h$x for the x values in user coordinates
##' h$y for the y values in user coordinates
##' These can be converted as in grconvertX(h$x, from="ndc", to="user")
setMethod(".addhandlerchanged",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gGraphicsRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            da <- getWidget(obj)
            ID <- gSignalConnect(da, "button-release-event", function(w, e) {
              if(!isFirstMouseClick(e))
                return(FALSE)
              coords <- drawableToNDC(w)
              h <- list(obj=obj,
                        action=action,
                        x=grconvertX(coords$x, from="ndc", to="user"),
                        y=grconvertY(coords$y, from="ndc", to="user"))
              handler(h, ...)
              return(FALSE)             # propagate
            })
          })

##' Draw a rectangle for rubber banding
daDrawRectangle <- function(da,  x0, x, y0, y) {

  x <- c(x0, x); y <- c(y0, y)
  x0 <- min(x); x <- max(x)
  y0 <- min(y); y <- max(y)
  
  allocation = da$allocation ## GetAllocation()
  da.w <- allocation$allocation$width
  da.h <- allocation$allocation$height 

#  da.w <- da$getAllocation()$width
#  da.h <- da$getAllocation()$height

  ## background style
  gcb <- gdkGCNew(da$window)
  gcb$copy(da["style"]$blackGc)
  gcb$setRgbFgColor(gdkColorParse("gray50")$color)
  gcb$setLineAttributes(line.width=1, line.style=GdkLineStyle["solid"],
                        cap.style=GdkCapStyle["butt"], join.style=GdkJoinStyle["miter"])
  ## foreground style
  gc <- gdkGCNew(da$window)
  gc$copy(da["style"]$blackGc)
  gc$setRgbFgColor(gdkColorParse("black")$color)
  gc$setRgbBgColor(gdkColorParse("gray50")$color)
  gc$setLineAttributes(line.width=1, line.style=GdkLineStyle["double-dash"],
                       cap.style=GdkCapStyle["butt"], join.style=GdkJoinStyle["miter"])
  gc$setDashes(c(8, 4))

  ## the entire rectangle to clear
  rect <- as.GdkRectangle(c(x=0, y=0, width=da.w, height=da.h))
  da$setData("lastRect", rect)

  for (i in 1:2) {
    ## draw in background color first
    tmp.gc <- if (i == 1) gcb else gc
    gdkDrawRectangle(da$window, gc=tmp.gc, filled=FALSE, x=x0, y=y0, width=x-x0, height=y-y0)
  }

  gdkWindowProcessAllUpdates()
  while (gtkEventsPending()) gtkMainIterationDo(blocking=FALSE)
   

}


##' find width and height from allocation, which surprisingly seems to change from time to time
daGetWidthHeight <- function(da) {
  allocation <- da$getAllocation()
  ## now, do we have width, height?
  if("width" %in% names(allocation)) {
    return(c(width=allocation$width, height=allocation$height))
  } else if("allocation" %in% names(allocation)) {
    return(c(width=allocation$allocation$width, height=allocation$allocation$height))
  } else {
    stop("Can't get width.height allocation?")
  }
}

##' clear all rectangles that came from rubber banding
daClearRectangle <- function(da) {

  last <- da$getData("lastRect")
  if(!is.null(last)) 
    da$window$invalidateRect(last, FALSE)
  
  gdkWindowProcessAllUpdates()
  while (gtkEventsPending()) gtkMainIterationDo(blocking=FALSE)
}

##' convert rectangle on drawable into NDC coordinates
drawableToNDC <- function(da) {
  ## convert to normalized device coordinates
  e <- da$getData("env")
  x.pixel <- sort(c(e$x0, e$x))
  y.pixel <- sort(c(e$y0, e$y))

  wh <- daGetWidthHeight(da)
  da.w <- wh[1]
  da.h <- wh[2]

  ndc <- list(x=x.pixel/da.w, y= 1- rev(y.pixel/da.h))
  return(ndc)
}
          


          
##' copy graphic to clipboard for cut-and-paste
##'
##' I can't seem to bind this to ctrl-c (or some such), as I can't get key-press-event
##' to work on tihs widget. Here as an example:
##' http://ruby-gnome2.sourceforge.jp/hiki.cgi?tut-gtk2-agtkw-draww
##' da['can-focus'] <- TRUE; da$addEvents(GdkEventMask["key-press-mask"]) were tried
## @param da either the drawable object (from a callback say) or the ggraphics object.
copyToClipboard <- function(da) {
  if(!is(da, "GtkDrawingArea"))
    da <- getWidget(da)                 # ggraphics object
  da.w <- da$getAllocation()$width
  da.h <- da$getAllocation()$height
  buf <- gdkPixbufGetFromDrawable(src=da$window, src.x=0, src.y=0,
                                  dest.x=0, dest.y=0, width=da.w, height=da.h)
  gtkClipboardGet("CLIPBOARD")$setImage(buf)
}

##################################################
##
## dev.print and dev.copy2eps have a test on the device that needs Cairo added to it
devPrintHack = function (device = postscript, ...) 
{
  current.device <- dev.cur()
  nm <- names(current.device)[1]
  if (nm == "null device") 
    stop("no device to print from")
  if (!(nm %in% c("Cairo", "X11", "GTK", "gnome", "windows", "quartz"))) 
    stop("can only print from screen device")
  oc <- match.call()
  print(oc)
  oc[[1]] <- as.name("dev.copy")
  oc$device <- device
  din <- par("din")
  w <- din[1]
  h <- din[2]
  if (missing(device)) {
    if (is.null(oc$file)) 
      oc$file <- ""
    hz0 <- oc$horizontal
    hz <- if (is.null(hz0)) 
      ps.options()$horizontal
    else eval.parent(hz0)
    paper <- oc$paper
    if (is.null(paper)) 
      paper <- ps.options()$paper
    if (paper == "default") 
      paper <- getOption("papersize")
    paper <- tolower(paper)
    switch(paper, a4 = {
      wp <- 8.27
      hp <- 11.69
    }, legal = {
      wp <- 8.5
      hp <- 14
    }, executive = {
      wp <- 7.25
      hp <- 10.5
    }, {
      wp <- 8.5
      hp <- 11
    })
    wp <- wp - 0.5
    hp <- hp - 0.5
    if (!hz && is.null(hz0) && h < wp && wp < w && w < hp) {
      hz <- TRUE
    }
    else if (hz && is.null(hz0) && w < wp && wp < h && h < 
             hp) {
      hz <- FALSE
    }
    else {
      h0 <- ifelse(hz, wp, hp)
      if (h > h0) {
        w <- w * h0/h
        h <- h0
      }
      w0 <- ifelse(hz, hp, wp)
      if (w > w0) {
        h <- h * w0/w
        w <- w0
      }
    }
    if (is.null(oc$pointsize)) {
      pt <- ps.options()$pointsize
      oc$pointsize <- pt * w/din[1]
    }
    if (is.null(hz0)) 
      oc$horizontal <- hz
    if (is.null(oc$width)) 
      oc$width <- w
    if (is.null(oc$height)) 
      oc$height <- h
  }
  else {
    devname <- deparse(substitute(device))
    if (devname %in% c("png", "jpeg", "bmp") && is.null(oc$width) && 
        is.null(oc$height)) 
      warning("need to specify one of 'width' and 'height'")
    if (is.null(oc$width)) 
      oc$width <- if (!is.null(oc$height)) 
        w/h * eval.parent(oc$height)
      else w
    if (is.null(oc$height)) 
      oc$height <- if (!is.null(oc$width)) 
        h/w * eval.parent(oc$width)
      else h
  }
  dev.off(eval.parent(oc))
  dev.set(current.device)
}

dev.copy2eps.hack = function (...) 
{
  current.device <- dev.cur()
  nm <- names(current.device)[1]
  if (nm == "null device") 
    stop("no device to print from")
  if (!(nm %in% c("Cairo","X11", "GTK", "gnome", "windows", "quartz"))) 
    stop("can only print from screen device")
  oc <- match.call()
  
  
  oc[[1]] <- as.name("dev.copy")
  oc$device <- postscript
  oc$onefile <- FALSE
  oc$horizontal <- FALSE
  if (is.null(oc$paper)) 
    oc$paper <- "special"
  din <- par("din")
  w <- din[1]
  h <- din[2]
  if (is.null(oc$width)) 
    oc$width <- if (!is.null(oc$height)) 
      w/h * eval.parent(oc$height)
    else w
  if (is.null(oc$height)) 
    oc$height <- if (!is.null(oc$width)) 
      h/w * eval.parent(oc$width)
    else h
  if (is.null(oc$file)) 
    oc$file <- "Rplot.eps"
  dev.off(eval.parent(oc))
  dev.set(current.device)
} 

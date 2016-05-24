setClass("gWindowRGtk",
         contains="gContainerRGtk",
         prototype=prototype(new("gContainerRGtk"))
         )

## constructor
setMethod(".gwindow",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   title="Window", visible=TRUE,
                   width = NULL, height = NULL,
                   parent = NULL,
                   handler=NULL, action = NULL,
                   ...
                   ) {

            force(toolkit)
            
            window <- gtkWindowNew("toplevel", show = FALSE)
            window$SetTitle(title)

            ## set default size give 400 x 280 default
            if(is.null(width)) width <- 400
            if(is.null(height)) height = .7*width
            window$SetDefaultSize(width, height)

            ## set location -- renamed to parent
            location <- parent
            if(!is.null(location)) {
              if(inherits(location,"guiContainer") ||
                 inherits(location,"guiComponent")) {
                ## a gWidget.
                widget <- getToolkitWidget(location)
                if(!inherits(widget,"GtkWindow"))
                  widget <- getGtkWindow(widget)
                window$SetTransientFor(widget)
                window$SetPosition(GtkWindowPosition["center-on-parent"])
                window$SetDestroyWithParent(TRUE)
                ## windows fixes
                window$setSkipTaskbarHint(TRUE)
                window$setSkipPagerHint(TRUE)
              } else {
                ## check that location is a numeric pair
                if(is.numeric(location) && length(location) >= 2) {
                  location <- as.integer(location)
                  window$Move(location[1],location[2])
                }
              }
            }

            ## make object
            obj <- as.gWidgetsRGtk2(window)

            
            if (!is.null(handler)) {
              ## handler can't refer to h$obj, as it is already <invalid>
              ## by the time it gets here.
              id <- addhandlerdestroy(obj, handler=handler, action=action)
            }

            if(visible)
              window$Show()

            return(obj)
          })

as.gWidgetsRGtk2.GtkWindow <- function(widget,...) {
  window <- widget
  obj <- new("gWindowRGtk",block=window, widget=window,
    toolkit=guiToolkit("RGtk2"))

  if(!is.null(tag(obj,"menubargroup"))) {
    ## already a gwindow. Move on
    return(obj)
  }

  ## may or may not have child. 
  child <- window$GetChild()
  
  ## if there, save child, then put into contentPane
  if(!is.null(child)) 
    window$Remove(child)   # put into cpg
  
  (mbg <- ggroup(spacing=0)); svalue(mbg) <- 0
  (tbg <- ggroup(spacing=0)); svalue(tbg) <- 0
  (ibg <- ggroup(spacing=0, horizontal=FALSE)); svalue(ibg) <- 0
  (cpg <- ggroup(spacing=0)); svalue(cpg) <- 0
  (sbg <- ggroup(spacing=0)); svalue(sbg) <- 0
  
  tag(obj,"menubargroup") <- mbg
  tag(obj,"toolbargroup") <- tbg
  tag(obj, "infobargroup") <- ibg
  tag(obj,"contentPane") <- cpg
  tag(obj,"statusbargroup") <- sbg
  
  tbl <- gtkTable(rows=4, columns=1, homogeneous=FALSE)
  tag(obj,"table") <- tbl
  tbl$SetColSpacings(0)
  tbl$SetRowSpacings(0)
  
  tbl$Attach(getBlock(mbg), 0,1,0,1, yoptions = c("fill"))
  tbl$Attach(getBlock(tbg), 0,1,1,2, yoptions = c("fill"))
  tbl$Attach(getBlock(ibg), 0,1,2,3, xoptions=c("shrink", "fill"), yoptions = c("shrink"))
  tbl$AttachDefaults(getBlock(cpg), 0,1,3,4)
  ## size grip issue if no statusbar
  tmp <- getBlock(cpg);  tmp['border-width'] = 13

  tbl$Attach(getBlock(sbg), 0,1,5,6, yoptions = c("fill"))
  
  window$Add(tbl)

  ## give back child if there
  if(!is.null(child)) add(cpg, child, expand=TRUE)
  
  return(obj)
  
}

##################################################
## Methods

## Old method, when gwindow did not have a ggroup packed in.
## setMethod(".add",
##           signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWindowRGtk", value="gWidgetRGtk"),
##           function(obj, toolkit, value, ...) {
##             getWidget(obj)$Add(value)
##           })


## methods
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWindowRGtk"),
          function(obj, toolkit, index=NULL, drop=NULL, ..) {
            getWidget(obj)$GetTitle()
          })

setMethod(".svalue<-",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWindowRGtk"),
          function(obj, toolkit, index=NULL,..., value) {
            ## set the title
            getWidget(obj)$SetTitle(value)
            return(obj)
          })

## no visible() method
setMethod(".visible<-",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWindowRGtk"),
          function(obj, toolkit, ..., value) {
            value = as.logical(value)
            if(value == TRUE)
              getWidget(obj)$Show()
            else
              getWidget(obj)$Hide()

            return(obj)
          })

##' focus will raise window
##'
##' @param object gwindow object
##' @param toolkit name of toolkit
##' @param ... ignored
##' @return NULL called for side effect of raising window
setMethod(".focus",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWindowRGtk"),
          function(obj, toolkit, ...) {
            w <- getBlock(obj)
            w$present()
          })

setReplaceMethod(".focus",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWindowRGtk"),
                 function(obj, toolkit, ..., value) {
                   if(as.logical(value)) {
                     w <- getBlock(obj)
                     w$present()
                   }
                   return(obj)
                 })
                
setMethod(".size",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWindowRGtk"),
          function(obj, toolkit, ...) {
            theSize = getWidget(obj)$GetSize()
            return(unlist(theSize[2:3]))
          })

setMethod(".update",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWindowRGtk"),
          function(object, toolkit, ...) {
            w <- getWidget(object)
            w$setSizeRequest(-1, -1)
            invisible()
          })

## Add and delete. Special methods for [menu|tool|status]bars
##  add
setMethod(".add",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWindowRGtk", value="gWidgetRGtk"),
          function(obj, toolkit, value, ...) {
            ## should fix expand=TRUE here
#            .add(obj, toolkit, getBlock(value),...)
            add(tag(obj,"contentPane"), value, expand=TRUE, fill="both") # no ...
          })

setMethod(".add",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWindowRGtk", value="RGtkObject"),
          function(obj, toolkit, value, ...) {
            ## should fix expand=TRUE here
            theArgs=list(...)
            theArgs$expand=TRUE
            gp <- tag(obj,"contentPane")
            do.call("add",list(obj=gp,value=value,theArgs))
#            add(tag(obj,"contentPane"), value, ...)
          })

## menubar
setMethod(".add",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWindowRGtk", value="gMenuRGtk"),
          function(obj, toolkit, value, ...) {
            add(tag(obj,"menubargroup"), value, expand=TRUE)
          })
## toolbar
setMethod(".add",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWindowRGtk", value="gToolbarRGtk"),
          function(obj, toolkit, value, ...) {
            add(tag(obj,"toolbargroup"), value, expand=TRUE)
          })
## statusbar
setMethod(".add",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWindowRGtk", value="gStatusbarRGtk"),
          function(obj, toolkit, value, ...) {
            add(tag(obj,"statusbargroup"), value, expand=TRUE)
            tmp <- getBlock(tag(obj, "contentPane"))
            tmp['border-width'] <- 0
          })

## delete
setMethod(".delete",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWindowRGtk", widget="gWidgetRGtk"),
          function(obj, toolkit, widget, ...) {
            delete(tag(obj,"contentPane"), widget, ...)
          })
## menubar
setMethod(".delete",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWindowRGtk", widget="gMenuRGtk"),
          function(obj, toolkit, widget, ...) {
            delete(tag(obj,"menubargroup"), widget, ...)
          })
## toolbar
setMethod(".delete",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWindowRGtk", widget="gToolbarRGtk"),
          function(obj, toolkit, widget, ...) {
            delete(tag(obj,"toolbargroup"), widget, ...)
          })
## statusbar
setMethod(".delete",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWindowRGtk", widget="gStatusbarRGtk"),
          function(obj, toolkit, widget, ...) {
            delete(tag(obj,"statusbargroup"), widget, ...)
            widget <- getWidget(obj)
            tmp <- getBlock(tag(obj, "contentPane"))
            tmp['border-width'] <- 13
          })


## dispatches
setMethod(".dispose",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWindowRGtk"),
          function(obj, toolkit, ...) {
            obj@widget$Destroy()
          })

##################################################
## handlers
## THis intercepts the windowmanager delete-event, destroy does not
setMethod(".addhandlerunrealize",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWindowRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            theArgs = list(...)
            gtktry(connectSignal(obj@widget,
                          signal="delete-event",
                          f = function(...) {
                            val = handler(...)
                            if(is.logical(val))
                              return(val)
                            else
                              return(FALSE) # do delete
                          },
                          data=list(obj=if(!is.null(theArgs$actualobj))
                            theArgs$actualobj
                          else
                            obj, action=action,...),
                          user.data.first = TRUE,
                          
                          after=FALSE),
                silent=TRUE)
          })

setMethod(".addhandlerdestroy",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWindowRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            .addHandler(obj, toolkit, signal="destroy", handler, action, ...)
          })

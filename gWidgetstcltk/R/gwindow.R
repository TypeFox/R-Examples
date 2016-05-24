## constructor
setMethod(".gwindow",
          signature(toolkit="guiWidgetsToolkittcltk"),
          function(toolkit,
                   title="Window", visible=TRUE,
                   width = NULL, height = NULL, parent=NULL,
                   handler=NULL, action = NULL,
                   ...
                   ) {

            force(toolkit)

            ## don't draw until asked
            tclServiceMode(FALSE)

            win <- tktoplevel()
            tktitle(win) <- title
            tkwm.state(win,"withdrawn") # was at beginneing
            
            ## enable autoresizing
            tkwm.geometry(win,"")


            
            ## how to set location???
            location <- parent          # renamed
            if(!is.null(location)) {
              if(is(location,"guiWidget") ||
                 is(location, "gWindowtcltk") ||
                 is(location, "tkwin")) {
                location <- getToolkitWidget(location)
                curgeo <- tclvalue(tkwm.geometry(location))
                ## widthXheight+xpos+ypos
                pos <- unlist(strsplit(curgeo, "\\+"))
                sz <- unlist(strsplit(pos[1],"x"))
                xpos = as.numeric(pos[2]); ypos=as.numeric(pos[3])
                tkwm.geometry(win,paste("+",xpos+30,"+",ypos+30,sep="")) # shift

                tkwm.transient(win, location) # set transient
                tkbind(location,"<Destroy>",function(...) tkdestroy(win))
              } else if(is.numeric(location) && length(location) == 2) {
                tkwm.geometry(win, location[1], location[2])
              }

            }

            ## pack a frame inside for theme issues:
            ## tkdocs.com:

##             Strictly speaking, we could just put the other parts of
##             our interface directly into the main root window, without
##             the intervening content frame. However, the main window
##             isn't itself part of the "themed" widgets, so its background color wouldn't
##             match the themed widgets we will put inside it. Using a
##             "themed" frame widget to hold the content ensures that the
##             background is correct.
            
            ## pack in frame for adding to
            contentPane <- ttkframe(win, padding=c(3,3,12,12))
            tkgrid(contentPane, row=1, column = 0, sticky="nwes")
            tkgrid.columnconfigure(win, 0, weight = 1)
            tkgrid.rowconfigure(win, 1, weight = 1)

            ## pack in toolbar
            tb <- ttkframe(win)
            tkgrid(tb, row=0, column = 0, sticky = "nswe")

            ## pack in statusbar
            sb <- ttkframe(win)
            tkconfigure(sb, borderwidth = 1, relief="sunken")
            tkgrid(sb, row=2, column = 0, sticky="we")
            
            ## debugging code
            ## just to see the frame
            ## tkconfigure(contentPane, borderwidth=4, relief="solid")
            ## tkconfigure(tb, borderwidth=4, relief="solid")

            ## size the frame object
            ## set default size? only minsize here
            if(!is.null(width)) {
              if(is.null(height)) height = .7*width
              tkconfigure(contentPane, width=as.integer(width), height=as.integer(height))
              tkgrid.propagate(contentPane,FALSE) ## make frame size carry forward
            }

            
            obj <- new("gWindowtcltk",block=win, widget=contentPane, toolkit=toolkit,
              ID=getNewID(),e=new.env())

            obj@e$parentContainer <- NULL
            tag(obj,"tb") <- tb
            tag(obj,"sb") <- sb
            
            
            if (!is.null(handler)) {
              id <- addhandlerdestroy(obj, handler=handler, action=action)
            }

            tclServiceMode(TRUE)

            if(visible) {
              tkwm.state(win,"normal")
            }

            return(obj)
          })
##################################################
## Methods 
## getToolkitWidget returns window -- not frame


## general add
setMethod(".add",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWindowtcltk", value="gWidgettcltk"),
          function(obj, toolkit, value, ...) {

            ## add parent, children
            childComponents <- obj@e$childComponents
            if(is.null(childComponents))
              childComponents <- list()
            obj@e$childComponents <- c(childComponents, value)
            value@e$parentContainer <- obj

            
            ## pack into frame
            tkpack(getBlock(value),
                   expand=TRUE, fill="both")
            return(TRUE)

            ## --- IGNORED --
            ## adding widget to window means pack
            theArgs = list(...)
            packArgs = list(getBlock(value))
            if(!is.null(theArgs$expand) && theArgs$expand) {
             packArgs$expand=TRUE
              packArgs$fill = "both"
              packArgs$side="top"
            } else {
              packArgs$side="top"
            }
            ## override with anchor argument
            if(!is.null(theArgs$anchor)) {
              an = theArgs$anchor
              if(an[1] == 1)
                packArgs$side = "right"
              else if(an[1] == -1)
                packArgs$side = "left"
              else if(an[2] == 1)
                packArgs$side = "top"
              else
                packArgs$side = "bottom"
            }

            
            #do.call("tkpack", packArgs)
            packArgs$side <- NULL       # clera out for test
            do.call("tkgrid", packArgs)

            
          })

## return window -- not frame
setMethod(".getToolkitWidget",
          signature(obj="gWindowtcltk", toolkit="guiWidgetsToolkittcltk"),
          function(obj, toolkit) obj@block)

## add toolbar, menubar, statusbar
## menubar -- in gmenu

## toolbar
setMethod(".add",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWindowtcltk", value="gToolbartcltk"),
          function(obj, toolkit, value, ...) {

            tkpack(getBlock(value), anchor="w")

            ##             ## put before all others.
##             ## get children, check then put in. XXX
##             ## XXX  -- not working
##             g <- getWidget(obj)
##             slaves <- unlist(strsplit(tclvalue(tkpack("slaves",g))," "))
##             args <- list(getBlock(value),
##                          side="top",anchor="w",expand=FALSE, fill="x")
##             if(length(slaves))
##               args$before = slaves[1]
##             do.call("tkpack",args)

##            tag(obj,"toolbar") <- getBlock(value)
          })
## statusbar
setMethod(".add",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWindowtcltk", value="gStatusbartcltk"),
          function(obj, toolkit, value, ...) {

            tkpack(getBlock(value), anchor="w")

            ##            ## put after all others
##            ## XXX Get children, put last -- NOT WORKING!!
##             g = getWidget(obj)
##             slaves = unlist(strsplit(tclvalue(tkpack("slaves",g))," "))
##             args <- list(getBlock(value),
##                          side="top",anchor="w",expand=FALSE, fill="x")
##             if(length(slaves))
##               args$after <- slaves[length(slaves)]
##             do.call("tkpack",args)

##             tag(obj,"statusbar") <- getBlock(value)
            
          })


## methods

## svalue refers to title
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWindowtcltk"),
          function(obj, toolkit, index=NULL, drop=NULL, ..) {
            ## return title
            val <- tcl("wm","title",getBlock(obj))
            tclvalue(val)
          })

setMethod(".svalue<-",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWindowtcltk"),
          function(obj, toolkit, index=NULL,..., value) {
            ## set the title
            tcl("wm","title",getBlock(obj), as.character(value))
            return(obj)
          })


setMethod(".size", 
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWindowtcltk"),
          function(obj, toolkit, ...) {
            widget <- getBlock(obj)
            width <- tclvalue(tkwinfo("width",widget))
            height <- tclvalue(tkwinfo("height",widget))
            return(as.numeric(c(width=width, height=height)))
          })

setReplaceMethod(".size",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWindowtcltk"),
          function(obj, toolkit, ...,value) {
            tkwm.minsize(getBlock(obj), value[1], value[2])
            return(obj)
          })

setMethod(".dispose",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWindowtcltk"),
          function(obj, toolkit, ...) {
            tcl("after",5,function() {
              tkdestroy(getBlock(obj))
            })
          })

setReplaceMethod(".visible",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWindowtcltk"),
          function(obj, toolkit, ...,value) {
            if(as.logical(value)) {
              tkwm.state(obj@block,"normal")
            } else {
              tkwm.state(obj@block,"withdrawn")
            }
            return(obj)
            })
          
##' update will cause window to resize to natural size
##'
##' @param object gwindow object
##' @param toolkit name of toolkit
##' @param ... ignored
##' @return NULL
setMethod(".update",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWindowtcltk"),
          function(object, toolkit, ...) {
            w <- getBlock(object)
            tkwm.geometry(w, "")
            invisible()
          })


##' focus will raise window
##'
##' @param object gwindow object
##' @param toolkit name of toolkit
##' @param ... ignored
##' @return NULL called for side effect of raising window
setMethod(".focus",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWindowtcltk"),
          function(obj, toolkit, ...) {
            w <- getBlock(obj)
            tkraise(w)
          })

setReplaceMethod(".focus",
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gWindowtcltk"),
                 function(obj, toolkit, ..., value) {
                   if(as.logical(value)) {
                     w <- getBlock(obj)
                     tkraise(w)
                   }
                   return(obj)
                 })
                
##################################################
## handlers


setMethod(".addhandlerunrealize",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWindowtcltk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            win <- getBlock(obj)
            h <- list(obj = obj, action=action,...)
            tkwm.protocol(win, "WM_DELETE_WINDOW",
                          function(...) {
                            val <- handler(h,...)
                            ## FALSE -- destroy, TRUE -- keep
                            if(is.null(val)  || !is.logical(val) || !val)
                              tkdestroy(win) ## revers
                          })
          })

## no ID
setMethod(".addhandlerdestroy",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWindowtcltk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            ## need to tkbind explicitly here
            f <- function() {
              h <- list(obj=obj, action=action)
              handler(h)
            }
            tkbind(getWidget(obj), "<Destroy>", f)
          })

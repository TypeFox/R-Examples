setClass("gPanedgrouptcltk",
         contains="gContainertcltk",
         prototype=prototype(new("gContainertcltk"))
         )

## TODO: method obj[1 or 2 ] <- replacewidget
setMethod(".gpanedgroup",
          signature(toolkit="guiWidgetsToolkittcltk"),
          function(toolkit,
                   widget1, widget2, horizontal=TRUE, container=NULL, ...) {
            ## add a paned group
            
            force(toolkit)
            
            
            if(is.null(container)) {
              message(gettext("No NULL containers in tcltk. Creating a new window\n"))
              container=gwindow()
            } else if(is.logical(container) && container) {
              container = gwindow()
            }
            
            if(!is(container,"guiWidget")) {
              container = gwindow()
            }
            
            ## process args
            if(horizontal)
              orient = "horizontal"
            else
              orient = "vertical"
            
            
            tt <- getWidget(container)
            ##            pg <- tkwidget(tt,"panedwindow", orient=orient)
            pg <- ttkpanedwindow(tt, orient=orient)
            tkpack(pg, expand=TRUE, fill="both")
            
            
            ## make object -- note block is pg so that add works correctly
            ## as it calls getBlock(container)
            obj = new("gPanedgrouptcltk", block=pg, widget=pg,
              toolkit=toolkit,ID=getNewID(), e = new.env())
            
            tag(obj,"horizontal") <- horizontal
            
            if(!missing(widget1) || !missing(widget2)) {
              gwCat(gettext("In tcltk, you use the gpanedgroup object in the container argument of a constructor\n"))
            }
            
            return(obj)
          })


## add -- use this rather than at construction time
setMethod(".add",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gPanedgrouptcltk", value="gWidgettcltk"),
          function(obj, toolkit, value, ...) {
            ## add parent, children
            childComponents <- obj@e$childComponents
            if(is.null(childComponents))
              childComponents <- list()
            obj@e$childComponents <- c(childComponents, value)
            value@e$parentContainer <- obj

            theArgs = list(...)
#            argList = list(getWidget(obj),"add",getBlock(value))
            argList = list(getWidget(obj),"insert","end",getBlock(value))

            ## args to position
            sticky = "n"
            if(!is.null(theArgs$anchor)) {
              sticky = xyToAnchor(theArgs$anchor)
            }
            if(!is.null(theArgs$expand) && theArgs$expand) {
              if(tag(obj,"horizontal"))
                sticky = "news"
              else
                sticky = "news"
            }
            argList$sticky = sticky

            ## for ttk
            argList$sticky <- NULL
            
            do.call("tcl", argList) ## tcl(tt,"add",widget,...)
            
          })

## delete means we can readd -- in this case we actually dispose, as
## the widget doesn't get added back?
setMethod(".delete",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gPanedgrouptcltk",
                    widget="gWidgettcltk"),
          function(obj, toolkit, widget, ...) {
            ## call forget

            tcl(getWidget(obj),"forget",getBlock(widget))
          })



## svalue refers to sash position between 0 and 1
## sashpos
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gPanedgrouptcltk"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {
            sashpos <- as.numeric(tclvalue(tcl(getWidget(obj),"sashpos",0)))
            theSize <- size(obj)
            
            if(tag(obj,"horizontal"))
              return(sashpos/theSize[1])
            else
              return(sashpos/theSize[2])
          })

## svalue sets position
setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gPanedgrouptcltk"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   if(0 <= value && value <= 1) {
                     theSize <- size(obj)
                     if(tag(obj,"horizontal"))
                       pos <- floor(value *  theSize[1])
                     else
                       pos <- floor(value *  theSize[2])

                     tcl(getWidget(obj),"sashpos", 0, as.integer(pos))
                   }
                   return(obj)
                 })

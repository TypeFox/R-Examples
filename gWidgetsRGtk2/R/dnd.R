
                                        # functions to handle DND

## helper, like rawToChar. From R.oo/R/ASCII.R
# Alternatively one can do like this. Idea by Peter Dalgaard,
# Dept. of Biostatistics, University of Copenhagen, Denmark.

## ASCII <- c("\000", sapply(1:255, function(i) parse(text=paste("\"\\",
##                                                      structure(i,class="octmode"), "\"", sep=""))[[1]]) );

## intToChar = function(i, ...) {
##   ASCII[i %% 256 + 1];
## }

## Brian Ripley says the above will fail as of 2.8.0 version of R.
## So we try this instead

intToChar <- function(x) rawToChar(as.raw(x), multiple = TRUE)

## A little buggy right now: drop target had drag-data-received called 2 times
## action argument in addhandler isn't handled properly
## a gross hack to allow objects to be dropped.

TARGET.TYPE.TEXT   = 80                 # 
TARGET.TYPE.PIXMAP = 81                 # NOT IMPLEMENTED
TARGET.TYPE.OBJECT = 82

gWidgetTargetTypes = list(
  text = gtkTargetEntry("text/plain", 0, TARGET.TYPE.TEXT),
  pixmap = gtkTargetEntry("image/x-pixmap", 0, TARGET.TYPE.PIXMAP),
  object = gtkTargetEntry("text/plain", 0, TARGET.TYPE.OBJECT)
 )
  

## Part of gross hack to allow objects to be dropped
## hide this list for storing drop information. This is typically a pointer to an RGtkObject
## .gWidgetDropTargetList <- list()
.gWidgetDropTargetList <- new.env()
.gWidgetDropTargetListKey = ".gWidgetDropTargetListKey" # goes in front


DropList <- setRefClass("DropList",
                        fields=list(
                          l="list"
                          ),
                        methods=list(
                          initialize=function(...) {
                            initFields(l=list())
                            callSuper(...)
                          },
                          set_key=function(key, value) {
                            l[[key]] <<- value
                          },
                          get_key=function(key, remove=TRUE) {
                            l[[key]]
                            if(remove)
                              l[[key]] <<- NULL
                          }))

                          


## function used by RGtkObject and gWidgetRGtk
addDropSource = function(obj, toolkit, targetType="text", handler=NULL, action=NULL, ...) {

  ##  ver = getRGtk2Version()   ## too slow!
  x <- read.dcf(system.file("DESCRIPTION", package="RGtk2"))
  version <- x[1,'Version']
  ver <- strsplit(version,"\\.")[[1]]
  names(ver) <- c("major","minor","mini?")
  
  tmp = gtkDragSourceSet(getWidget(obj),
    if(ver['major'] == "2" && as.numeric(ver['minor']) < 10) {
      GdkModifierType[c("button1-mask","button3-mask")]
    } else {
      GdkModifierType["button1-mask"] | GdkModifierType["button3-mask"]
    },
    list(gWidgetTargetTypes[[targetType]]), #    targets,
    GdkDragAction["copy"])
  
  
  ## uses handler in a closure
  sourceHandler = function(h, widget, context, selection,
    targetType, eventTime) {
    ## what gets set in selection gets passed on to drop target
    if(targetType == TARGET.TYPE.PIXMAP) {
      ## this is untested!
      selection$Set(selection$Target(), 8,
                    paste(svalue(h$obj),collapse="\n"))
    } else if(targetType == TARGET.TYPE.OBJECT) {
      ## this is tricky! we want to store an object, but selection
      ## seemingly only likes to store text. So instead we store the
      ## name of a component in an invisible list we've snuck into the
      ## globalenvironment.
      ## This assumes the object you want to sneak into your
      ## DND is in action argument
      if(!is.null(action)) {
        key = Paste(.gWidgetDropTargetListKey,tempfile())                  # why not?
        ##        .gWidgetDropTargetList[[key]] <<- action
#        tmplst = getFromNamespace(".gWidgetDropTargetList",
#          "gWidgetsRGtk2")
        tmplst <- .gWidgetDropTargetList[["gWidgetsRGtk2"]]

        ## the tmplst is empty!!
#        tmplst[[key]] <- action
 #       .gWidgetDropTargetList[["gWidgetsRGtk2"]] <- tmplst
        .gWidgetDropTargetList[["gWidgetsRGtk2"]] <- list(key=key, action=action)

##        assignInNamespace(".gWidgetDropTargetList", tmplst,
##                          "gWidgetsRGtk2")
        
        selection$SetText(key)
      }
    } else {
      ## it is TEXT type
      if(is.null(handler)) {
        value = svalue(h$obj)
      } else {
        value  = gtktry(handler(h), silent=TRUE)
        if(inherits(value,"try-error")) {
          gwCat(sprintf("Error: handler returns: %s\n",value))
        }
      }
      ## what gets set here is passed to drop target
      selection$SetText(str=value)
    }
    return(TRUE)
            }
  ## this gets drag-data-get signal
  ## action isn't working! (For "object" action is passed in already)
  
  theArgs = list(...)
  
  id = .addHandler(obj,toolkit,"drag-data-get",sourceHandler,actualobj=theArgs$actualobj)#action=action)
  invisible(id)
  
}


setMethod(".adddropsource",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="RGtkObject"),
          function(obj, toolkit, targetType="text",
                   handler=NULL, action=NULL, ...) {
            addDropSource(obj, toolkit, targetType, handler, action, ...)
          })


setMethod(".adddropsource",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit, targetType="text",
                   handler=NULL, action=NULL, ...) {
            addDropSource(obj, toolkit, targetType, handler, action, ...)
          })

## motino
setMethod(".adddropmotion",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit,  handler=NULL, action=NULL, ...) {
            .addHandler(obj,toolkit, signal="drag-motion",handler, action, ...)
          })
setMethod(".adddropmotion",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="RGtkObject"),
          function(obj, toolkit,  handler=NULL, action=NULL, ...) {
            .addHandler(obj,toolkit, signal="drag-motion",handler, action, ...)
          })

## target -- how to add for RGtkObjects?
addDropTarget = function(obj, toolkit, targetType="text", handler=NULL, action=NULL,
  actualobj = NULL,...) {
  ## acutalobj is used by glabel to put target onto evb, use obj for svalue() etc.
  
  gtkDragDestUnset(getWidget(obj))

  tmp = gtkDragDestSet(getWidget(obj),
    c("GTK_DEST_DEFAULT_ALL"),
                                        #targets,
    list(gWidgetTargetTypes[[targetType]]),
    GdkDragAction[c("copy")]
    )
  id = NA
  
  ## gets handler for closure
  drophandler = function(
    h,
    widget,context, x, y, selection,
    targetType, eventTime
    ) {

    ## override
    if(!is.null(h$actualobj))
      h$obj = h$actualobj
    
    ## we would like to filter by target type, but it doesn't
    ## work right for us, so instead we hack this in. If the
    ## text has the key in from do onething, otherwise do the other
    
    ## get dropdata
    if(targetType == TARGET.TYPE.OBJECT ||
       targetType == TARGET.TYPE.TEXT) {
      dropdata = selection$GetText()
      if(is.integer(dropdata)) 
        dropdata = Paste(intToChar(dropdata))
      else
        dropdata = rawToChar(dropdata)
      
      ## is this an actino thingy, or not?
      if(length(grep(Paste("^",.gWidgetDropTargetListKey), dropdata)) > 0) {
        ## Dropdata is key, look up value in .gWidgetDropTargetList

        
        lst <- .gWidgetDropTargetList[["gWidgetsRGtk2"]]

        sourceAction <- lst$action
        ## clear
        .gWidgetDropTargetList[["gWidgetsRGtk2"]] <- list(key="", action=NULL)

        ## XXXX
##         ## It is an action thing. An object was dropped, not a text value
##         sourceAction = .gWidgetDropTargetList[[dropdata]]
##         ##        .gWidgetDropTargetList[[dropdata]] <<- NULL
## ##        tmplst = getFromNamespace(".gWidgetDropTargetList", "gWidgetsRGtk2")
##         tmplst <- .gWidgetDropTargetList[["gWidgetsRGtk2"]]
##         tmplst[[dropdata]] <- NULL
## ##        assignInNamespace(".gWidgetDropTargetList", tmplst, "gWidgetsRGtk2")
##         .gWidgetDropTargetList[["gWidgetsRGtk2"]] <-  tmplst 
        ## what to do with handler?
        if(!is.null(handler)) {

          
          h$dropdata = sourceAction; h$x = x; h$y = y
          out = gtktry(
            handler(h, widget=widget, context=context, x=x, y=y, selection=selection,
                    targetType=targetType,
                    eventTime=eventTime),
            silent=TRUE)
          if(inherits(out,"try-error")) {
            gwCat(sprintf("Error: handler has issue: %s\n",out))
          }
        } else{
          gwCat(gettext("No default handler when action object is passed in\n"))
        }
      } else {
        ## this is text case
                  dropdata = gsub(Paste("^",.gWidgetDropTargetListKey),"", dropdata)
                  ## set drop data into object passed to handlers
                  if(!is.null(handler)) {             # handler = function(h,...)
                    h$dropdata = dropdata; h$x = x; h$y = y
                    out = gtktry(
                      handler(h ,widget=widget, context=context, x=x, y=y,
                              selection=selection,
                              targetType=targetType,
                              eventTime=eventTime),
                      silent=TRUE)
                    if(inherits(out,"try-error")) {
                      gwCat(sprintf("Error: handler has issue: %s\n",out))
                    }
                  } else {
                    svalue(h$obj) <- dropdata
                  }
                }
                return(TRUE)
              } else {
                gwCat(gettext("Nothing defined for this Target type\n"))
              }
            }
            
            ## Why is pixmap stuff not working? --later
            ##    else if(targetType == TARGET.TYPE.PIXMAP) {
            ##      dropdata = selection$GetPixbuf()
            ##      if(!is.null(handler)) {
            ##        h$dropdata = dropdata; h$x = x; h$y = y
            ##        handler(h,widget=widget, context=context, x=x, y=y, selection=selection,
            ##                targetType=targetType,
            ##                eventTime=eventTime)
            ##      } else {
            ##        cat("No default handler for pixbuf data.\n")
            ##      }        
            ##    }
            ##    else {
            ##      ## TARGET.TYPE.TEXT
            ##      dropdata = selection$GetData()
            ##      if(is.integer(dropdata)) 
            ##        dropdata = Paste(intToChar(dropdata))
            ##      else
            ##        dropdata = rawToChar(dropdata)
            
            ##      ## set drop data into object passed to handlers
            ##      if(!is.null(handler)) {             # handler = function(h,...)
            ##        h$dropdata = dropdata; h$x = x; h$y = y
            ##        handler(h,widget=widget, context=context, x=x, y=y, selection=selection,
            ##                targetType=targetType,
            ##                eventTime=eventTime)
            ##      } else {
            ##        svalue(h$obj) <- dropdata
            ##      }
            ##      return(FALSE)
            
            ##    }
            ## }
            
            ## now add drop handler and return id
            id = .addHandler(obj,toolkit,"drag-data-received", drophandler, action=action, actualobj=actualobj)
            invisible(id)
          }

setMethod(".adddroptarget",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit, targetType="text", handler=NULL, action=NULL, ...) {
            addDropTarget(obj, toolkit, targetType, handler, action, ...)
          })
setMethod(".adddroptarget",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="RGtkObject"),
          function(obj, toolkit, targetType="text", handler=NULL, action=NULL, ...) {
            addDropTarget(obj, toolkit, targetType, handler, action, ...)
          })
            

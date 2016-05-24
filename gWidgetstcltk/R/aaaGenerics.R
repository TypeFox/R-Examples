require(methods)
require(digest)
require(tcltk)


MSG = function(...) message("DEBUG: ",...,"\n")
missingMsg = function(x) {
  if(missing(x)) x = "XXX"
  message("This method ",x," needs to be written\n")
}


## toolkit class
## register classes here for toolkits
## Not needed as in gwidgets
## setClass("guiWidgetsToolkittcltk",
##          contains="guiWidgetsToolkit",
##          prototype=prototype(new("guiWidgetsToolkit"))
##          )




##################################################
## put S3 classes from tcltk into S4 classes
## got these from apropos("New") -> try(class(do.call(i,list())))


oldClasses =c("tkwin", "tclVar", "tclObj")
setClass("tcltkObject")
lapply(oldClasses, function(i) {
  setOldClass(i)
  setIs(i,"tcltkObject")
})


setOldClass("try-error")                # for handling try-errors


## a base class which is virtual


##################################################
## A virtual class to hold either RGTK or these guys

## A virtual class for our newly defined objects
## this one contains the ID for the object.
## this may better be done within the NAMESPACE


id.env <- new.env()
id.env[['n']] <- 0L
getNewID = function() {                 # get new one, incremented
  n <- id.env[['n']]
  id.env[['n']] <- n + 1
  return(n+1)
}
         

setClass("gWidgettcltk",
         representation(ID="numeric",
                        e="environment"
                        ),
         )


setClassUnion("guiWidgetORgWidgettcltkORtcltkObject",
              c("guiWidget","gWidgettcltk","tcltkObject"))

## subclss
setClass("gComponenttcltk",
         representation(
                        block="guiWidgetORgWidgettcltkORtcltkObject",
                        widget="guiWidgetORgWidgettcltkORtcltkObject",
                        toolkit="guiWidgetsToolkit"
                        ),
         contains="gWidgettcltk",
         )
setClass("gContainertcltk",
         representation(
                        block="guiWidgetORgWidgettcltkORtcltkObject",
                        widget="guiWidgetORgWidgettcltkORtcltkObject",
                        toolkit="guiWidgetsToolkit"
                   ),
         contains="gWidgettcltk",
         )

setClass("gComponentR5tcltk",
         representation(R5widget="ANY"),
         contains="gComponenttcltk",
         )

## make tcltk S3 object S4 objects

oldclasses = c("tkwin", "tclVar")
for(i in oldclasses) {
  setOldClass(i)
  setIs(i,"guiWidgetORgWidgettcltkORtcltkObject")
}





##################################################
### Common methods.    Specific to a class are put into the file for that class

## we have two definitions. For instance, "svalue" and ".svalue". The "svalue" method dispatches on the object to the .svalue method. This allows us to use svalue instead of .svalue when defining the methods/constructors inside this package.


setMethod("svalue",signature(obj="gWidgettcltk"),
          function(obj, index=NULL, drop=NULL, ...) {
            .svalue(obj, obj@toolkit, ..., index=index, drop=drop)
          })



## svalue
## need method for character and AsIs
setMethod("svalue",signature(obj="character"),
          function(obj, index=NULL, drop=NULL, ...)  {
            ifelse(length(obj) == 1,
                   return(getObjectFromString(obj)),
                   return(obj)
                   )
          })
## method for Any is just a pass through
setMethod("svalue",signature(obj="ANY"),
          function(obj, index=NULL, drop=NULL, ...)  {
            return(obj)
          })


setMethod(".svalue",signature(toolkit = "guiWidgetsToolkittcltk", obj="character"),
          function(obj, toolkit, index=NULL, drop=NULL,  ...)  {
            ifelse(length(obj) == 1,
                   return(getObjectFromString(obj)),
                   return(NA)
                   )
          })

## svalue<- -- objec specific
setReplaceMethod("svalue",signature(obj="gWidgettcltk"),
          function(obj, index=NULL, ...,value) {
            .svalue(obj, obj@toolkit, index=index, ...) <- value
            obj
          })

                   
                 
## [
setMethod("[",
          signature(x="gWidgettcltk"),
          function(x,i,j,...,drop=TRUE) {
            
            return(.leftBracket(x, x@toolkit,i,j,...,drop=TRUE))
          })

## [<-
setReplaceMethod("[",signature(x="gWidgettcltk"),
          function(x,i,j,...,value) {
            if(missing(i) && missing(j))
              .leftBracket(x, x@toolkit,...) <- value
            else if(missing(j))
              .leftBracket(x, x@toolkit,i,...) <- value
            else 
              .leftBracket(x, x@toolkit,i,j,...) <- value
            return(x)
          })

## size ## return size -- not implemented
setMethod("size",signature(obj="gWidgettcltk"),
          function(obj, ...) {
            .size(obj, obj@toolkit,...)
          })

setMethod(".size", 
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
          function(obj, toolkit, ...) {
            width <- tclvalue(tkwinfo("width",getWidget(obj)))
            height <- tclvalue(tkwinfo("height",getWidget(obj)))

            return(as.numeric(c(width=width, height=height)))
          })

## size<-
setReplaceMethod("size",signature(obj="gWidgettcltk"),
          function(obj, ..., value) {
            .size(obj, obj@toolkit,...) <- value
            return(obj)
          })

setReplaceMethod(".size", 
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
                 function(obj, toolkit, ..., value) {
                   width <- value[1]
                   if(length(value) > 1)
                     height <- value[2]
                   else
                     height <- 0
                   if(height > 0)
                     tkconfigure(getWidget(obj), width=width, height=height)
                   else
                     tkconfigure(getWidget(obj), width=width)

                   return(obj)
                 })

## size for components is funny. For many width is characters, height
## is lines of text
setReplaceMethod(".size", 
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gComponenttcltk"),
                 function(obj, toolkit, ..., value) {
                   ## width in characters, height in lines
                   ## convert Pixels to each

                   ## XXX TODO: Hack in iterative process to fix size. This doesn't match
                   ## size(obj) <- width; size(obj)[1] - width == 0 (or even close)

                   ## simple way is
                   width <- value[1]
                   height <- NULL
                   if(length(value) > 1)
                     height <- value[2]
                   
                   ## ## set width
                   ## f <- function(lamda) {
                   ##   tkconfigure(getWidget(obj), width=ceiling(width*lamda/widthOfChar))
                   ##   act_width <- as.numeric(tkwinfo("width", getWidget(obj)))
                   ##   abs(act_width - width)
                   ## }
                   ## nlm(f, 1)#, stepmax=.05, steptol=.01)

                   ## if(!is.null(height)) {
                   ##   f <- function(char_height) {
                   ##     tkconfigure(getWidget(obj), height=ceiling(width/char_height))
                   ##     act_width <- as.numeric(tkwinfo("height", getWidget(obj)))
                   ##     abs(act_width - width)
                   ##   }
                   ##   nlm(f, heightOfChar, steptol=1)

                   ## }



                   width <- ceiling(width/widthOfChar)
                   if(!is.null(height))
                     height <- ceiling(height/heightOfChar)

                   if(!is.null(height))
                     tkconfigure(getWidget(obj), width=width, height=height)
                   else
                     tkconfigure(getWidget(obj), width=width)

                   return(obj)
                 })

## this works if container has no children (gwindow) but fails otherwise.
setReplaceMethod(".size", 
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gContainertcltk"),
                 function(obj, toolkit, ..., value) {
                   ## pixels for tkframe etc
                   width <- value[1]
                   if(length(value) > 1)
                     height <- value[2]
                   else
                     height <- 0
                   if(height > 0)
                     tkconfigure(getWidget(obj), width=width, height=height)
                   else
                     tkconfigure(getWidget(obj), width=width)

                   return(obj)
                 })



## visible
setMethod("visible",signature(obj="gWidgettcltk"),
          function(obj, set=NULL, ...) {
            .visible(obj,obj@toolkit, set=set, ...)
          })

##' get visibility
setMethod(".visible",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
          function(obj, toolkit, set=TRUE, ...) {
            widget <- getBlock(obj)
            if(is.null(set)) {
              ## return logical
              as.logical(tkwinfo("viewable", widget))
            } else {
              .visible(obj, toolkit) <- set
              return(NA)
            }
          })

## is widget viewable
setMethod(".visible",
                 signature(toolkit="guiWidgetsToolkittcltk",obj="tkwin"),
                 function(obj, toolkit, set=TRUE, ...) {
                   w <- getWidget(obj)
                   as.logical(tkwinfo("viewable", w))
                 })


## visible<-
setReplaceMethod("visible",signature(obj="gWidgettcltk"),
          function(obj, ..., value) {
            .visible(obj, obj@toolkit, ...) <- value
            return(obj)
          })

setReplaceMethod(".visible",
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
                 function(obj, toolkit, ..., value) {
                   message("visible<- not implemented\n")
                   return(obj)
                 })
setReplaceMethod(".visible",
                 signature(toolkit="guiWidgetsToolkittcltk",obj="tkwin"),
                 function(obj, toolkit, ..., value) {
                   ## visible not implemented
                   
                   return(obj)
                 })


## enabled -- TRUE If state is normal
##' enabled different for ttk widget
enabled_ttkwidget <- function(x, ...) {
  !as.logical(tcl(x, "instate", "disabled"))
}
enabled_tkwidget <- function(x, ...) {
  as.character(tkcget(x, "-state")) == "normal"
}

setMethod("enabled",signature(obj="gWidgettcltk"),
          function(obj, ...) {
            .enabled(obj, obj@toolkit,...)
          })
setMethod(".enabled",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
          function(obj, toolkit, ...) {
            w <- getWidget(obj)
            if(isTtkWidget(w))
              enabled_ttkwidget(w)
            else
              enabled_tkwidget(w)
          })

## enabled<-
setReplaceMethod("enabled",signature(obj="gWidgettcltk"),
          function(obj, ..., value) {
            .enabled(obj, obj@toolkit,...) <- value
            return(obj)
          })

setReplaceMethod(".enabled",
                 signature(toolkit="guiWidgetsToolkit",obj="gWidgettcltk"),
                 function(obj, toolkit, ..., value) {
                   .enabled(obj,guiToolkit("tcltk"), ...) <- value
                   return(obj)
                 })

setenabled_ttkwidget <- function(x, value) {
  tcl(x, "state", ifelse(value, "!disabled", "disabled"))
}
setenabled_tkwidget <- function(x, value) {
  tkconfigure(x, state=ifelse(as.logical(value), "normal", "disabled"))
}
setReplaceMethod(".enabled",
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
                 function(obj, toolkit, ..., value) {
                   w <- getWidget(obj)

                   if(isTtkWidget(w))
                     setenabled_ttkwidget(w, as.logical(value))
                   else
                     setenabled_tkwidget(w, as.logical(value))

                   ## recurse into childComponents
                   childComponents <- obj@e$childComponents
                   
                   if(!is.null(childComponents))
                     lapply(childComponents,function(i) enabled(i) <- value)
                            
                   return(obj)
                 })

## editable or readonly
## I want a editable<- method for gdf, gcombobox, glabel
setMethod(".editable",
          signature(toolkit="guiWidgetsToolkittcltk", obj="gWidgettcltk"),
          function(obj, toolkit) {
            widget <- getWidget(obj)            
            as.character(tkcget(widget, "-state")) != "readonly"
          })

setReplaceMethod(".editable",
                 signature(toolkit="guiWidgetsToolkittcltk", obj="gWidgettcltk", value="logical"),
                 function(obj, toolkit, ..., value) {
                   widget <- getWidget(obj)
                   tkconfigure(widget, "state"=ifelse(value, "normal", "readonly"))
                   return(obj)
                 })



## focus
focus_ttkwidget <- function(x, ...) as.logical(tcl(x, "instate", "focus"))
focus_tkwidget <- function(x, ...) {
  if(is.tkwin(x))
    x <- x$ID
  tl <- tkwinfo("toplevel", x)
  cur <- as.character(tcl("focus", displayof=tl))
  return(cur == x)
}

setMethod("focus",signature(obj="gWidgettcltk"),
          function(obj, ...) {
            .focus(obj, obj@toolkit,...)
          })

setMethod(".focus",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
          function(obj, toolkit, ...) {
            w <- getWidget(obj)
            if(isTtkWidget(w))
             focus_ttkwidget(w)
            else
              focus_tkwidget(w)
          })

## focus<-
setReplaceMethod("focus",signature(obj="gWidgettcltk"),
          function(obj, ..., value) {
            .focus(obj, obj@toolkit,...) <- value
            return(obj)
          })

setfocus_ttkwidget <- function(x, value) if(value) tcl(x, "state", "focus")
setfocus_tkwidget <- function(x, value) if(value) tkfocus(x)

setReplaceMethod("focus",signature(obj="tcltkObject"),
          function(obj, ..., value) {
            w <- getWidget(obj)
            if(isTtkWidget(w))
              setfocus_ttkwidget(w, value)
            else
              setfocus_tkwidget(w, value)
            return(obj)
          })


setReplaceMethod(".focus",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
          function(obj, toolkit, ..., value) {
            focus(obj@widget, toolkit, ...) <- value
            return(obj)
          })
                 

setReplaceMethod(".focus",
          signature(toolkit="guiWidgetsToolkittcltk",obj="tcltkObject"),
          function(obj, toolkit, ..., value) {
            value = as.logical(value)
            if(as.logical(value))
              tkfocus(getBlock(obj))

            return(obj)
          })

## default Widget is initially focused. SHould have binding for this
## defaultWidget
setMethod("defaultWidget",signature(obj="gWidgettcltk"),
          function(obj, ...) {
            .defaultWidget(obj, obj@toolkit,...)
          })

setMethod(".defaultWidget",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
          function(obj, toolkit, ...)
          focus(obj)
          )

## defaultWidget<-
setReplaceMethod("defaultWidget",signature(obj="gWidgettcltk"),
                 function(obj, ..., value) {
                   .defaultWidget(obj, obj@toolkit,...) <- value
                   return(obj)
                 })

setReplaceMethod("defaultWidget",signature(obj="tcltkObject"),
          function(obj, ..., value) {
            .defaultWidget(obj, toolkit=guiToolkit("tcltk"),...) <- value
            return(obj)
          })


setReplaceMethod(".defaultWidget",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
          function(obj, toolkit, ..., value) {
            widget <- getWidget(obj)
            .defaultWidget(widget, toolkit, ...) <- value
            return(obj)
          })

setReplaceMethod(".defaultWidget",
          signature(toolkit="guiWidgetsToolkittcltk",obj="tcltkObject"),
          function(obj, toolkit, ..., value) {
            value = as.logical(value)
            if(value)
              tkfocus(obj)
            return(obj)
          })


## isExtant



## enabled -- TRUE If state is normal
setMethod("isExtant",signature(obj="gWidgettcltk"),
          function(obj, ...) {
            .isExtant(obj, obj@toolkit,...)
          })
setMethod(".isExtant",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
          function(obj, toolkit, ...) {
            w <- getWidget(obj)
            as.logical(as.numeric(tkwinfo("exists", w)))
          })


## tooltip<-
setReplaceMethod("tooltip",signature(obj="gWidgettcltk"),
          function(obj, ..., value) {
            .tooltip(obj, obj@toolkit,...) <- value
            return(obj)
          })

setReplaceMethod(".tooltip",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
          function(obj, toolkit, ..., value) {
            widget <- getWidget(obj)
            tk2tip(widget, paste(value, collapse="\n"))
            return(obj)
          })

setReplaceMethod("tooltip",signature(obj="tcltkObject"),
          function(obj, ..., value) {
            ## set the tip.
            tk2tip(obj, paste(value, collapse="\n"))
            return(obj)
          })




## font
## The .font method is not imported from gWidgets, or exported from gWidgetstcltk. Add this bac if you are going to use this method

setMethod("font",signature(obj="gWidgettcltk"),
          function(obj, ...) {
            warning("font() not defined. Set fonts with font<-")
            return()
            .font(obj, obj@toolkit,...)
          })

## font<-
setReplaceMethod("font",signature(obj="gWidgettcltk"),
          function(obj, ..., value) {
            .font(obj, obj@toolkit,...) <- .fixFontMessUp(value)
            return(obj)
          })
setReplaceMethod(".font",
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
                 function(obj, toolkit, ..., value) {
                   .font(obj@widget, toolkit, ...) <- value
                   return(obj)
                 })

.font.styles = list(
  families = c("normal","sans","serif","monospace"),
  weights = c("normal","oblique","italic"),
  styles = c("ultra-light","light","normal","bold","ultra-bold","heavy"),
  colors = c("black","blue","red","green","brown","yellow","pink")
)

## common
merge.list <- function(x,y, overwrite=TRUE) {
  for(i in names(y)) {
    if(is.null(x[[i]]) || overwrite)
      x[[i]] <- y[[i]]
  }
  x
}

## ... passed into tkfont.create as font name
fontlistFromMarkup <- function(markup,...) {

  if(!is.list(markup))
    markup <- lapply(markup,function(x) x)
  
  fontList <- list(...)
  if(!is.null(markup$family))
    fontList <- merge(fontList, list(family=switch(markup$family,
                                       "normal"="times",
                                       "sans" = "helvetica",
                                       "serif" = "times",
                                       "monospace"="courier",
                                       markup$family)))
  if(!is.null(markup$style))
    fontList <- merge(fontList, list(slant=switch(markup$style,
                                       "normal"="roman",
                                       "oblique"="roman",
                                       "italic"="italic",
                                       "roman")))
  if(!is.null(markup$weight))
    fontList <- merge(fontList, list(weight=switch(markup$weight,
                                       "heavy"="bold",
                                       "ultra-bold"="bold",
                                       "bold"="bold",
                                       "normal"="normal",
                                       "light"="normal",
                                       "ultra-light" = "normal",
                                       "normal")))
  
  if(!is.null(markup$size))
    if(is.numeric(markup$size))
      fontList <- merge(fontList, list(size=markup$size))
    else
      fontList <- merge(fontList,list(size = switch(markup$size,
                                        "xxx-large"=24,
                                        "xx-large"=20,
                                        "x-large"=18,
                                        "large"=16,
                                        "medium"=12,
                                        "small"=10,
                                        "x-small"=8,
                                        "xx-small"=6,
                                        as.integer(markup$size))))
  return(fontList)
}




setReplaceMethod(".font",
                 signature(toolkit="guiWidgetsToolkittcltk",obj="tcltkObject"),
                 function(obj, toolkit, ..., value) {
                   fname <- paste(as.character(date()),rnorm(1), sep="") ## some random string
                   theArgs <- fontlistFromMarkup(value, fname)

                   ## now call
                   ## font with ttk is different -- fix XXX
                   theFont = do.call("tkfont.create",theArgs)
                   ret <- try(tkconfigure(getWidget(obj), font=fname), silent=TRUE)
                   ## colors are different
                   if("color" %in% names(value))
                       try(tkconfigure(getWidget(obj), foreground=value['color']), silent=TRUE) 
                   ## all done
                   return(obj)
                   
                 })



## tag, tag<-
## In RGtk2 we used the getData() and setData() methods. In tcltk we use the
## crummy implementation from rJava -- a list which grows without bound




## ## create namespace object
## tags = list()
## assignInNamespace("tags",list(),"gWidgetstcltk")

## ## clear out tags for this ID. Not exported. Is this used?
## Tagsclear = function(obj) {

##   id = obj@ID
  
##   tags = getFromNamespace("tags",ns="gWidgetstcltk")
##   allKeys = names(tags)

##   inds = grep(paste("^",id,"::",sep=""),allKeys)
##   if(length(inds) == 0)
##     return(NA)

##   ## else
##   tags[[inds]] <- NULL
##   assignInNamespace("tags",tags,ns="gWidgetstcltk")
## }


setMethod("tag",signature(obj="gWidgettcltk"),
          function(obj,i,drop=TRUE, ...) {
            if(missing(drop)) drop <- TRUE
            .tag(obj, obj@toolkit,i, drop=drop,...)
          })
## dispatch in *this* toolkit, not present in obj
setMethod("tag",signature(obj="tcltkObject"),
          function(obj,i,drop=TRUE, ...) {
            if(missing(drop)) drop <- TRUE            
            .tag(obj, guiToolkit("tcltk"),i, drop=drop,...)
          })

setMethod(".tag", signature(toolkit="guiWidgetsToolkittcltk",obj="guiWidget"),
          function(obj, toolkit, i, drop=TRUE, ...) {
            if(missing(i)) i = NULL
            if(missing(drop)) drop <- TRUE                        
            .tag(obj@widget,toolkit,  i, drop=drop,  ...)
          })
setMethod(".tag", signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
          function(obj, toolkit, i, drop=TRUE, ...) {
            if(missing(i)) i = NULL
            if(missing(drop)) drop <- TRUE                                    

            if(is.null(i))
              return(as.list(obj@e))
            else
              return(obj@e[[i]])
            
          })

## tag <-
setReplaceMethod("tag",signature(obj="gWidgettcltk"),
          function(obj, i, replace=TRUE, ..., value) {
            .tag(obj, obj@toolkit,i,replace, ...) <- value
            return(obj)
          })
## dispatch in *this* toolkit, not present in obj
setReplaceMethod("tag",signature(obj="tcltkObject"),
          function(obj,i, replace=TRUE, ..., value) {
            .tag(obj, guiToolkit("tcltk"),i, replace, ...) <- value
            return(obj)
          })

## objects can be in many different flavors: guiWIdget, gWidgettcltk, tcltkObject
setReplaceMethod(".tag", signature(toolkit="guiWidgetsToolkittcltk",obj="guiWidget"),
          function(obj, toolkit, i, replace=TRUE, ..., value) {
            if(missing(i)) i = NULL
            .tag(obj@widget,toolkit,  i, replace, ...) <- value
            return(obj)
          })

setReplaceMethod(".tag", signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
          function(obj, toolkit, i, replace=TRUE, ..., value) {
            if(missing(i)) i = NULL
            

            obj@e[[i]] <- value
            return(obj)

          })

##################################################
## id -- define for "ANY" as well
setMethod("id",signature(obj="gWidgettcltk"),
          function(obj, ...) {
            tag(obj,".tcltkID")
          })
setMethod("id",signature(obj="tcltkObject"),
          function(obj, ...) {
            tag(obj, ".tcltkID", ...)
            return(obj)
          })
setMethod("id",signature(obj="ANY"),
          function(obj, ...) {
            if(!is.null(theID<- attr(obj,"id"))) {
              return(theID)
            } else {
              if(is.character(obj)) {
                return(obj[1])
              } else {
                dps = deparse(substitute(obj))
                attr(obj,"id") <- dps
                return(dps)
              }
            }
          })


setMethod(".id", signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
          function(obj, toolkit, ...) {
            tag(obj,".tcltkID", ...)
          })
setMethod(".id", signature(toolkit="guiWidgetsToolkittcltk",obj="tcltkObject"),
          function(obj, toolkit,  ...) {
            return(tag(obj,".tcltkID"))
          })


## id<-
setReplaceMethod("id",signature(obj="gWidgettcltk"),
          function(obj, ..., value) {
            tag(obj,".tcltkID", ...) <- value
            return(obj)
          })
## dispatch in *this* toolkit, not present in obj
setReplaceMethod("id",signature(obj="tcltkObject"),
          function(obj, ..., value) {
            tag(obj, ".tcltkID", ...) <- value
            return(obj)
          })
setReplaceMethod("id",signature(obj="ANY"),
          function(obj, ..., value) {
            attr(obj,"id") <- value
            return(obj)
          })


## we need a .id to handle dispatch from guiWidgets, otherwise, we use id()
setReplaceMethod(".id", signature(toolkit="guiWidgetsToolkittcltk",
                                  obj="gWidgettcltk"),
          function(obj, toolkit, ..., value) {
            id(obj, ...) <- value
            return(obj)
          })



## add method is biggie
## we have several levels of classes here guiWidget -- gWidgetRGkt -- tcltkObject, when
## we get down to that level we can finally add
setMethod("add",signature(obj="gWidgettcltk"),
          function(obj, value, ...) {
            .add(obj, obj@toolkit,value,...)
          })
setMethod(".add",
          signature(toolkit="guiWidgetsToolkittcltk",
                    obj="guiWidget", value="ANY"),
          function(obj, toolkit, value, ...) {
            gwCat(gettext("Can't add without a value\n"))
          })
setMethod(".add",
          signature(toolkit="guiWidgetsToolkittcltk",
                    obj="gWidgettcltk", value="try-error"),
          function(obj, toolkit, value, ...) {
            gmessage(paste("Error:",obj))
          })
## pushdonw
setMethod(".add",
          signature(toolkit="guiWidgetsToolkittcltk",
                    obj="guiWidget", value="guiWidgetORgWidgettcltkORtcltkObject"),
          function(obj, toolkit, value, ...) {
            .add(obj@widget, toolkit, value, ...)
          })

## for gWindow
setMethod(".add",
          signature(toolkit="guiWidgetsToolkittcltk",
                    obj="gContainertcltk", value="guiWidget"),
          function(obj, toolkit, value, ...) {
            .add(obj, toolkit, value@widget, ...)
          })

## for gContainer
setMethod(".add",
          signature(toolkit="guiWidgetsToolkittcltk", obj="gContainertcltk",value="gWidgettcltk"),
          function(obj, toolkit, value, ...) {

            ## add parent, children
            childComponents <- obj@e$childComponents
            if(is.null(childComponents))
              childComponents <- list()
            obj@e$childComponents <- c(childComponents, value)
            value@e$parentContainer <- obj

            ## inherit enabled from parent
            try(.enabled(value,toolkit) <- .enabled(obj,toolkit),silent=TRUE)
            
            theArgs = list(...)
            ## passed to do.call. Populate this
            argList = list(getBlock(value))

            ## expand, fill, anchor
            ## XXX make expand option default to TRUE
            expand <- getWithDefault(theArgs$expand,
                                     getWithDefault(getOption("gw:tcltkDefaultExpand", FALSE)))

            ## fill
            horizontal <- obj@horizontal
            fill <- getWithDefault(theArgs$fill, ifelse(horizontal, "y", "x")) # FALSE, x, y, both=TRUE
            if(is.logical(fill)) {
              if(fill)
                fill <- "both"
              else
                fill <- NULL
            }

            ## the default anchor. -1,0 or W makes layouts nicer looking IMHO
            defaultAnchor <- getWithDefault(getOption("gw:tcltkDefaultAnchor"), c(-1, 0))
            anchor <- xyToAnchor(getWithDefault(theArgs$anchor, defaultAnchor))

            ## expand: if TRUE then can either anchor or fill. If 
            if(!expand) {
              fill <- NULL
            }


            argList$expand <- expand
            argList$fill <- fill
            argList$anchor <- anchor

            
            if(obj@horizontal)
              argList$side = "left"
            else
              argList$side = "top"

            ## call tkpack
            do.call("tkpack",argList)

            tcl("update","idletasks")
            
            if(!is.null(widget <- .tag(obj,toolkit, "scrollable.widget"))) {
              ## get scrollbars to add to end etc.
              tcl("event","generate",getWidget(obj),"<Configure>")
              tkxview.moveto(widget,1)
              tkyview.moveto(widget,1)
            }
          })

##' Add method to incorporate tk widget into gui:
##'
##' @example
##' g = ggroup(cont=gwindow())
##' library(tkrplot)
##' l = tkrplot(getToolkitWidget(g), function() hist(rnorm(100)))
##' add(g, l)
setMethod(".add",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gContainertcltk", value="tkwin"),
          function(obj, toolkit, value, ...) {
            tkpack(value, expand=TRUE, fill="both")
          })




## addSPring, addSpace
setMethod("addSpring",signature(obj="gWidgettcltk"),
          function(obj, ...) {
            .addSpring(obj, obj@toolkit,...)
          })

setMethod(".addSpring",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gContainertcltk"),
          function(obj, toolkit, ...) {

            tt <- getBlock(obj)
            blankLabel <- ttklabel(tt,text=" ")

            if(obj@horizontal) {
              ## doesn't work!
              tkpack(blankLabel, expand=TRUE, fill="x", side="left")
            } else {
              tkpack(blankLabel, expand=TRUE, fill="y", side="top")
            }
            invisible()
          })

setMethod("addSpace",signature(obj="gWidgettcltk"),
          function(obj, value, ...) {
            .addSpace(obj,obj@toolkit,value,...)
          })

setMethod(".addSpace",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gContainertcltk"),
          function(obj, toolkit, value, ...) {
            theArgs = list(...)
            horizontal = ifelse(is.null(theArgs$horizontal),
              TRUE,
              as.logical(theArgs$horizontal))

            tt <- getBlock(obj)
            value = as.integer(value)
            if(horizontal)
              tkpack(ttklabel(tt, text=""),side="left",padx=value)
            else
              tkpack(ttklabel(tt, text=""),side="top", pady=value)
            invisible()
          })

## delete -- get down to two tcltkObjects
setMethod("delete",signature(obj="gWidgettcltk"),
          function(obj, widget, ...) {
            .delete(obj, obj@toolkit,widget,...)
          })

## push down to tcltk vs tcltk. Can be 9 possiblities!
setMethod(".delete",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gContainertcltk",widget="guiWidget"),
          function(obj, toolkit, widget, ...) {
            .delete(obj, toolkit, widget@widget, ...)
          })
setMethod(".delete",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gContainertcltk",widget="gWidgettcltk"),
          function(obj, toolkit, widget, ...) {
            ## call remove on container
            tkpack.forget(getBlock(widget))
          })

## dispose -- delete the parent window, or something else
setMethod("dispose",signature(obj="gWidgettcltk"),
          function(obj, ...) {
            .dispose(obj, obj@toolkit,...)
          })

setMethod(".dispose",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
          function(obj, toolkit, ...) {
            tcl("after",5,function() {
              tt <- getTopParent(getBlock(obj))
              tkgrab.release(tt)
              tkdestroy(tt)
            })                          # wait a pause 
          })




## update
setMethod("update",signature(object="gWidgettcltk"),
          function(object, ...) {
            .update(object, object@toolkit, ...)
          })

setMethod(".update",
          signature(toolkit="guiWidgetsToolkittcltk",object="gComponenttcltk"),
          function(object, toolkit, ...) {

            missingMsg(".update");return()
            
            object@widget$QueueDraw()
          })

##
##
##################################################


##################################################
## handlers. Also in aaaHandlers
##
## basic handler for adding with a signal. Not exported.
setGeneric("addhandler", function(obj, signal, handler, action=NULL, ...)
           standardGeneric("addhandler"))
setMethod("addhandler",signature(obj="guiWidget"),
          function(obj, signal, handler, action=NULL, ...) {
            .addhandler(obj@widget, obj@toolkit, signal, handler, action, ...)
          })
setMethod("addhandler",signature(obj="gWidgettcltk"),
          function(obj, signal, handler, action=NULL, ...) {
            .addhandler(obj, obj@toolkit, signal, handler, action, ...)
          })
setMethod("addhandler",signature(obj="tcltkObject"),
          function(obj, signal, handler, action=NULL, ...) {
            .addhandler(obj, guiToolkit("tcltk"), signal, handler, action, ...)
          })

## method for dispatch
setGeneric(".addhandler",
           function(obj, toolkit,
                  signal, handler, action=NULL, ...)
           standardGeneric(".addhandler"))


setMethod(".addhandler",
          signature(toolkit="guiWidgetsToolkittcltk",obj="guiWidget"),
          function(obj, toolkit,
                   signal, handler, action=NULL, ...) {
            .addhandler(obj@widget, force(toolkit), signal, handler, action, ...)
          })

setMethod(".addhandler",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
          function(obj, toolkit,
                   signal, handler, action=NULL, ...) {
            .addHandler(obj, force(toolkit), signal, handler, action, ...)
          })




## Make upcase for Handler
setGeneric(".addHandler",
           function(obj, toolkit,
                  signal, handler, action=NULL, ...)
           standardGeneric(".addHandler"))


setMethod(".addHandler",
          signature(toolkit="guiWidgetsToolkittcltk",obj="guiWidget"),
          function(obj, toolkit,
                   signal, handler, action=NULL, ...) {
            .addhandler(obj@widget, force(toolkit), signal, handler, action, ...)
          })




## removew handler
## removehandler
setMethod("removehandler", signature("gWidgettcltk"),
          function(obj, ID=NULL, ...) {
            .removehandler(obj, obj@toolkit, ID, ...)
          })
setMethod("removehandler", signature("tcltkObject"),
          function(obj, ID=NULL, ...) {
            .removehandler(obj, guiToolkit("tcltk"), ID, ...)
          })

## in aaaHandlers
## setMethod(".removehandler",
##           signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
##           function(obj, toolkit, ID=NULL, ...) {

##             ## ID here has two components
##             type = ID[1]
##             handlerID=as.character(ID[2])
##             ID = as.character(obj@ID)

##             ## remove from list
##             allHandlers = getFromNamespace("allHandlers",ns="gWidgetstcltk")

##             ## is this a idleHandler
##             if(type == "addIdleListener") {
##               t = allHandlers[[ID]][[type]][[handlerID]]$timer
##               t = .jcall(t,"V","stopTimer")
##             }
##             allHandlers[[ID]][[type]][[handlerID]]<-NULL
##               ## now store the hash
##             assignInNamespace("allHandlers",allHandlers,ns="gWidgetstcltk")
##           })


## blockhandler
setMethod("blockhandler", signature("gWidgettcltk"),
          function(obj, ID=NULL, ...) {
            .blockhandler(obj, obj@toolkit, ID, ...)
          })
setMethod("blockhandler", signature("tcltkObject"),
          function(obj, ID=NULL, ...) {
            .blockhandler(obj, guiToolkit("tcltk"), ID, ...)
          })

## setMethod(".blockhandler",
##           signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
##           function(obj, toolkit, ID=NULL, ...) {
##             .blockhandler(getWidget(obj),toolkit,ID,...)
##           })

## setMethod(".blockhandler",
##           signature(toolkit="guiWidgetsToolkittcltk",obj="tcltkObject"),
##           function(obj, toolkit, ID=NULL, ...) {
##             gwCat(gettext("define block handler\n"))
##           })

## unblock handler
setMethod("unblockhandler", signature("gWidgettcltk"),
          function(obj, ID=NULL, ...) {
            .unblockhandler(obj, obj@toolkit, ID, ...)
          })
setMethod("unblockhandler", signature("tcltkObject"),
          function(obj, ID=NULL, ...) {
            .unblockhandler(obj, guiToolkit("tcltk"), ID, ...)
          })

## setMethod(".unblockhandler",
##           signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
##           function(obj, toolkit, ID=NULL, ...) {
##             .blockhandler(getWidget(obj),toolkit,ID,...)
##           })

## setMethod(".unblockhandler",
##           signature(toolkit="guiWidgetsToolkittcltk",obj="tcltkObject"),
##           function(obj, toolkit, ID=NULL, ...) {
##             cat("define unblock handler\n")
##           })



## addhandlerchanged
setMethod("addhandlerchanged",signature(obj="gWidgettcltk"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerchanged(obj, obj@toolkit, handler, action, ...)
          })
setMethod("addhandlerchanged",signature(obj="tcltkObject"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerchanged(obj, guiToolkit("tcltk"), handler, action, ...)
          })
setMethod("addhandlerchanged",signature(obj="ANY"),
          function(obj, handler=NULL, action=NULL, ...) {
            warning("No method addhandlerchanged for object of class",class(obj),"\n")
          })

setMethod(".addhandlerchanged",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
          function(obj, toolkit,
                   handler, action=NULL, ...) {
            .addHandler(obj, toolkit, signal="<KeyPress>",
                        handler=handler, action=action, ...)
          })


## expose: expose-event or realize
setMethod("addhandlerexpose",signature(obj="gWidgettcltk"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerexpose(obj,obj@toolkit,handler, action, ...)
          })
setMethod("addhandlerexpose",signature(obj="tcltkObject"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerexpose(obj, guiToolkit("tcltk"), handler, action, ...)
          })

setMethod(".addhandlerexpose",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
          function(obj, toolkit,
                   handler, action=NULL, ...) {
            .addHandler(obj, toolkit, signal="<Expose>",
                        handler=handler, action=action, ...)
          })

setMethod(".addhandlerexpose",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gComponenttcltk"),
          function(obj, toolkit,
                   handler, action=NULL, ...) {
            .addHandler(obj,toolkit, signal="<Realize>",
                        handler=handler, action=action, ...)
          })

## unrealize: unrealize
setMethod("addhandlerunrealize",signature(obj="gWidgettcltk"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerunrealize(obj, obj@toolkit,handler, action, ...)
          })
setMethod("addhandlerunrealize",signature(obj="tcltkObject"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerunrealize(obj, guiToolkit("tcltk"),handler, action, ...)
          })

setMethod(".addhandlerunrealize",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
          function(obj, toolkit,
                   handler, action=NULL, ...) {
            .addHandler(obj, toolkit, signal="<Unmap>",
                        handler=handler, action=action, ...)
          })

## destroy: destroy
setMethod("addhandlerdestroy",signature(obj="gWidgettcltk"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerdestroy(obj, obj@toolkit,handler, action, ...)
          })
setMethod("addhandlerdestroy",signature(obj="tcltkObject"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerdestroy(obj, guiToolkit("tcltk"),handler, action, ...)
          })

setMethod(".addhandlerdestroy",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
          function(obj, toolkit,
                   handler, action=NULL, ...) {
            .addHandler(obj, toolkit, signal="<Destroy>",
                        handler=handler, action=action, ...)
          })

## keystroke: changed
setMethod("addhandlerkeystroke",signature(obj="gWidgettcltk"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerkeystroke(obj, obj@toolkit,handler, action, ...)
          })
setMethod("addhandlerkeystroke",signature(obj="tcltkObject"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerkeystroke(obj, guiToolkit("tcltk"),handler, action, ...)
          })

## setMethod(".addhandlerkeystroke",
##           signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
##           function(obj, toolkit,
##                    handler, action=NULL, ...) {
##             .addHandler(obj, toolkit, signal="<KeyPress>",
##                         handler=handler, action=action, ...)
##           })

##' for gedit, gtext
##'
##' This uses the FUN argument for .addHandler to pass in a special function This way the arguments
##' that tcltk passes in can be used
##' 
##' %K The keysym corresponding to the event, substituted as a textual string. Valid only for
##'     KeyPress and KeyRelease events.
##' This handler can not be blocked or removed!
##' 
setMethod(".addhandlerkeystroke",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
          function(obj,toolkit, handler=NULL, action=NULL,...) {
            .addHandler(obj,toolkit,"<KeyRelease>",handler,action,
                        FUN = function(K) {
                          h = list(obj = obj, action = action, key=K)
                          runHandlers(obj, "<KeyRelease>", h, ...)
                          ## handler(h)
                        })
          })


## clicked: clicked
setMethod("addhandlerclicked",signature(obj="gWidgettcltk"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerclicked(obj, obj@toolkit,handler, action, ...)
          })
setMethod("addhandlerclicked",signature(obj="tcltkObject"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerclicked(obj, guiToolkit("tcltk"),handler, action, ...)
          })

setMethod(".addhandlerclicked",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
          function(obj, toolkit,
                   handler, action=NULL, ...) {
            .addHandler(obj, toolkit, signal="<Button-1>",
                        handler=handler, action=action, ...)
          })

## doubleclick: no default
setMethod("addhandlerdoubleclick",signature(obj="gWidgettcltk"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerdoubleclick(obj,obj@toolkit,handler, action, ...)
          })
setMethod("addhandlerdoubleclick",signature(obj="tcltkObject"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerdoubleclick(obj,guiToolkit("tcltk"),handler, action, ...)
          })

setMethod(".addhandlerdoubleclick",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
          function(obj, toolkit,
                   handler, action=NULL, ...) {
            .addHandler(obj, toolkit, signal="<Double-Button-1>",
                        handler=handler, action=action, ...)
          })

## rightclick: button-press-event -- handle separately
setMethod("addhandlerrightclick",signature(obj="gWidgettcltk"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerrightclick(obj,obj@toolkit,handler, action, ...)
          })
setMethod("addhandlerrightclick",signature(obj="tcltkObject"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerrightclick(obj,guiToolkit("tcltk"),handler, action, ...)
          })

setMethod(".addhandlerrightclick",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
          function(obj, toolkit,
                   handler, action=NULL, ...) {

            if(windowingsystem() == "aqua" ||
               grepl("^mac",.Platform$pkgType)) {
              id <- lapply(c("<Control-1>", "<Button-2>", "<Button-3>"), function(i) {
                id <- .addHandler(obj, toolkit, signal=i,
                            handler=handler, action=action, ...)
                list(obj=obj, id=id, signal=i)    # for remove/block/unblock
              })
            } else {
              id <- .addHandler(obj, toolkit, signal="<Button-3>", 
                          handler=handler, action=action, ...)
            }
            invisible(id)
          })


## focus
setMethod("addhandlerfocus",signature(obj="gWidgettcltk"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerfocus(obj,obj@toolkit,handler, action, ...)
          })
setMethod("addhandlerfocus",signature(obj="tcltkObject"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerfocus(obj,guiToolkit("tcltk"),handler, action, ...)
          })

setMethod(".addhandlerfocus",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
          function(obj, toolkit,
                   handler, action=NULL, ...) {
            .addHandler(obj, toolkit, signal="<FocusIn>",
                        handler=handler, action=action, ...)
          })

##' blur: blur should be focus out but is mouse out and focus out, so called twice!
setMethod("addhandlerblur",signature(obj="gWidgettcltk"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerblur(obj,obj@toolkit,handler, action, ...)
          })
setMethod("addhandlerblur",signature(obj="tcltkObject"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerblur(obj,guiToolkit("tcltk"),handler, action, ...)
          })

##' blur is on mouse motion here and change in focus. Handle is called twice!!!
setMethod(".addhandlerblur",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
          function(obj, toolkit,
                   handler, action=NULL, ...) {
            IDs <- lapply(c("<FocusOut>", "<Leave>"), function(i)
                          .addHandler(obj, toolkit, signal=i,
                                      handler=handler, action=action, ...)
                          )
            IDs
          })



## mousemotion
setMethod("addhandlermousemotion",signature(obj="gWidgettcltk"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlermousemotion(obj,obj@toolkit,handler, action, ...)
          })
setMethod("addhandlermousemotion",signature(obj="tcltkObject"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlermousemotion(obj,guiToolkit("tcltk"),handler, action, ...)
          })

setMethod(".addhandlermousemotion",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
          function(obj, toolkit,
                   handler, action=NULL, ...) {
            .addHandler(obj, toolkit, signal="<Motion>",
                        handler=handler, action=action, ...)
          })

## idle
setMethod("addhandleridle",signature(obj="gWidgettcltk"),
          function(obj, handler=NULL, action=NULL, interval=1000, ...) {
            .addhandleridle(obj, obj@toolkit,
                            handler=handler, action=action, interval=interval, ...)
          })
setMethod("addhandleridle",signature(obj="tcltkObject"),
          function(obj, handler=NULL, action=NULL, interval=1000, ...) {
            .addhandleridle(obj, guiToolkit("tcltk"),
                            handler=handler, action=action, interval=interval, ...)
          })


## addpopumenu
## ## this does not get exported
.addPopupMenu = function(obj,   menulist, action=NULL,...) {
  editPopupMenu <- getWidget(gmenu(menulist, popup=TRUE, action=action,container=obj, ...))
            
  RightClick <- function(x,y) # x and y are the mouse coordinates
    {
      V = getWidget(obj)
      rootx <- as.integer(tkwinfo("rootx",V))
      rooty <- as.integer(tkwinfo("rooty",V))
      xTxt <- as.integer(x)+rootx
      yTxt <- as.integer(y)+rooty
      tcl("tk_popup",editPopupMenu,xTxt,yTxt)
              }
  tkbind(getWidget(obj), "<Button-1>",RightClick)
}
  
.add3rdMousePopupMenu = function(obj,  menulist, action=NULL, ...) {

   
  editPopupMenu <- getWidget(gmenu(menulist, popup=TRUE, action=action,container=obj, ...))
            
  RightClick <- function(x,y) # x and y are the mouse coordinates
    {
      V = getWidget(obj)
      rootx <- as.integer(tkwinfo("rootx",V))
      rooty <- as.integer(tkwinfo("rooty",V))
      xTxt <- as.integer(x)+rootx
      yTxt <- as.integer(y)+rooty
      tcl("tk_popup",editPopupMenu,xTxt,yTxt)
    }
  W <- getWidget(obj)
  if(isMac()) {
    tkbind(W, "<Button-2>", RightClick)
    tkbind(W, "<Control-1>", RightClick)
  } else {
    tkbind(W, "<Button-3>",RightClick)
  }
}

setMethod("addpopupmenu",signature(obj="gWidgettcltk"),
          function(obj, menulist, action=NULL, ...) {
            .addpopupmenu(obj, obj@toolkit,menulist, action, ...)
          })
setMethod("addpopupmenu",signature(obj="tcltkObject"),
          function(obj, menulist, action=NULL, ...) {
            .addpopupmenu(obj, guiToolkit("tcltk"), menulist, action, ...)
          })


  
### 
setMethod(".addpopupmenu",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
          function(obj, toolkit, menulist, action=NULL, ...) {
            .addPopupMenu(obj, menulist, action=action, ...)
})


## add3rdmousepopupmenu
setMethod("add3rdmousepopupmenu",signature(obj="gWidgettcltk"),
          function(obj, menulist, action=NULL, ...) {
            .add3rdmousepopupmenu(obj, obj@toolkit,menulist, action, ...)
          })

setMethod("add3rdmousepopupmenu",signature(obj="tcltkObject"),
          function(obj, menulist, action=NULL,...) {
            .add3rdmousepopupmenu(obj, guiToolkit("tcltk"),menulist, action,...)
          })

setMethod(".add3rdmousepopupmenu",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gWidgettcltk"),
          function(obj, toolkit, menulist,action=NULL, ...) {
            .add3rdMousePopupMenu(obj,  menulist, action, ...)
          })
setMethod(".add3rdmousepopupmenu",
          signature(toolkit="guiWidgetsToolkittcltk",obj="tcltkObject"),
          function(obj, toolkit, menulist, action=NULL, ...) {
            .add3rdMousePopupMenu(obj, menulist, action, ...)
          })


## "dotmethods" defined in dnd.R
## adddropsource
setMethod("adddropsource",signature(obj="gWidgettcltk"),
          function(obj, targetType="text", handler=NULL, action=NULL, ...) {
            .adddropsource(obj, obj@toolkit,targetType=targetType,
                           handler=handler, action=action, ...)
          })
setMethod("adddropsource",signature(obj="tcltkObject"),
          function(obj, targetType="text", handler=NULL, action=NULL, ...) {
            .adddropsource(obj, guiToolkit("tcltk"),targetType=targetType,
                           handler=handler, action=action, ...)
          })

## adddropmotion
setMethod("adddropmotion",signature(obj="gWidgettcltk"),
          function(obj,  handler=NULL, action=NULL, ...) {
            .adddropmotion(obj, obj@toolkit,
                           handler=handler, action=action, ...)
          })
setMethod("adddropmotion",signature(obj="tcltkObject"),
          function(obj, handler=NULL, action=NULL, ...) {
            .adddropmotion(obj, guiToolkit("tcltk"),
                           handler=handler, action=action, ...)
          })

## adddroptarget
setMethod("adddroptarget",signature(obj="gWidgettcltk"),
          function(obj, targetType="text", handler=NULL, action=NULL, ...) {
            .adddroptarget(obj, obj@toolkit,targetType=targetType,
                           handler=handler, action=action, ...)
          })

setMethod("adddroptarget",signature(obj="tcltkObject"),
          function(obj, targetType="text", handler=NULL, action=NULL, ...) {
            .adddroptarget(obj, guiToolkit("tcltk"),targetType=targetType,
                           handler=handler, action=action, ...)
          })


## R Methods
setMethod("dim", "gWidgettcltk", function(x) .dim(x,x@toolkit))
setMethod(".dim",
          signature(toolkit="guiWidgetsToolkittcltk",x="gWidgettcltk"),
          function(x,toolkit) {
            gwCat(sprintf("Define dim for x of class: %s", class(x)[1]))
            return(NULL)
})
setMethod("length", "gWidgettcltk", function(x) .length(x,x@toolkit))
setMethod(".length",
#          signature(toolkit="guiWidgetsToolkittcltk"),
          signature(toolkit="ANY",x="ANY"),
          function(x,toolkit) {
#            gwCat(sprintf("Define length for x of class:%s\n"),class(x)[1])
            #return(NULL)
#            message("calling length for class", class(x)[1])
            return(NA)
})
          
setMethod("dimnames", "gWidgettcltk", function(x) .dimnames(x,x@toolkit))
setReplaceMethod("dimnames",
                 signature(x="gWidgettcltk"),
                 function(x,value) {
                   .dimnames(x,x@toolkit) <- value
                   return(x)
                 })
## as of 2.5.0 this became primiive
if(as.numeric(R.Version()$major) <= 2 &
   as.numeric(R.Version()$minor) <= 4.1) {
  setGeneric("names")
  setGeneric("names<-")
}

setMethod("names", "gWidgettcltk", function(x) .names(x,x@toolkit))
setReplaceMethod("names",
                 signature(x="gWidgettcltk"),
                 function(x,value) {
                   .names(x,x@toolkit) <- value
                   return(x)
                 })

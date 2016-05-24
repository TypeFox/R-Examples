## toolkit class
## register classes here for toolkits
## setClass("guiWidgetsToolkitRGtk2",
##          contains="guiWidgetsToolkit",
##          prototype=prototype(new("guiWidgetsToolkit"))
##          )



##################################################
## put S3 classes from RGtk2 into S4 classes
## got these from apropos("New") -> try(class(do.call(i,list())))
oldClasses =
c(
  "AtkNoOpObjectFactory",
  "AtkObjectFactory",
  "AtkRelationSet",
  "AtkStateSet",
  "GBoxed",
  "GObject",
  "GScanner",
  "GdkDragContext",
  "GdkPixbufLoader",
  "GdkRegion",
  "GtkAboutDialog",
  "GtkAccelGroup", 
  "GtkAccelLabel",
  "GtkAction",
  "GtkActionGroup",
  "GtkAdjustment",
  "GtkAlignment",
  "GtkArrow",
  "GtkAspectFrame",
  "GtkBin",
  "GtkBox",
  "GtkButton",
  "GtkButtonBox",
  "GtkCList",
  "GtkCTree",
  "GtkCalendar",
  "GtkCellRenderer",
  "GtkCellRendererCombo",
  "GtkCellRendererPixbuf",
  "GtkCellRendererProgress",
  "GtkCellRendererText",
  "GtkCellRendererToggle",
  "GtkCellView",
  "GtkCheckButton",
  "GtkCheckMenuItem",
  "GtkColorButton",
  "GtkColorSelection",
  "GtkColorSelectionDialog",
  "GtkCombo",
  "GtkComboBox",
  "GtkComboBoxEntry",
  "GtkContainer",
  "GtkCurve",
  "GtkDialog",
  "GtkDrawingArea",
  "GtkEntry",
  "GtkEntryCompletion",
  "GtkEventBox",
  "GtkExpander",
  "GtkFileFilter",
  "GtkFileSelection",
  "GtkFileChooserWidget",
  "GtkFixed",
  "GtkFontButton",
  "GtkFontSelection",
  "GtkFontSelectionDialog",
  "GtkFrame",
  "GtkGammaCurve",
  "GtkHBox",
  "GtkHButtonBox",
  "GtkHPaned",
  "GtkHRuler",
  "GtkHScale",
  "GtkHScrollbar",
  "GtkHSeparator",
  "GtkHandleBox",
  "GtkIMContext",
  "GtkIMContextSimple",
  "GtkIMMulticontext",
  "GtkIconFactory",
  "GtkIconSet",
  "GtkIconSource",
  "GtkIconTheme",
  "GtkIconView",
  "GtkImage",
  "GtkImageMenuItem",
  "GtkInfoBar",
  "GtkInputDialog",
  "GtkInvisible",
  "GtkItem",
  "GtkLabel",
  "GtkLayout",
  "GtkList",
  "GtkListItem",
  "GtkMenu",
  "GtkMenuBar",
  "GtkMenuItem",
  "GtkMenuShell",
  "GtkMisc",
  "GtkNotebook",
  "GtkObject",
  "GtkOptionMenu",
  "GtkPaned",
  "GtkProgress",
  "GtkProgressBar",
  "GtkRadioAction",
  "GtkRadioButton",
  "GtkRange",
  "GtkRcStyle",
  "GtkRuler",
  "GtkScale",
  "GtkScrollbar",
  "GtkScrolledWindow",
  "GtkSeparator",
  "GtkSeparatorMenuItem",
  "GtkSeparatorToolItem",
  "GtkSizeGroup",
  "GtkSocket",
  "GtkSpinButton",
  "GtkStatusbar",
  "GtkStyle",
  "GtkTable",
  "GtkTearoffMenuItem",
  "GtkTextAttributes",
  "GtkTextBuffer",
  "GtkTextChildAnchor",
  "GtkTextTag",
  "GtkTextTagTable",
  "GtkTextView",
  "GtkTipsQuery",
  "GtkToggleAction",
  "GtkToggleButton",
  "GtkToggleToolButton",
  "GtkToolButton",
  "GtkToolItem",
  "GtkToolbar",
  "GtkTooltips",
  "GtkTreeModelSort",
  "GtkTreePath",
  "GtkTreeStore",
  "GtkTreeView",
  "GtkTreeViewColumn",
  "GtkUIManager",
  "GtkVBox",
  "GtkVButtonBox",
  "GtkVPaned",
  "GtkVRuler",
  "GtkVScale",
  "GtkVScrollbar",
  "GtkVSeparator",
  "GtkViewport",
  "GtkWidget",
  "GtkWindow",
  "GtkWindowGroup",
  ##
  "PangoAttrList",
  "PangoCairoFcFontMap",
  "PangoCoverage",
  "PangoFcFontMap",
  "PangoFontDescription",
  "PangoFontMap",
  "PangoGlyphString",
  "PangoItem",
  "GObject",
  "RGtkDataFrame",
  ## add in others that come up
  "GtkTreeSelection"
  )


setOldClass("RGtkObject")
lapply(oldClasses, function(i) {
  setOldClass(i)
  setIs(i,"RGtkObject")
})

setOldClass("try-error")                # for handling try-errors


## a base class which is virtual


##################################################
## A virtual class to hold either RGTK or these guys

## A virtual class for our newly defined objects
setClass("gWidgetRGtk")


## this worked in 2.4.0 but not later, remove
##setIs("guiWidget","guiWidgetORgWidgetRGtkORRGtkObject")
##setIs("gWidgetRGtk","guiWidgetORgWidgetRGtkORRGtkObject")
##setIs("RGtkObject","guiWidgetORgWidgetRGtkORRGtkObject")

setClassUnion("guiWidgetORgWidgetRGtkORRGtkObject",
              c("guiWidget","gWidgetRGtk","RGtkObject"))

## subclss
## This behaviour changed in R as of 2.7.0. We throw in the towel here and use "ANY and not the preferred guiWidgetORgWidgetRGtkORRGtkObject

##.Rversion <- R.Version()
##if(.Rversion$major == "2" && .Rversion$minor == "7.0" ) {
setClass("gComponentRGtk",
         representation(
                        block="ANY",
                        widget="ANY",
                        toolkit="guiWidgetsToolkit"
                        ),
         contains="gWidgetRGtk",
         )
setClass("gContainerRGtk",
         representation(
                        block="ANY",
                        widget="ANY",
                        toolkit="guiWidgetsToolkit"
                   ),
         contains="gWidgetRGtk",
         )
## } else {
##   setClass("gComponentRGtk",
##            representation(
##                           block="guiWidgetORgWidgetRGtkORRGtkObject",
##                           widget="guiWidgetORgWidgetRGtkORRGtkObject",
##                           toolkit="guiWidgetsToolkit"
##                           ),
##            contains="gWidgetRGtk",
##            )
##   setClass("gContainerRGtk",
##            representation(
##                           block="guiWidgetORgWidgetRGtkORRGtkObject",
##                           widget="guiWidgetORgWidgetRGtkORRGtkObject",
##                           toolkit="guiWidgetsToolkit"
##                           ),
##            contains="gWidgetRGtk",
##            )
  
## }




##################################################
### Common methods.    Specific to a class are put into the file for that class

## we have two definitions. For instance, "svalue" and ".svalue". The "svalue" method dispatches on the object to the .svalue method. This allows us to use svalue instead of .svalue when defining the methods/constructors inside this package.


setMethod("svalue",signature(obj="gWidgetRGtk"),
          function(obj, index=NULL, drop=NULL, ...) {
            .svalue(obj, obj@toolkit, index=index, drop=drop, ...)
          })



## svalue
## need method for character
setMethod("svalue",signature(obj="character"),
          function(obj, index=NULL, drop=NULL, ...)  {
            ifelse(length(obj) == 1,
                   return(getObjectFromString(obj)),
                   return(obj)
                   )
          })
setMethod(".svalue",signature(toolkit = "guiWidgetsToolkitRGtk2", obj="character"),
          function(obj, toolkit, index=NULL, drop=NULL,  ...)  {
            ifelse(length(obj) == 1,
                   return(getObjectFromString(obj)),
                   return(NA)
                   )
          })

## svalue<- -- objec specific
setReplaceMethod("svalue",signature(obj="gWidgetRGtk"),
          function(obj, index=NULL, ...,value) {
            .svalue(obj, obj@toolkit, index=index, ...) <- value
            obj
          })

## [
setMethod("[",
          signature(x="gWidgetRGtk"),
          function(x,i,j,...,drop=TRUE) {
            
            return(.leftBracket(x, x@toolkit,i,j,...,drop=TRUE))
#            if(missing(i) && missing(j))
#              .leftBracket(x@widget, toolkit,,,...,drop=TRUE)
#            else if(missing(j))
#              .leftBracket(x@widget, toolkit,i,,...,drop=TRUE)
#            else
#              .leftBracket(x@widget, toolkit,i,j,...,drop=TRUE)
          })

## [<-
setReplaceMethod("[",signature(x="gWidgetRGtk"),
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
setMethod("size",signature(obj="gWidgetRGtk"),
          function(obj, ...) {
            return()
            .size(obj, obj@toolkit,...)
          })

setMethod(".size",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit, ...) {
            ## send to gWidgetRGtk2
            .size(obj@widget, toolkit, ...)
          })
setMethod(".size",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="RGtkObject"),
          function(obj, toolkit, ...) {
            ## get from SizeAllocation()
            val <- obj$GetAllocation()
            return(c(width=val$width, height=val$height))
          })

## not needed?
setMethod(".size",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="GtkWindow"),
          function(obj, toolkit, ...) {
            print("returning size information not yet implemented.")
          })

## size<-
setReplaceMethod("size",signature(obj="gWidgetRGtk"),
          function(obj, ..., value) {
            .size(obj, obj@toolkit,...) <- value
            return(obj)
          })

setReplaceMethod(".size", 
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
                 function(obj, toolkit, ..., value) {
                   if(length(value) >= 2) {
                     width <- value[1]; height <- value[2]
                   } else if(names(value) == "height") {
                     width <- -1; height <- value
                   } else {
                     width <- value; height <- -1
                   }
                   widget = obj@widget
                   widget$SetSizeRequest(width,height)
#                   widget$SetDefaultSize(width,height)

                   return(obj)
                 })

## visible
setMethod("visible",signature(obj="gWidgetRGtk"),
          function(obj, set=NULL, ...) {
            .visible(obj,obj@toolkit, set=set, ...)
          })

setMethod(".visible",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit, set=TRUE, ...) {
            widget <- getWidget(obj)
            if(is.null(set))
              widget['visible']
            else
              if(as.logical(set))
                widget$Show()
              else
                widget$Hide()

          })


## visible<-
setReplaceMethod("visible",signature(obj="gWidgetRGtk"),
          function(obj, ..., value) {
            .visible(obj, obj@toolkit, ...) <- value
            return(obj)
          })

setReplaceMethod(".visible",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
                 function(obj, toolkit, ..., value) {
                   visible(obj, value)
                   return(obj)
                 })

## isExtant -- can we see the window
setMethod("isExtant",signature(obj="gWidgetRGtk"),
          function(obj, ...) {
            .isExtant(obj,obj@toolkit, ...)
          })

setMethod(".isExtant",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit, ...) {
            widget = getWidget(obj)
            if(is.null(widget)) return(FALSE)
            !inherits(widget,"<invalid>") # test to see if destroyed
          })


## enabled -- not implemeneted, don't   know how to find sensitive. Would need to keep in
##            in the widget using tag or somesuch
setMethod("enabled",signature(obj="gWidgetRGtk"),
          function(obj, ...) {
            warning("enable not defined, try enabled<-()")
            return(NA)
            .enabled(obj, obj@toolkit,...)
          })

setMethod(".enabled",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit, ...) {
            widget <- getWidget(obj)
            if("sensitive" %in% names(widget))
              widget['sensitive']
            else
              TRUE
          })

setMethod(".enabled",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="GtkWindow"),
          function(obj, toolkit, ...) {
            print("returning enabled information not yet implemented.")
          })

## enabled<-
setReplaceMethod("enabled",signature(obj="gWidgetRGtk"),
          function(obj, ..., value) {
            .enabled(obj, obj@toolkit,...) <- value
            return(obj)
          })

setReplaceMethod(".enabled",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
                 function(obj, toolkit, ..., value) {
                   widget = getWidget(obj)
                   .enabled(widget, toolkit, ...) <- value
                   return(obj)
                 })
setReplaceMethod(".enabled",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="RGtkObject"),
                 function(obj, toolkit, ..., value) {
                   obj$SetSensitive(as.logical(value))
                   return(obj)
                 })


## editable -- can the widget be edited
setMethod("editable",signature(obj="gWidgetRGtk"),
          function(obj, ...) {
            return()
            .editable(obj, obj@toolkit,...)
          })

setMethod(".editable",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit, ...) {
            message("no default editable method")
          })

## editable<-
setReplaceMethod("editable",signature(obj="gWidgetRGtk"),
                 function(obj, ..., value) {
                   .editable(obj, obj@toolkit,...) <- value
                   return(obj)
                 })

setReplaceMethod(".editable", 
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
                 function(obj, toolkit, ..., value) {
                   message("no default editable<- method")
                   return(obj)
                 })

## focus
setMethod("focus",signature(obj="gWidgetRGtk"),
          function(obj, ...) {
            .focus(obj, obj@toolkit,...)
          })

setMethod(".focus",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit, ...) focus(obj) <- TRUE)

## focus<-
setReplaceMethod("focus",signature(obj="gWidgetRGtk"),
          function(obj, ..., value) {
            .focus(obj, obj@toolkit,...) <- value
            return(obj)
          })

setReplaceMethod(".focus",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit, ..., value) {
            focus(obj@widget, toolkit, ...) <- value
            return(obj)
          })

setReplaceMethod("focus",signature(obj="RGtkObject"),
          function(obj, ..., value) {
            .focus(obj, toolkit=guiToolkit("RGtk2"),...) <- value
            return(obj)
          })


## window
setReplaceMethod(".focus",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="GtkWindow"),
          function(obj, toolkit, ..., value) {
            value = as.logical(value)
            if(value)
              obj$GetWindow()$Raise()
            else
              obj$GetWindow()$Lower()
            return(obj)
          })
## other objects
setReplaceMethod(".focus",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="RGtkObject"),
          function(obj, toolkit, ..., value) {
            value = as.logical(value)
            if(value) {
              obj$GrabFocus()
              obj$GetParentWindow()$Raise()
            } else {
              obj$GetParentWindow()$Lower()
            }
            return(obj)
          })


## tooltip<-
setReplaceMethod("tooltip",signature(obj="gWidgetRGtk"),
          function(obj, ..., value) {
            .tooltip(obj, obj@toolkit,...) <- value
            return(obj)
          })

setReplaceMethod(".tooltip",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit, ..., value) {
            tooltip(obj@widget, toolkit, ...) <- value
            return(obj)
          })

setReplaceMethod("tooltip",signature(obj="RGtkObject"),
          function(obj, ..., value) {
            ## set the tip.
            obj$setTooltipText(paste(value, collapse="\n"))
            ## deprecated
            ## tooltipGroup <- try(gtkTooltips(), silent = TRUE)
            ## if(inherits(tooltipGroup, "try-error"))
            ##   return(obj)
            ## ## some widgets don't allow a tooltip (glabel, ...)
            ## ## right check is widget.flags()&gtk.NO_WINDOW
            ## try(tooltipGroup$setTip(obj, tip.text = value), silent=TRUE)
            return(obj)
          })

## default Widget
## defaultWidget
setMethod("defaultWidget",signature(obj="gWidgetRGtk"),
          function(obj, ...) {
            .defaultWidget(obj, obj@toolkit,...)
          })

setMethod(".defaultWidget",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit, ...)
          getWidget(obj)['has-default']
          )

## defaultWidget<-
setReplaceMethod("defaultWidget",signature(obj="gWidgetRGtk"),
                 function(obj, ..., value) {
                   .defaultWidget(obj, obj@toolkit,...) <- value
                   return(obj)
                 })

setReplaceMethod("defaultWidget",signature(obj="RGtkObject"),
          function(obj, ..., value) {
            .defaultWidget(obj, toolkit=guiToolkit("RGtk2"),...) <- value
            return(obj)
          })


setReplaceMethod(".defaultWidget",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit, ..., value) {
            widget <- getWidget(obj)
            .defaultWidget(widget, toolkit, ...) <- value
            return(obj)
          })

setReplaceMethod(".defaultWidget",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="RGtkObject"),
          function(obj, toolkit, ..., value) {
            value = as.logical(value)
            obj['can-default'] <- value
            obj$grabDefault()
#            obj['receives-default'] <- value
            return(obj)
          })


## fonts
.font.styles = list(
  family = c("normal","sans","serif","monospace"),
  style = c("normal","oblique","italic"),
  weight = c("ultra-light","light","normal","bold","ultra-bold","heavy"),
  colors = c("black","blue","red","green","brown","yellow","pink")
)  
## font sizes ## old defs for .PangoScale are no longer valid as of 10.4
fontSizes <- c(
               "xx-large"= PANGO_SCALE_XX_LARGE,
               "x-large" = PANGO_SCALE_X_LARGE,
               "large"   = PANGO_SCALE_LARGE,
               "medium"  = PANGO_SCALE_MEDIUM,
               "small"   = PANGO_SCALE_SMALL,
               "x-small" = PANGO_SCALE_X_SMALL,
               "xx-small" = PANGO_SCALE_XX_SMALL
               )

setMethod("font",signature(obj="gWidgetRGtk"),
          function(obj, ...) {
            warning("font() not defined. Set fonts with font<-")
            return()
            .font(obj, obj@toolkit,...)
          })

setMethod(".font",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="GtkWindow"),
          function(obj, toolkit, ...) {
            print("returning font information not yet implemented.")
          })
## font<-
setReplaceMethod("font",signature(obj="gWidgetRGtk"),
          function(obj, ..., value) {
            .font(obj, obj@toolkit,...) <- .fixFontMessUp(value)
            return(obj)
          })
setReplaceMethod(".font",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
                 function(obj, toolkit, ..., value) {
                   .font(obj@widget, toolkit, ...) <- value
                   return(obj)
                 })
setReplaceMethod(".font",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="RGtkObject"),
                 function(obj, toolkit, ..., value) {

                   ## value might be a vector, we use a list -- from .fixFontMessUp
                   if(!is.list(value)) {
                     tmp = value
                     value = list()
                     for(i in names(tmp)) value[[i]] = tmp[i]
                   }
                   
                   
                   string = ""

                   
                   ## do family, weight, style
                   for(i in c("family", "weight", "style")) {
                     if(!is.null(value[[i]])) {
                       x <- .font.styles[[i]]
                       ind <- charmatch(value[[i]], x)
                       if(!is.na(ind)) {
                         string <- paste(string, x[ind[1]], sep=" ")
                         if(i == "family")
                           string <- paste(string,",", sep="")
                       }
                     }
                   }
                   
                   ## size can be integer or name -- relative to 12pt
                   
                   if(!is.null(value$size)) {
                     ## is it numeric or character?
                     warn <- getOption("warn"); options(warn=2) # hack to avoid warning -- we want an error here
                     out <- try(as.integer(value[['size']]), silent=TRUE)
                     options(warn=warn)
                     if(!inherits(out, "try-error"))
                       string <- Paste(string," ",out)
                     else if (!is.na(ind <- charmatch(value[['size']], names(fontSizes)))) # fuzzy match?
                       string <- Paste(string, " ", paste(ceiling(12*fontSizes[ind[1]]),"px", sep=""))
                   }
                   string <- gsub(",$","",string) # strip , if present
                   
                   if(string != "") {
                     fontDescr = pangoFontDescriptionFromString(string)
                     obj$ModifyFont(fontDescr)
                   }

                   ## colors
                   if(!is.null(value$color))
                     obj$modifyFg(GtkStateType[1], value[['color']])

                   
                   return(obj)
                 })
## tag, tag<-
setMethod("tag",signature(obj="gWidgetRGtk"),
          function(obj,i,drop=TRUE, ...) {
            if(missing(drop)) drop <- TRUE
            .tag(obj, obj@toolkit,i, drop=drop,...)
          })
## dispatch in *this* toolkit, not present in obj
setMethod("tag",signature(obj="RGtkObject"),
          function(obj,i,drop=TRUE, ...) {
            if(missing(drop)) drop <- TRUE            
            .tag(obj, guiToolkit("RGtk2"),i, drop=drop,...)
          })

setMethod(".tag", signature(toolkit="guiWidgetsToolkitRGtk2",obj="guiWidget"),
          function(obj, toolkit, i, drop=TRUE, ...) {
            if(missing(i)) i = NULL
            if(missing(drop)) drop <- TRUE                        
            .tag(obj@widget,toolkit,  i, drop=drop,  ...)
          })
setMethod(".tag", signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit, i, drop=TRUE, ...) {
            if(missing(i)) i = NULL
            if(missing(drop)) drop <- TRUE                                    
            .tag( obj@block,toolkit,  i, drop=drop,  ...)
          })
setMethod(".tag", signature(toolkit="guiWidgetsToolkitRGtk2",obj="RGtkObject"),
          function(obj, toolkit, i, drop=TRUE, ...) {
            if(missing(i)) i = NULL

            lst = obj$GetData(".tagKey")
            if(is.null(i))
              return(lst)
            if(drop) {
              if(length(i) == 1)
                return(lst[[i]])
              else
                return(sapply(i, function(j) lst[j]))
            } else {
              return(sapply(i, function(j) lst[j]))
            }
          })

## tag <-
setReplaceMethod("tag",signature(obj="gWidgetRGtk"),
          function(obj, i, replace=TRUE, ..., value) {
            .tag(obj, obj@toolkit,i,replace, ...) <- value
            return(obj)
          })
## dispatch in *this* toolkit, not present in obj
setReplaceMethod("tag",signature(obj="RGtkObject"),
          function(obj,i, replace=TRUE, ..., value) {
            .tag(obj, guiToolkit("RGtk2"),i, replace, ...) <- value
            return(obj)
          })

## objects can be in many different flavors: guiWIdget, gWidgetRGtk2, RGtkObject
setReplaceMethod(".tag", signature(toolkit="guiWidgetsToolkitRGtk2",obj="guiWidget"),
          function(obj, toolkit, i, replace=TRUE, ..., value) {
            if(missing(i)) i = NULL
            .tag(obj@widget,toolkit,  i, replace, ...) <- value
            return(obj)
          })

setReplaceMethod(".tag", signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit, i, replace=TRUE, ..., value) {
            if(missing(i)) i = NULL
            .tag( obj@block, toolkit,  i, replace, ...) <- value
            return(obj)
          })

setReplaceMethod(".tag", signature(toolkit="guiWidgetsToolkitRGtk2",obj="RGtkObject"),
          function(obj, toolkit, i, replace=TRUE, ..., value) {
            if(missing(i) || is.null(i)) {
              warning("Need to specify a key to the 'i' argument of tag<-")
            } else {
              allData = obj$GetData(".tagKey")
              if(is.null(allData)) allData = list()
              
              if(as.logical(replace)) {
                allData[[i]] <- value
              } else {
                allData[[i]] <- c(allData[[i]], value)
              }
              obj$SetData(".tagKey", allData)
            }
            return(obj)
          })


##################################################
## id -- define for "ANY" as well
setMethod("id",signature(obj="gWidgetRGtk"),
          function(obj, ...) {
            tag(obj,".gtkID")
          })
setMethod("id",signature(obj="RGtkObject"),
          function(obj, ...) {
            tag(obj, ".gtkID", ...)
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


setMethod(".id", signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit, ...) {
            tag(obj,".gtkID", ...)
          })
setMethod(".id", signature(toolkit="guiWidgetsToolkitRGtk2",obj="RGtkObject"),
          function(obj, toolkit,  ...) {
            return(tag(obj,".gtkID"))
          })


## id<-
setReplaceMethod("id",signature(obj="gWidgetRGtk"),
          function(obj, ..., value) {
            tag(obj,".gtkID", ...) <- value
            return(obj)
          })
## dispatch in *this* toolkit, not present in obj
setReplaceMethod("id",signature(obj="RGtkObject"),
          function(obj, ..., value) {
            tag(obj, ".gtkID", ...) <- value
            return(obj)
          })
setReplaceMethod("id",signature(obj="ANY"),
          function(obj, ..., value) {
            attr(obj,"id") <- value
            return(obj)
          })


## we need a .id to handle dispatch from guiWidgets, otherwise, we use id()
setReplaceMethod(".id", signature(toolkit="guiWidgetsToolkitRGtk2",
                                  obj="gWidgetRGtk"),
          function(obj, toolkit, ..., value) {
            id(obj, ...) <- value
            return(obj)
          })



## add method is biggie
## we have several levels of classes here guiWidget -- gWidgetRGkt -- RGtkObject, when
## we get down to that level we can finally add
setMethod("add",signature(obj="gWidgetRGtk"),
          function(obj, value, ...) {
            .add(obj, obj@toolkit,value,...)
          })
setMethod(".add",
          signature(toolkit="guiWidgetsToolkitRGtk2",
                    obj="guiWidget", value="ANY"),
          function(obj, toolkit, value, ...) {
            cat(gettext("Can't add without a value\n"))
          })
setMethod(".add",
          signature(toolkit="guiWidgetsToolkitRGtk2",
                    obj="gWidgetRGtk", value="try-error"),
          function(obj, toolkit, value, ...) {
            gmessage(paste("Error:",value), title="Error adding oject",
                     icon="error")
          })
## pushdonw

setMethod(".add",
          signature(toolkit="guiWidgetsToolkitRGtk2",
                    obj="guiWidget", value="guiWidgetORgWidgetRGtkORRGtkObject"),
          function(obj, toolkit, value, ...) {
            .add(obj@widget, toolkit, value, ...)
          })

## for gWindow
setMethod(".add",
          signature(toolkit="guiWidgetsToolkitRGtk2",
                    obj="gContainerRGtk", value="guiWidget"),
          function(obj, toolkit, value, ...) {
            .add(obj, toolkit, value@widget, ...)
          })

setMethod(".add",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gContainerRGtk", value="gWidgetRGtk"),
          function(obj, toolkit, value, ...) {
            .add(obj, toolkit, value@block, ...)
          })




## addSPring, addSpace
setMethod("addSpring",signature(obj="gWidgetRGtk"),
          function(obj, ...) {
            .addSpring(obj, obj@toolkit,...)
          })

setMethod(".addSpring",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gContainerRGtk"),
          function(obj, toolkit, ...) {
            obj@widget$PackStart(gtkHBoxNew(),TRUE,TRUE,0) # expand and fill set to TRUE
          })

setMethod("addSpace",signature(obj="gWidgetRGtk"),
          function(obj, value, ...) {
            .addSpace(obj,obj@toolkit,value,...)
          })

setMethod(".addSpace",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gContainerRGtk"),
          function(obj, toolkit, value, ...) {
            theArgs = list(...)
            horizontal = ifelse(is.null(theArgs$horizontal),
              TRUE,
              as.logical(theArgs$horizontal))

            if(horizontal) {
              tmp = ggroup(); size(tmp) <- c(value,1)
            } else {
              tmp = ggroup(); size(tmp) <- c(1, value)
            }
            add(obj, tmp)
          })

## delete -- get down to two RGtkObjects
setMethod("delete",signature(obj="gWidgetRGtk"),
          function(obj, widget, ...) {
            .delete(obj, obj@toolkit,widget,...)
          })

## push down to RGtk vs RGtk. Can be 9 possiblities!
setMethod(".delete",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gContainerRGtk",widget="guiWidget"),
          function(obj, toolkit, widget, ...) {
            .delete(obj, toolkit, widget@widget, ...)
          })
setMethod(".delete",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gContainerRGtk",widget="gWidgetRGtk"),
          function(obj, toolkit, widget, ...) {
            .delete(obj@widget, toolkit, widget, ...)
          })

setMethod(".delete",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="RGtkObject",widget="gWidgetRGtk"),
          function(obj, toolkit, widget, ...) {
            .delete(obj, toolkit, getBlock(widget), ...)
          })
setMethod(".delete",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="RGtkObject",widget="guiWidget"),
          function(obj, toolkit, widget, ...) {
            .delete(obj, toolkit, widget@widget, ...)
          })
setMethod(".delete",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk",widget="RGtkObject"),
          function(obj, toolkit, widget, ...) {
            .delete(obj@widget, toolkit, widget, ...)
          })
setMethod(".delete",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="RGtkObject",widget="RGtkObject"),
          function(obj, toolkit, widget, ...) {
            ## Remove after checking
            if(!inherits(obj, "<invalid>") &&
               !inherits(widget, "<invalid>") &&
               widget$getParent() == obj) {
              obj$Remove(widget)
            }
            return(TRUE)
          })

## dispose -- delete the parent window, or something else
setMethod("dispose",signature(obj="gWidgetRGtk"),
          function(obj, ...) {
            .dispose(obj, obj@toolkit,...)
          })

setMethod(".dispose",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit, ...) {
            widget = getWidget(obj)
            if(inherits(widget,"GtkWindow")) {
              widget$Destroy()
              return(TRUE)
            } else {
              widget = widget$GetParentWindow()
              if(inherits(widget,"<invalid>"))
                return(FALSE)
              else
                widget$Destroy()
              return(TRUE)
            }
          })




## update
setMethod("update",signature(object="gWidgetRGtk"),
          function(object, ...) {
            .update(object, object@toolkit, ...)
          })

setMethod(".update",
          signature(toolkit="guiWidgetsToolkitRGtk2",object="gComponentRGtk"),
          function(object, toolkit, ...) {
            object@widget$QueueDraw()
          })

##
##
##################################################


##################################################
## handlers
##
## basic handler for adding with a signal. Now exported.
setGeneric("addhandler", function(obj, signal, handler, action=NULL, ...)
           standardGeneric("addhandler"))
setGeneric("addHandler", function(obj, signal, handler, action=NULL, ...)
           standardGeneric("addHandler"))

setMethod("addhandler",signature(obj="guiWidget"),
          function(obj, signal, handler, action=NULL, ...) {
            .addHandler(obj@widget, obj@toolkit, signal, handler, action, ...)
          })
setMethod("addhandler",signature(obj="gWidgetRGtk"),
          function(obj, signal, handler, action=NULL, ...) {
            .addHandler(obj, obj@toolkit, signal, handler, action, ...)
          })
setMethod("addhandler",signature(obj="RGtkObject"),
          function(obj, signal, handler, action=NULL, ...) {
            .addHandler(obj, guiToolkit("RGtk2"), signal, handler, action, ...)
          })
setMethod("addHandler",signature(obj="gWidgetRGtk"),
          function(obj, signal, handler, action=NULL, ...) {
            .addHandler(obj@widget, obj@toolkit, signal, handler, action, ...)
          })

setMethod(".addhandler",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit, signal, handler, action=NULL, ...) {
            .addHandler(obj, obj@toolkit, signal, handler, action, ...)
          })

## method for dispatch
setGeneric(".addHandler",
           function(obj, toolkit,
                  signal, handler, action=NULL, ...)
           standardGeneric(".addHandler"))


setMethod(".addHandler",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="guiWidget"),
          function(obj, toolkit,
                   signal, handler, action=NULL, ...) {
            .addHandler(obj@widget, toolkit, signal, handler, action, ...)
          })

setMethod(".addHandler",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit,
                   signal, handler, action=NULL, ...) {
            
            ## need to return logical when an event, not always,
            ## but gives trouble if not
            modifyHandler = function(...) {
              val <- handler(...)
              if(!is.logical(val))
                return(TRUE)
              else
                return(val)
            }

            theArgs = list(...)

            ## fix value passed into gWidgets handlers
            h <- list()
            h$obj <- obj
            if(!is.null(theArgs$actualobj)) {
              h$obj <- theArgs$actualobj
              theArgs$actualobj <- NULL
            }
            if(length(theArgs)) 
              for(i in names(theArgs))
                h[[i]] <- theArgs[[i]]

            h$action <- action

            callbackID <- gtktry(connectSignal(getWidget(obj), ### issue: getWidget(obj),
                                            signal=signal,
                                            f=modifyHandler,
                                            data=h,
                                            user.data.first = TRUE,
                                            after = FALSE), silent=FALSE)
            if(inherits(callbackID,"try-error")) {
              gwCat(sprintf("Couldn't add signal %s for object of class %s",
                      signal, class(obj)[1]),"\n")
              return(NA)
            } else {
              ## now put handler into object
              handler.ID = tag(obj, "handler.id")
              if(is.null(handler.ID))
                handler.ID =list()
              handler.ID[[length(handler.ID)+1]] = callbackID
              tag(obj, "handler.id") <- handler.ID
              
              ##
              ##            addhandlerdestroy(obj, handler=function(h,...)
              ##                              removehandler(h$obj,h$action),
              ##                              action=ID)
              ## return ID
              invisible(callbackID)
            }
          })

setMethod(".addHandler",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="RGtkObject"),
          function(obj, toolkit, signal, handler, action=NULL, ..., after=FALSE) {

            theArgs = list(...)

            ## need to return FALSE to propogate to next handler on events
            modifyHandler = function(h,...) {
              handler(h,...)
              return(FALSE)
            }

            ## fix value passed into gWidgets handlers
            h <- list()
            h$obj <- obj
            if(!is.null(theArgs$actualobj)) {
              h$obj <- theArgs$actualobj
              theArgs$actualobj <- NULL
            }
            ## pass in extra values, eg addHandlerClicked(obj, f=FUN, extra=value) will h$extra key

            if(length(theArgs)) 
              for(i in names(theArgs))
                h[[i]] <- theArgs[[i]]
              
            h$action <- action
            
            
            callbackID <- gtktry(gSignalConnect(obj,
                                            signal=signal,
                                            f=modifyHandler,
                                            data=h,
                                            user.data.first = TRUE,
                                            after = after),
                                 silent=TRUE)
            ## can't' stuff in handler IDS
            if(inherits(callbackID,"try-error")) {
              gwCat(sprintf("Couldn't connect signal: %s for object of class %s\n",
                    signal, class(obj)[1]))
              return(NA)
            } else {
              ## store ID into list
              lst <- obj$getData("handler.id")
              if(is.null(lst)) lst <- list()
              lst <- c(lst,callbackID)
              obj$setData("handler.id", lst)
              
              invisible(callbackID)
            }
          })

  ## removehandler
setMethod("removehandler", signature("gWidgetRGtk"),
          function(obj, ID=NULL, ...) {
            .removehandler(obj, obj@toolkit, ID, ...)
          })
setMethod("removehandler", signature("RGtkObject"),
          function(obj, ID=NULL, ...) {
            .removehandler(obj, guiToolkit("RGtk2"), ID, ...)
          })

### JV: Need to consolidate this and the next, The difference is the callbackIDS?
setMethod(".removehandler",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit, ID=NULL, ...) {
            if(missing(ID))
              callbackIDs = .tag(obj,toolkit,"handler.id")
            else
              callbackIDs = ID

            ## timeout handler?
            if(class(callbackIDs) == "GTimeoutId") {
              out <- sapply(callbackIDs, function(i) {
                class(i) <- "GTimeoutId" # sapply strips attributes!
                gSourceRemove(i)
              })
              return(out)
            }

            
            if(!is.null(callbackIDs)) {
              if(!is.list(callbackIDs)) {
                callbackIDs = list(callbackIDs)
              }
              
              widget = obj@widget
              retval = logical(length(callbackIDs))
              for(i in 1:length(callbackIDs)) {
                if(is.list(callbackIDs[[i]])) # recurse if a list
                  for(i in callbackIDs[[i]]) .removehandler(obj, toolkit, i)
                isCallbackID = gtktry(checkPtrType(callbackIDs[[i]],"CallbackID"),silent=TRUE)
                if(!inherits(isCallbackID,"try-error")) {
                  retval[i] = gtktry(gSignalHandlerDisconnect(widget, callbackIDs[[1]]), silent=TRUE)

#                  retval[i] = gtktry(gtkObjectDisconnectCallbackHack(widget, callbackIDs[[i]]),
#                          silent=TRUE)
                } else {
                  gwCat("DEBUG: ID not of callbackID\n")
                  retval[i] = FALSE
                }
              }
              for(i in rev(which(retval==TRUE)))
                callbackIDs[[i]] <- NULL
              .tag(obj,toolkit, "handler.id", replace=FALSE) <- callbackIDs
              return(retval)
            } else {
              return(FALSE)
            }
          })

## for RGtkObject
setMethod(".removehandler",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="RGtkObject"),
          function(obj, toolkit, ID=NULL, ...) {
            if(missing(ID))
              callbackIDs = .tag(obj,toolkit,"handler.id")
            else
              callbackIDs = ID
            
            if(!is.null(callbackIDs)) {
              if(!is.list(callbackIDs)) {
                callbackIDs = list(callbackIDs)
              }
              widget = obj
              retval = c()
              for(i in 1:length(callbackIDs)) {
                if(is.list(callbackIDs[[i]])) # recurse if a list
                  for(i in callbackIDs[[i]]) .removehandler(obj, toolkit, i)
                isCallbackID = gtktry(checkPtrType(callbackIDs[[i]],"CallbackID"),silent=TRUE)
                if(!inherits(isCallbackID,"try-error")) {
                  retval[i] = widget$disconnectCallback(callbackIDs[[i]]) #gtkObjectDisconnectCallbackHack(widget, callbackIDs[[i]])
                } else {
                  gwCat("DEBUG: ID not of callbackID\n")
                  retval[i] = FALSE
                }
              }

              for(i in rev(which(retval==TRUE)))
                callbackIDs[[i]] <- NULL
              .tag(obj,toolkit, "handler.id", replace=FALSE) <- callbackIDs
              return(retval)
            } else {
              return(FALSE)
            }
          })


## blockhandler
setMethod("blockhandler", signature("gWidgetRGtk"),
          function(obj, ID=NULL, ...) {
            .blockhandler(obj, obj@toolkit, ID, ...)
          })
setMethod("blockhandler", signature("RGtkObject"),
          function(obj, ID=NULL, ...) {
            .blockhandler(obj, guiToolkit("RGtk2"), ID, ...)
          })

## caps
setMethod("blockHandler", signature("gWidgetRGtk"),
          function(obj, ID=NULL, ...) {
            .blockhandler(obj, obj@toolkit, ID, ...)
          })
setMethod("blockHandler", signature("RGtkObject"),
          function(obj, ID=NULL, ...) {
            .blockhandler(obj, guiToolkit("RGtk2"), ID, ...)
          })



setMethod(".blockhandler",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit, ID=NULL, ...) {
            .blockhandler(getWidget(obj),toolkit,ID,...)
          })

setMethod(".blockhandler",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="RGtkObject"),
          function(obj, toolkit, ID=NULL, ...) {
            if(is.null(ID))
              ID <- tag(obj,"handler.id")
            lapply(ID, function(i)
                     gSignalHandlerBlock(obj,i))
            return()
          })

## unblock handler
setMethod("unblockhandler", signature("gWidgetRGtk"),
          function(obj, ID=NULL, ...) {
            .unblockhandler(obj, obj@toolkit, ID, ...)
          })
setMethod("unblockhandler", signature("RGtkObject"),
          function(obj, ID=NULL, ...) {
            .unblockhandler(obj, guiToolkit("RGtk2"), ID, ...)
          })
## camelcase
setMethod("unblockHandler", signature("gWidgetRGtk"),
          function(obj, ID=NULL, ...) {
            .unblockhandler(obj, obj@toolkit, ID, ...)
          })
setMethod("unblockHandler", signature("RGtkObject"),
          function(obj, ID=NULL, ...) {
            .unblockhandler(obj, guiToolkit("RGtk2"), ID, ...)
          })

setMethod(".unblockhandler",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit, ID=NULL, ...) {
            .unblockhandler(getWidget(obj),toolkit,ID,...)
          })

setMethod(".unblockhandler",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="RGtkObject"),
          function(obj, toolkit, ID=NULL, ...) {
            if(is.null(ID))
              ID <- tag(obj,"handler.id")

            sapply(ID, function(i)
                   gSignalHandlerUnblock(obj,i))
            return()
          })


## addhandlerchanged
setMethod("addhandlerchanged",signature(obj="gWidgetRGtk"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerchanged(obj, obj@toolkit, handler, action, ...)
          })
setMethod("addhandlerchanged",signature(obj="RGtkObject"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerchanged(obj, guiToolkit("RGtk2"), handler, action, ...)
          })
setMethod("addhandlerchanged",signature(obj="ANY"),
          function(obj, handler=NULL, action=NULL, ...) {
            warning("No method addhandlerchanged for object of class",class(obj),"\n")
          })
## caps
setMethod("addHandlerChanged",signature(obj="gWidgetRGtk"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerchanged(obj, obj@toolkit, handler, action, ...)
          })
setMethod("addHandlerChanged",signature(obj="RGtkObject"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerchanged(obj, guiToolkit("RGtk2"), handler, action, ...)
          })

setMethod(".addhandlerchanged",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit,
                   handler, action=NULL, ...) {
            .addHandler(obj, toolkit, signal="changed",
                        handler=handler, action=action, ...)
          })


## expose: expose-event or realize
setMethod("addhandlerexpose",signature(obj="gWidgetRGtk"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerexpose(obj,obj@toolkit,handler, action, ...)
          })
setMethod("addhandlerexpose",signature(obj="RGtkObject"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerexpose(obj, guiToolkit("RGtk2"), handler, action, ...)
          })

setMethod(".addhandlerexpose",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit,
                   handler, action=NULL, ...) {
            .addHandler(obj, toolkit, signal="expose-event",
                        handler=handler, action=action, ...)
          })

setMethod(".addhandlerexpose",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gComponentRGtk"),
          function(obj, toolkit,
                   handler, action=NULL, ...) {
            .addHandler(obj,toolkit, signal="realize",
                        handler=handler, action=action, ...)
          })

## unrealize: unrealize
setMethod("addhandlerunrealize",signature(obj="gWidgetRGtk"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerunrealize(obj, obj@toolkit,handler, action, ...)
          })
setMethod("addhandlerunrealize",signature(obj="RGtkObject"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerunrealize(obj, guiToolkit("RGtk2"),handler, action, ...)
          })

setMethod(".addhandlerunrealize",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit,
                   handler, action=NULL, ...) {
            .addHandler(obj, toolkit, signal="unrealize",
                        handler=handler, action=action, ...)
          })

## destroy: destroy
setMethod("addhandlerdestroy",signature(obj="gWidgetRGtk"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerdestroy(obj, obj@toolkit,handler, action, ...)
          })
setMethod("addhandlerdestroy",signature(obj="RGtkObject"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerdestroy(obj, guiToolkit("RGtk2"),handler, action, ...)
          })

setMethod(".addhandlerdestroy",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit,
                   handler, action=NULL, ...) {
            ## See
            ## http://www.moeraki.com/pygtktutorial/pygtk2tutorial/sec-SteppingThroughHelloWorld.html
            ## for difference between "destroy" and "delete-event"
            .addHandler(obj, toolkit, signal="destroy",
                        handler=handler, action=action, ...)
          })

## keystroke: changed
setMethod("addhandlerkeystroke",signature(obj="gWidgetRGtk"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerkeystroke(obj, obj@toolkit,handler, action, ...)
          })
setMethod("addhandlerkeystroke",signature(obj="RGtkObject"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerkeystroke(obj, guiToolkit("RGtk2"),handler, action, ...)
          })

setMethod(".addhandlerkeystroke",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit,
                   handler, action=NULL, ...) {
            .addHandler(obj, toolkit, signal="changed",
                        handler=handler, action=action, ...)
          })

## clicked: clicked
setMethod("addhandlerclicked",signature(obj="gWidgetRGtk"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerclicked(obj, obj@toolkit,handler, action, ...)
          })
setMethod("addhandlerclicked",signature(obj="RGtkObject"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerclicked(obj, guiToolkit("RGtk2"),handler, action, ...)
          })
## caps
setMethod("addHandlerClicked",signature(obj="gWidgetRGtk"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerclicked(obj, obj@toolkit,handler, action, ...)
          })
setMethod("addHandlerClicked",signature(obj="RGtkObject"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerclicked(obj, guiToolkit("RGtk2"),handler, action, ...)
          })

setMethod(".addhandlerclicked",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit,
                   handler, action=NULL, ...) {
            .addHandler(obj, toolkit, signal="clicked",
                        handler=handler, action=action, ...)
          })

setMethod(".addhandlerclicked",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="RGtkObject"),
          function(obj, toolkit,
                   handler, action=NULL, ...) {
            .addHandler(obj, toolkit, signal="clicked",
                        handler=handler, action=action, ...)
          })


## doubleclick: no default
setMethod("addhandlerdoubleclick",signature(obj="gWidgetRGtk"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerdoubleclick(obj,obj@toolkit,handler, action, ...)
          })
setMethod("addhandlerdoubleclick",signature(obj="RGtkObject"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerdoubleclick(obj,guiToolkit("RGtk2"),handler, action, ...)
          })
## caps
setMethod("addHandlerDoubleclick",signature(obj="gWidgetRGtk"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerdoubleclick(obj, obj@toolkit,handler, action, ...)
          })
setMethod("addHandlerDoubleclick",signature(obj="RGtkObject"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerdoubleclick(obj, guiToolkit("RGtk2"),handler, action, ...)
          })


setMethod(".addhandlerdoubleclick",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit,
                   handler, action=NULL, ...) {
            warning("No default handler for double click")
          })


## rightclick: button-press-event -- handle separately
setMethod("addhandlerrightclick",signature(obj="gWidgetRGtk"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerrightclick(obj,obj@toolkit,handler, action, ...)
          })
setMethod("addhandlerrightclick",signature(obj="RGtkObject"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerrightclick(obj,guiToolkit("RGtk2"),handler, action, ...)
          })
## caps
setMethod("addHandlerRightclick",signature(obj="gWidgetRGtk"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerrightclick(obj, obj@toolkit,handler, action, ...)
          })
setMethod("addHandlerRightclick",signature(obj="RGtkObject"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerrightclick(obj, guiToolkit("RGtk2"),handler, action, ...)
          })


            
## use actualobj=obj to pass in a different obj to h$obj
         
## use actualobj=obj to pass in a different obj to h$obj
setMethod(".addhandlerrightclick",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit,
                   handler, action=NULL, ...) {
            .addhandlerrightclick(getWidget(obj), toolkit, handler, action, actualobj=obj,...)
          })

setMethod(".addhandlerrightclick",
#          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="RGtkObject"),
          function(obj, toolkit,
                   handler, action=NULL, ..., after=FALSE) {
            theArgs = list(...)

            ## fix value passed into gWidgets handlers
            h <- list()
            h$handler <- handler
            h$obj <- obj
            if(!is.null(theArgs$actualobj)) {
              h$obj <- theArgs$actualobj
              theArgs$actualobj <- NULL
            }
            ## pass in extra values, eg addHandlerClicked(obj, f=FUN, extra=value) will h$extra key

            if(length(theArgs)) 
              for(i in names(theArgs))
                h[[i]] <- theArgs[[i]]
              
            h$action <- action
            
            gtktry(connectSignal(getWidget(obj),
                              signal = "button-press-event",
                              f = function(h, w, eventButton,...) {
                                if(isRightMouseClick(eventButton)) {
                                  h$handler(h,w, eventButton, ...)
                                }
                                return(FALSE)         # stop propagation
                              },
                              data = h,
                              user.data.first = TRUE,
                              after = after
                              ),
                silent=TRUE)
          })


### Column click things

## click: no default
setMethod("addhandlercolumnclicked",signature(obj="gWidgetRGtk"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlercolumnclicked(obj,obj@toolkit,handler, action, ...)
          })
setMethod("addhandlercolumnclicked",signature(obj="RGtkObject"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlercolumnclicked(obj,guiToolkit("RGtk2"),handler, action, ...)
          })
setMethod("addHandlerColumnClicked",signature(obj="gWidgetRGtk"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlercolumnclicked(obj, obj@toolkit,handler, action, ...)
          })
setMethod("addHandlerColumnClicked",signature(obj="RGtkObject"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlercolumnclicked(obj, guiToolkit("RGtk2"),handler, action, ...)
          })



## doubleclick: no default
setMethod("addhandlercolumndoubleclick",signature(obj="gWidgetRGtk"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlercolumndoubleclick(obj,obj@toolkit,handler, action, ...)
          })
setMethod("addhandlercolumndoubleclick",signature(obj="RGtkObject"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlercolumndoubleclick(obj,guiToolkit("RGtk2"),handler, action, ...)
          })
setMethod("addHandlerColumnDoubleclick",signature(obj="gWidgetRGtk"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlercolumndoubleclick(obj, obj@toolkit,handler, action, ...)
          })
setMethod("addHandlerColumnDoubleclick",signature(obj="RGtkObject"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlercolumndoubleclick(obj, guiToolkit("RGtk2"),handler, action, ...)
          })



## rightclick
setMethod("addhandlercolumnrightclick",signature(obj="gWidgetRGtk"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlercolumnrightclick(obj,obj@toolkit,handler, action, ...)
          })
setMethod("addhandlercolumnrightclick",signature(obj="RGtkObject"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlercolumnrightclick(obj,guiToolkit("RGtk2"),handler, action, ...)
          })
setMethod("addHandlerColumnRightclick",signature(obj="gWidgetRGtk"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlercolumnrightclick(obj, obj@toolkit,handler, action, ...)
          })
setMethod("addHandlerColumnRightclick",signature(obj="RGtkObject"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlercolumnrightclick(obj, guiToolkit("RGtk2"),handler, action, ...)
          })




## focus -- on focus call this
setMethod("addhandlerfocus",signature(obj="gWidgetRGtk"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerfocus(obj,obj@toolkit,handler, action, ...)
          })
setMethod("addhandlerfocus",signature(obj="RGtkObject"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerfocus(obj,guiToolkit("RGtk2"),handler, action, ...)
          })

setMethod(".addhandlerfocus",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit,
                   handler, action=NULL, ...) {
            .addHandler(obj,toolkit,
                        signal="focus-in-event",
                        handler, action, ...)
          })



## blur -- leave focus
setMethod("addhandlerblur",signature(obj="gWidgetRGtk"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerblur(obj,obj@toolkit,handler, action, ...)
          })
setMethod("addhandlerblur",signature(obj="RGtkObject"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlerblur(obj,guiToolkit("RGtk2"),handler, action, ...)
          })

setMethod(".addhandlerblur",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit,
                   handler, action=NULL, ...) {
            ## blur handler should return FALSE
            f <- function(h,...) {
              handler(h,...)
              return(FALSE)
            }
            .addHandler(obj,toolkit,
                        signal="focus-out-event",
                        f, action, ...)
          })




##
## mousemotion -- like mouseover
setMethod("addhandlermousemotion",signature(obj="gWidgetRGtk"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlermousemotion(obj,obj@toolkit,handler, action, ...)
          })
setMethod("addhandlermousemotion",signature(obj="RGtkObject"),
          function(obj, handler=NULL, action=NULL, ...) {
            .addhandlermousemotion(obj,guiToolkit("RGtk2"),handler, action, ...)
          })

setMethod(".addhandlermousemotion",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit,
                   handler, action=NULL, ...) {
            .addHandler(obj,guiToolkit("RGtk2"),
                        signal="enter-notify-event",
                        handler, action, ...)
          })



## idle
setMethod("addhandleridle",signature(obj="gWidgetRGtk"),
          function(obj, handler=NULL, action=NULL, interval=1000, ...) {
            .addhandleridle(obj, obj@toolkit,
                            handler=handler, action=action, interval=interval, ...)
          })
setMethod("addhandleridle",signature(obj="RGtkObject"),
          function(obj, handler=NULL, action=NULL, interval=1000, ...) {
            .addhandleridle(obj, guiToolkit("RGtk2"),
                            handler=handler, action=action, interval=interval, ...)
          })

setMethod(".addhandleridle",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit,
                   handler=NULL, action=NULL, interval=1000, ...) {

            idlehandler = function(h,...) {
              if(!is.null(h$handler) && is.function(h$handler))
                h$handler(h,...)
              invisible(TRUE)
            }
            
##            ID = gtkAddTimeout(
            ID = gTimeoutAdd(
              interval,
              idlehandler,
              data=list(obj=obj, action=action, handler=handler)
              )
            
            ## tidy up when done
            .addhandlerdestroy(obj,toolkit, handler=function(h,...) {
#              gtkRemoveTimeout(h$action)
              gSourceRemove(h$action)
            },action=ID)
            
            invisible(ID)
            
          })


## addpopumenu
setMethod("addpopupmenu",signature(obj="gWidgetRGtk"),
          function(obj, menulist, action=NULL, ...) {
            .addpopupmenu(obj, obj@toolkit,menulist, action, ...)
          })
setMethod("addpopupmenu",signature(obj="RGtkObject"),
          function(obj, menulist, action=NULL, ...) {
            .addpopupmenu(obj, guiToolkit("RGtk2"), menulist, action, ...)
          })


## this does not get exported
addPopupMenuWithSignal = function(obj, toolkit,  menulist, action=NULL, signal="button-press-event", ...) {
  theArgs = list(...)                      
              
  f = function(h, ...) {
    mb = gmenu(h$action, popup = TRUE)
    event = gdkEventNew(GdkEventType["button-press"])
    mb = tag(mb,"mb")                   # the real menu bar
    gtkMenuPopupHack(mb, button = event$GetButton(),
                     activate.time=event$GetTime()
                     )
  }
  ## .addhandler not exported
  callbackID = .addHandler(obj,toolkit, signal = signal,handler=f, action=menulist)
  invisible(callbackID)
}

add3rdMousePopupMenuWithSignal = function(obj, toolkit,  menulist, action=NULL, signal="button-press-event", ...) {
  f = function(h, widget, event,...) {
    ## Mac use ctrl - button 1
    if(isRightMouseClick(event)) {
      mb = gmenu(h$action$menulist, popup = TRUE, action=h$action$passedaction)
      mb = tag(mb,"mb")                 # actual widget
      gtkMenuPopupHack(mb,button = event$GetButton(),
                       activate.time=event$GetTime()
                       )
      return(FALSE)
    } else {
      return(FALSE)
    }
  }
  callbackID = .addHandler(obj,toolkit, signal = "button-press-event",handler=f, action=list(menulist=menulist, passedaction=action))
  invisible(callbackID)
}


  
### need to deal with action
setMethod(".addpopupmenu",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit, menulist, action=NULL, ...) {
            addPopupMenuWithSignal(obj, toolkit, menulist, ...)
})


## add3rdmousepopupmenu
setMethod("add3rdmousepopupmenu",signature(obj="gWidgetRGtk"),
          function(obj, menulist, action=NULL, ...) {
            .add3rdmousepopupmenu(obj, obj@toolkit,menulist, action, ...)
          })

setMethod("add3rdmousepopupmenu",signature(obj="RGtkObject"),
          function(obj, menulist, action=NULL,...) {
            .add3rdmousepopupmenu(obj, guiToolkit("RGtk2"),menulist, action,...)
          })

setMethod(".add3rdmousepopupmenu",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gWidgetRGtk"),
          function(obj, toolkit, menulist,action=NULL, ...) {
            add3rdMousePopupMenuWithSignal(obj, toolkit,
                                           menulist, action, ...)
          })
setMethod(".add3rdmousepopupmenu",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="RGtkObject"),
          function(obj, toolkit, menulist, action=NULL, ...) {
            add3rdMousePopupMenuWithSignal(obj, toolkit,
                                           menulist, action, ...)
          })


## "dotmethods" defined in dnd.R
## adddropsource
setMethod("adddropsource",signature(obj="gWidgetRGtk"),
          function(obj, targetType="text", handler=NULL, action=NULL, ...) {
            .adddropsource(obj, obj@toolkit,targetType=targetType,
                           handler=handler, action=action, ...)
          })
setMethod("adddropsource",signature(obj="RGtkObject"),
          function(obj, targetType="text", handler=NULL, action=NULL, ...) {
            .adddropsource(obj, guiToolkit("RGtk2"),targetType=targetType,
                           handler=handler, action=action, ...)
          })

## adddropmotion
setMethod("adddropmotion",signature(obj="gWidgetRGtk"),
          function(obj,  handler=NULL, action=NULL, ...) {
            .adddropmotion(obj, obj@toolkit,
                           handler=handler, action=action, ...)
          })
setMethod("adddropmotion",signature(obj="RGtkObject"),
          function(obj, handler=NULL, action=NULL, ...) {
            .adddropmotion(obj, guiToolkit("RGtk2"),
                           handler=handler, action=action, ...)
          })

## adddroptarget
setMethod("adddroptarget",signature(obj="gWidgetRGtk"),
          function(obj, targetType="text", handler=NULL, action=NULL, ...) {
            .adddroptarget(obj, obj@toolkit,targetType=targetType,
                           handler=handler, action=action, ...)
          })

setMethod("adddroptarget",signature(obj="RGtkObject"),
          function(obj, targetType="text", handler=NULL, action=NULL, ...) {
            .adddroptarget(obj, guiToolkit("RGtk2"),targetType=targetType,
                           handler=handler, action=action, ...)
          })


## R Methods
setMethod("dim", "gWidgetRGtk", function(x) .dim(x,x@toolkit))
setMethod("length", "gWidgetRGtk", function(x) .length(x,x@toolkit))
setMethod(".length",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gWidgetRGtk"),
          function(x,toolkit) {
            return(NA)
            gwCat(sprintf("Define length for x of class: %s\n", class(x)[1]))
})
          
setMethod("dimnames", "gWidgetRGtk", function(x) .dimnames(x,x@toolkit))
setReplaceMethod("dimnames",
                 signature(x="gWidgetRGtk"),
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


setMethod("names", "gWidgetRGtk", function(x) .names(x,x@toolkit))
setReplaceMethod("names",
                 signature(x="gWidgetRGtk"),
                 function(x,value) {
                   .names(x,x@toolkit) <- value
                   return(x)
                 })



#### This may be useful to integrate gWidgets with glade
####
## S3 class for coercing to gWidget
as.gWidgetsRGtk2 <- function(widget,...) UseMethod("as.gWidgetsRGtk2")
as.gWidgetsRGtk2.default <- function(widget,...) {
  print(sprintf("No coercion to gWidget available for object of class %s",class(widget)))
  return(widget)
}

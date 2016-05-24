##################################################
### Code to handle toolkits
### A new toolkit should have:
### * a name gWidgetsXXXXX
### * a subclass of guiWidgetsToolkit
### * methods for .functions (.glabel, .svalue, ...)
### * optionally, methods for glabel(obj, ...) dispatching to .glabel(obj,toolkit, ...)

### toolkit code comes *before* others




##################################################
### Components



 
## gcheckbox
## the constructor

## gradio
## constructor


## gdroplist
## the class
## the constructor


## gcheckboxgroup
## the constructor

## gspinbutton
## the constructor

## gslider
## the constructor

## gedit
## the constructor

## gtext
## the constructor

## gaction

## gmenu
## the constructor

## gtoolbar
## the constructor

## gtable

## the constructor

## gdf
## the constructor

## gdfnotebook
## the constructor

## gtree
## the constructor


## gfile (gfilebrowse in base)
## the constructor

## gcalendar
## the constructor

## ggraphics

## the constructor


## gplotnotebook
## the constructor

## gimage

## the constructor

## gsvg

## the constructor


## gstatusbar

## the constructor

## ghtml function


## gseparator
## the constructor

##################################################
## these should be in gWidgets, based on others no toolkit specific code

## the constructor

## ghelp, ghelpbrowser

## the constructor

## ggenericwidget

## the constructor

## gformlayout

## the constructor

## gvarbrowser

## the constructor

##################################################
### Containers

## define a gWidget constructor
## not a generice          

## ggroup
 
## the constructor

## gframe
## the constructor

## gexpandgroup
## the constructor



## gnotebook
## the constructor

## glayout

## the constructor

## gpanedgroup
## the constructor

## icons
##################################################
##
## Methods


## add






## focus

## focus<-


## tooltip


## font

## font<-

## undo

## redo

## tag

## tag<-


## id

## id<-


## toOlkitprovidesgwidgetsdg
## Does thie toolkit provide the widget
## the constructor


##################################################


## ####
## put into RGtk2 only
## ## S3 class for coercing to gWidget
## as.gWidget <- function(obj,...) UseMethod("as.gWidget")
## as.gWidget.default <- function(obj,...) {
##   print(sprintf("No coercion to gWidget available for object of class %s",class(obj)))
##   return(obj)
## }

##' @include misc.R
NULL

##' Show info on Gtk widgets
##'
##' @return NULL
##' @export
showGtkWidgetInfo <- function() {
  if(!faux_require("RGtk2"))
    stop("This function require RGtk2")
  
  ## Main widgets
  ## entry widget
co <- gtkEntry()                        # object

nb <- gtkNotebook()                     # notebook
buffers <- list(                        # text buffers
                properties = gtkTextBuffer(),
                signals = gtkTextBuffer(),
                methods = gtkTextBuffer(),
                methodLookup = gtkTextBuffer()
                )
bnames <- gettext(c("Widget Properties","Widget Signals","Widget methods", "Lookup method"))
names(bnames) <- names(buffers)

## Main GUI
w <- gtkWindow(show=FALSE)
w$setTitle("Browse Gtk Object")
w$setDefaultSize(500,450)               # set size
g <- gtkVBoxNew(spacing=5, homogeneous=FALSE); w$add(g)
statusbar <- gtkStatusbar()
g$PackEnd(statusbar, expand=FALSE, fill=FALSE)
sg <- gtkHBoxNew(); g$packStart(sg, expand=FALSE, fill=FALSE)

sg$packStart(gtkHBox(), TRUE, TRUE,0) # push to right
sg$packStart(gtkLabel("gtk object:"),FALSE,FALSE)
sg$packStart(co, FALSE,FALSE,5)

## populate notebook
g$packStart(nb, TRUE,TRUE)
sapply(c("properties","signals","methods"),
       function(i) {
         tv <- gtkTextView(buffers[[i]])
         sw <- gtkScrolledWindow()
         sw$setPolicy("GTK_POLICY_AUTOMATIC","GTK_POLICY_AUTOMATIC")
         sw$add(tv)
         nb$appendPage(sw, gtkLabel(bnames[i]))
       })
nb$setCurrentPage(0)                    # start on first page

## method lookup is different
mlg <- gtkVBox()
nb$AppendPage(mlg, gtkLabel("method lookup"))

vg <- gtkHBox()
mlg$PackStart(vg, expand=FALSE, fill=FALSE)
l <- gtkLabel("Method name:"); l['xalign'] <- 1
vg$PackStart(l)
methodEntry <- gtkEntry()
meCompletion <- gtkEntryCompletionNew()
meCompletion$setModel(gtkListStore("gchararray"))
meCompletion$setTextColumn(0)
methodEntry$setCompletion(meCompletion)


vg$PackStart(methodEntry)
tv <- gtkTextView(buffers[["methodLookup"]])
sw <- gtkScrolledWindow()
sw$setPolicy("GTK_POLICY_AUTOMATIC","GTK_POLICY_AUTOMATIC")
sw$add(tv)
mlg$PackStart(sw, expand=TRUE, fill=TRUE)

## functions to get information
toLowerFirst <- function(x) {
  tmp <- unlist(strsplit(x,""))
  paste(tolower(tmp[1]),paste(tmp[-1],collapse=""),sep="")
}

lookupObject <- function(objectName, where=.GlobalEnv) {
  object <- get(objectName, envir=where)
  if(inherits(object, "try-error")) {
    statusbar$Push(1, paste("Can't find Gtk object:", objectName, sep=" "))
  } else {
    statusbar$Push(1,"")
  }
  return(object)
}

getPropertyText <- function(obj,
                            constructor=toLowerFirst(class(obj)[1])
                            ) {

  classes <- class(obj)     # vector
  properties <- obj$getPropInfo() #list

  out <- c(constructor,"","","Class inheritance:")
  for(i in 1:length(classes)) 
    out[length(out)+1] <-
      paste(paste(rep("\t",i-1),collapse=""), classes[i],sep="")

  out[length(out)+1] <- "Properties:"
  ## now turn properties into text
  for(i in names(properties)) {
    out[length(out)+1] <- i
    for(j in names(properties[[i]])) {
      out[length(out)+1] <- paste("\t",j,": ",
                                  as.character(properties[[i]][[j]]),
                                  sep="")
    }
  }

  return(out)
}


getSignalText <- function(obj,
                          constructor=toLowerFirst(class(obj)[1])
                          ) {
  signals <- gtkTypeGetSignals(constructor)
  out <- c()
  for(i in names(signals))
    out[length(out)+1] <- paste(i,":\t signal-id ",as.character(signals[[i]]),sep="")
  
  return(out)
}


getMethodsText <- function(obj,
                           constructor=toLowerFirst(class(obj)[1])
                           ) {
  out <- c(paste("Methods for", constructor))
  meths <- apropos(paste("^",toLowerFirst(constructor),sep=""),ignore.case=FALSE)
  out <- c(out,paste("\t",meths,sep=""))
  
  return(out)
}

possibleMethods <- function(obj) {
  ## return vector of possible method names
  l <- c()
  for(i in class(obj)) {
    meths <- apropos(paste("^",toLowerFirst(i), sep=""), ignore.case=FALSE)
    meths <- gsub(toLowerFirst(i), "", meths)
    meths <- meths[meths != ""]
    meths <- meths[!grepl("New",meths)]
    l <- c(l, meths)
  }
  return(l)
}


getMethodForObject <- function(obj, name, where=parent.frame()) {
  ## this is from RGtk2:.getAutoMethodByName
  classes <- c(attr(obj, "interfaces"), class(obj))
  sym <- paste(tolower(substring(classes, 1, 1)), # camelCase
               substring(classes, 2), 
               toupper(substring(name, 1, 1)), 
               substring(name, 2), 
               sep = "")
  which <- sapply(sym, exists, where)
  if (!any(which)) {
    msg <- paste("No such method", name, "for classes", 
                 paste(class(obj),  collapse = ",\n\t"))
  } else {
    methodName <- sym[which][1]
    msg <- paste(paste("Opening help page for", methodName, "\n"))
    ##
    print(help(methodName, help_type="html"))
  }
  return(msg)
}  



showGtkInfo <- function(obj,
                        constructor=toLowerFirst(class(obj)[1])) {

  if(missing(obj)) {
    constructor <- gsub("New$","",constructor)
    FUN <- match.fun(constructor)
    obj <- do.call(constructor, list(show=FALSE))
    if(inherits(obj,"try-error")) {
      statusbar$Push(paste("Can't get object for", constructor, sep=" "),0)
      return(FALSE)
    } else {
      statusbar$Push(1, "")
    }
  }


  ## update text buffers. SetText overwrites current text
  propText <- paste(getPropertyText(obj, constructor),collapse="\n")
  buffers[['properties']]$setText(propText)
  
  signalText <- paste(getSignalText(obj, constructor),collapse="\n")
  buffers[['signals']]$setText(signalText)

  methsText <- paste(getMethodsText(obj, constructor),collapse="\n")
  buffers[['methods']]$setText(methsText)

  ## update completion list for methods
  l <- possibleMethods(obj)
  store <- meCompletion$GetModel()
  store$Clear()
  for(i in l) {
    iter <- store$append()
    store$setValue(iter$iter, 0, i)
  }
}

newObject = function(widget, ...) {
  objectName <- widget$getText()
  object <- lookupObject(objectName)
  if(inherits(object, "try-error")) {
    showGtkInfo(constructor=objectName)
  } else {
    showGtkInfo(obj=object)
  }
  return(FALSE)
}

ID <- gSignalConnect(co,"activate", f=newObject)
ID <- gSignalConnect(co, "focus-out-event", f=newObject)

  
ID <- gSignalConnect(methodEntry, "activate", f=function(widget, user.data) {
  methodName <- widget$GetText()
  object <- lookupObject(user.data$GetText())
  if(inherits(object, "try-error"))
    return()
  txt <- getMethodForObject(object, methodName)
  txt <- paste(txt, collapse="\n")
  buffers[['methodLookup']]$setText(txt)
}, data=co)


## show GUI
w$showAll()

}

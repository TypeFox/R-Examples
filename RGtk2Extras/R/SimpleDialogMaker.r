# Written by Tom Taverner <t.taverner@gmail.com>
# for the U.S. Department of Energy (PNNL, Richland, WA, USA)
# Website: http://omics.pnl.gov/software
#
# Notice: This computer software was prepared by Battelle Memorial Institute,
# hereinafter the Contractor, under Contract No. DE-AC05-76RL0 1830 with the
# Department of Energy (DOE).  All rights in the computer software are reserved
# by DOE on behalf of the United States Government and the Contractor as
# provided in the Contract.
#
# NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY, EXPRESS OR
# IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
#
# This notice including this sentence must appear on any copies of this computer
# software.

#__Generic dialog maker__
#
#This set of functions is intended to assist rapid development of simple dialogs to
#serve as front ends to existing R functions. It is a riff off (or rip off) of John 
#Verzani's traitr toolkit where flexibility is sacrificed for ease of use.
#
#To make a dialog, all that is required is to list the dialog items for each 
#argument of the function. Then, the dialog is set for the function.
#
#Calling run.dialog() on the function will then open the dialog. When "OK" is 
# pressed, the dialog will call the function with the specified arguments.

# Changed eval(parse to "safe.eval" to avoid malicious insertions of code

if(.Platform$OS.type == "windows") {
  PLATFORM_OS_TYPE <- "windows"
} else if (.Platform$OS.type == "unix"){
  if (Sys.info()["sysname"] == "Darwin")
    PLATFORM_OS_TYPE <- "mac"
  else{ 
    PLATFORM_OS_TYPE <- "unix"    }
}

library(RGtk2)
#library(gWidgetsRGtk2)

my.gvarbrowser <- function(...){
    if(require(gWidgetsRGtk2))
       return(gvarbrowser(...))
    warning("gWidgetsRGtk2 not installed")
}

my.getToolkitWidget <- function(...){
    if(require(gWidgetsRGtk2))
       return(getToolkitWidget(...))
    warning("gWidgetsRGtk2 not installed")
}

my.svalue <- function(...){
    if(require(gWidgetsRGtk2))
       return(svalue(...))
    warning("gWidgetsRGtk2 not installed")
}

my.getwd <- function() gsub("(.*)/$", "\\1", getwd()) # for "C:/" we want to remove the end "/"

.ITEMS = c("choiceItem", "rangeItem", "fileItem", "integerItem", "radiobuttonItem", "objectItem",
      "numericItem", "stringItem", "trueFalseItem", "variableSelectorItem", "listItem", "dataframeItem", "variableStringItem", "buttonItem")
.DATA_OBJECTS <- c("data.frame", "matrix")


# we use our own version of gfile, my.gfile

# Utility function for choosing file with filters in platform independent way
my_choose_files <- function(fn ="*.*", multi=FALSE, filters=NULL, type = "open", caption=NULL){
  captions <- list(open = "Select File", save= "Save File", selectdir = "Select Directory")
  stopifnot(type%in%names(captions))
  if(is.null(caption)) caption <- captions[[type]]
  if(identical(.Platform$OS.type, "windows")){
    if(identical(type, "selectdir")) {
      f <- choose.dir(fn)
    } else {
      f <- choose.files(fn, multi = multi, filters=filters, caption=caption)
    }
  } else {
      f <- my.gfile(fn, multi = multi, type = type)
  }
  return(f)
}    

  # This function binds its arguments together into a ragged column-wise array
ragged.cbind <- function(...){
  items <- list(...)
  if(length(items) == 1) items <- items[[1]]
  for(item in items)
    stopifnot(all(sapply(item, is.atomic)))
  dd <- function(x){
    if(!is.null(dim(x))) return(dim(x))
    return(c(length(x), 1))
  }
  num.rows <- max(sapply(items, function(x) dd(x)[1]))
  rv <- NULL
  for(ii in 1:length(items)){
    item <- items[[ii]]
    tmp <- array(NA, c(num.rows, dd(item)[2]))
    if(!is.null(dim(item))){ 
      for(jj in 1:dd(item)[2]) tmp[1:dd(item)[1],jj] <- item[,jj]
      colnames(tmp) <- colnames(item)
    } else {
      if(dd(item)[1] > 0)
         tmp[1:dd(item)[1],] <- unlist(item)
      colnames(tmp) <- names(items)[[ii]]
    }                                
    rv <- cbind(rv, tmp)
  }
  return(rv)
}

  # exists("a$b$c")
object.exists <- function(text, envir=parent.env(environment())){
  if(!nchar(text) || !is.character(text) || length(text) != 1) return(FALSE)
  xx <- strsplit(text, "$", fixed=T)[[1]]
  if(!nchar(xx) || !length(xx) || !exists(xx[1])) return(FALSE)
  x1 <- get(xx[1], envir)
  if(length(xx) > 1) for (item in xx[-1]) {
    if(!item%in%names(x1)) return(FALSE)
    x1 <- x1[[item]]
  }
  return(TRUE)
}
  # Safely evaluates "a$b$c"
safe.eval <- function(text, envir=parent.env(environment())){
  xx <- strsplit(text, "$", fixed=T)[[1]]
  tryCatch({
    x1 <- get(xx[1], envir)
    if(length(xx) > 1) for (item in xx[-1]) x1 <- x1[[item]]
  }, error=  function(e) stop(paste("error in safe.eval:", as.character(e))))
  return(x1)
}

get.value <- function(dlg.item, ...) {
  #if(is.function(dlg.item$getData("get.value")))
  dlg.item$getData("get.value")(dlg.item, ...)
}

get.name <- function(dlg.item) dlg.item$getData("name")

set.value <- function(dlg.item, ...){
  if(is.function(dlg.item$getData("set.value")))
   dlg.item$getData("set.value")(dlg.item, ...)
}

signal.connect <- function(obj, f, data=NULL, signal="default"){
   signals <- obj$getData("signals")
   if(is.null(signals)) signals <- list()
   if(is.null(signals[[signal]])) signals[[signal]] <- list()
   signals[[signal]][[length(signals[[signal]])+1]] <- list(f=f, data=data)
   obj$setData("signals", signals)
}

get.signal <- function(dlg.item, signal="default") 
   dlg.item$getData("signals")[[signal]]

signal <- function(obj, signal="default"){  
   sapply(get.signal(obj, signal), function(x) {
    do.call(x$f, as.list(c(obj, x$data)))
   })
}

# Puts a frame around the item, if the item has attribute "suppress frame" it doesn't
myFrame <- function(item, label){
  name <- item$getData("name")
  suppress.frame <- !is.null(item$getData("suppress.frame"))
  if(suppress.frame){
    fram <- gtkHBoxNew() 
    fram$packStart(item, TRUE, TRUE, 5)
  } else {
    fram <- gtkFrameNew(label=label)
    fram$setBorderWidth(2)
    vbox <- gtkVBoxNew()
    hbox <- gtkHBoxNew()
    hbox$packStart(item, TRUE, TRUE, 5)
    vbox$packStart(hbox, TRUE, TRUE, 5)
    fram$add(vbox)      
  }  
  #fram$setData("item", item)
  #fram$setData("name", name)
  return(fram)
}

stringItem <- function(value, name="", multi=F, ...){
  if(!multi){
    item <- gtkEntryNew()
    item$setText(value)
    item$setData("name", name)
    item$setData("get.value", gtkEntryGetText)
    item$setData("set.value", function(item, value, ...) item$setText(value))
    gSignalConnect(item, "activate", function(...) {
      signal(item)
    })    
  } else {    
    item <- gtkVBoxNew()    
    sw <- gtkScrolledWindowNew()
    txt <- gtkTextViewNew()
    sw$add(txt)
   	sw$setPolicy("automatic", "automatic")
    buf <- txt$getBuffer()
    item$add(sw)
    item$setData("name", name)    
    item$setData("get.value", function(item)
       buf$getText(buf$getStartIter()$iter, buf$getEndIter()$iter))
    item$setData("set.value", function(item, value, ...) buf$setText(value))
    item$setData("expand", TRUE)  # expand in column    
    gSignalConnect(txt, "key-press-event", function(obj, evt) {
      if(evt[["keyval"]] == GDK_Return && (as.flag(evt[["state"]] & GdkModifierType['control-mask'])))
         signal(item)
      return(FALSE)
    }) 
  }
  return(item)
}

# This creates a widget that allows selection of multiple ranges of values.
# By Iago Conde
variableStringItem  <- function(value, name="", show.arrows, ...){
 hbox0 <- gtkHBoxNew(FALSE, 0)     
 
 button.add <- gtkButtonNewWithLabel(">>")

 gSignalConnect(button.add, "clicked", function(...) signal(hbox0, "add"))

 hbox0$packStart(button.add, FALSE, FALSE, 0)
 ent <- gtkEntryNew()
 ent$setText(value)
 hbox0$packEnd(ent,TRUE,TRUE,0)
 hbox0$setData("ent", ent)

 hbox0$setData("get.value", function(hbox0,selected=FALSE){
    texto <- hbox0$getData("ent")$getText()
    return(texto)
  }) 

  hbox0$setData("name",name)
  hbox0$setData("set.value", function(item, value, ...){  
    ent <- item$getData("ent")$getText()
    new_string <- paste(value, collapse = "; ")
    if(nchar(ent)) new_string <- paste(ent, new_string)
    item$getData("ent")$setText(new_string)    
  })

  return(hbox0)
}
    
choiceItem <- function(value=NULL, values=NULL, item.labels=NULL, name="", ...){

  if(is.null(values)) values <- value
  if("value"%in%names(value))
    default.value <- value[["value"]]
  else 
    default.value <- value[1]

  value <- default.value

  stopifnot(value%in%values)  
  if(is.null(item.labels)) item.labels <- values
  stopifnot(length(item.labels) == length(values))
  
  item <- gtkComboBoxNewText()
  MAX_WIDTH = 150
  sr <- item$getSizeRequest()
  item$setSizeRequest(MAX_WIDTH, -1)      
  for (ll in item.labels) item$appendText(ll)  
  item$setData("values", values)
  item$setData("labels", item.labels)

  if(length(values)) item$setActive(which(values %in% value)-1)
  changed.id <- gSignalConnect(item, "changed", function(...) {
    signal(item)
  }) 
    # Emit the signal if the box is popped up, to refresh
  #gSignalConnect(item, "popup", function(...) {
  #  signal(item)
  #})   

  item$setData("changed.id", changed.id)     

  item$setData("name", name)  
  
  item$setData("get.value", function(x, selected=TRUE, ...){
    item.labels <- item$getData("labels")
    values <- item$getData("values")
    if(selected){
      stopifnot(length(item.labels) == length(values))
      idx <- which(item.labels%in%x$getActiveText())
      return(values[idx])
    } else {
      return(values)
    }
  })

  item$setData("set.value", function(item, values, selected=NULL, ...){
    old.values <- item$getData("values")
    gSignalHandlerBlock(item, item$getData("changed.id"))
    if(length(old.values) > 0)
      for (ii in 1:length(old.values))
        item$removeText(0)
	  item$setData("values", values)
    item$setData("labels", values) # remove any labels
    for (v in values) {
      item$appendText(v)
      sr <- item$getSizeRequest()
      if(sr$width > MAX_WIDTH) item$setSizeRequest(MAX_WIDTH, -1)      
    }
	  if(is.null(selected) && length(values)) {
      item$setActive(0)
    } else if (!is.null(selected) && length(selected)==1 && selected%in%values){
      item$setActive(which(values%in%selected)-1)
    }
    gSignalHandlerUnblock(item, item$getData("changed.id"))        
  })

  return(item)
} 

buttonItem <- function(value, name = "", ...){
  item <- gtkHBoxNew()
  button <- gtkButtonNewWithLabel(value)
  item$packEnd(button, FALSE, FALSE, 0)

  item$setData("name", name)
  item$setData("value", NULL)

  gSignalConnect(button, "clicked", function(...){
    signal(item, "clicked")
  })

  item$setData("get.value", function(x) {
    item$getData("value")
  })

  item$setData("set.value", function(item, value, ...){
    item$setData("value", value)
  })

  item$setData("suppress.frame", TRUE)  

  return(item)
}


  # Runs a dialog from a buttonItem given the list
dialog_button <- function(item, user.data){
  rv <- run.dialog(list, dlg.list=user.data)$retval
  if(!is.null(rv)) set.value(item, rv)
}


radiobuttonItem <- function(value=NULL, values=NULL, 
    item.labels=NULL, name="", by.index=NULL, ...){
    
  if("value"%in%names(value))
    default.value <- value[["value"]]
  else 
    default.value <- value[1]

  item <- gtkVBoxNew()
  if(is.null(values)) values <- value
  if(is.null(by.index)) by.index <- FALSE            
  value <- default.value  
  stopifnot(value%in%values)
  if(is.null(item.labels)) 
    item.labels <- values  
  stopifnot(length(item.labels) == length(values))
  item$setData("values", values)
  
  stopifnot(length(values) > 1)
  the.choice <- which(values%in%value)
  grb <- gtkRadioButtonNewWithLabel(NULL, label=item.labels[1])
  if(the.choice == 1) grb$setActive(TRUE)
  item$add(grb)
  for(jj in 2:length(values)){
    grb <- gtkRadioButtonNewWithLabel(group=grb$getGroup(),label=item.labels[jj])
    if(jj == the.choice) grb$setActive(TRUE)          
    item$add(grb)
  }
  sapply(grb$getGroup(), function(x) gSignalConnect(x, "toggled", after=T, function(...) {
    if(x$getActive()) signal(item)    
    }))

  item$setData("name", name)
  item$setData("grb", grb)  
  item$setData("get.value", function(x){
    idx <- which(rev(sapply(x$getData("grb")$getGroup(), gtkToggleButtonGetActive)))
    return(values[idx])
  })
  
  item$setData("set.value", function(x, value, ...){  
    values <- x$getData("values")
    stopifnot(value%in%values || length(value) == 1)
    grb <- rev(x$getData("grb")$getGroup())
    idx <- which(values==value)
    grb[[idx]]$setActive(TRUE)
  })
  return(item)        
}


trueFalseItem <- function(value, name="", label = NULL, ...){
  stopifnot(is.logical(value))  
  if(is.null(label)) label <- name  
  item <- gtkCheckButtonNewWithLabel(label=label)
  item$setActive(value)
  gSignalConnect(item, "toggled", after=T, function(item) {
    signal(item, "default")
  })
  item$setData("name", name)
  item$setData("suppress.frame", TRUE)  
  item$setData("get.value", gtkToggleButtonGetActive)
  item$setData("set.value", function(item,value, ...){
    item$setActive(value)
    #signal(item, "default")
  })

  return(item)
}

numericItem <- function(value, name="", ...){
  stopifnot(!is.na(as.numeric(value)))
  item <- gtkEntryNew()
  item$setText(value)
  item$setData("name", name)
  item$setData("get.value", function(x) {
      as.numeric(x$getText())
    })
  item$setData("set.value", function(item, value, propagate=TRUE, ...) {
    item$setText(value)
    if(propagate) signal(item, "default")
  })
  return(item)
}

integerItem <- function(value, name="", from=-.Machine$integer.max, to=.Machine$integer.max, ...){
  stopifnot(!is.na(as.integer(value)))
  
  if("from"%in%names(value)) from <- value["from"]  
  if("to"%in%names(value)) to <- value["to"] 
  if("value"%in%names(value)) value <- value["value"]   
  
  item_adj <- gtkAdjustment(value=value, lower=from, upper=to, step.incr=1, page.incr=10)
  item <- gtkSpinButton(item_adj, 1.0, 0)
  item$setValue(value)
    
  item$setData("name", name)
  item$setData("get.value", function(x) {
     x$getValue()
    })
  item$setData("set.value", function(item, value, propagate=TRUE, ...) {
    item$setValue(value)
    if(propagate) signal(item, "default")
  })
  return(item)
}

rangeItem <- function(value, name="", from=0, to=100, by=1, ...){
  stopifnot(!is.na(as.numeric(value)))
  if("from"%in%names(value)) from <- value["from"]  
  if("to"%in%names(value)) to <- value["to"]  
  if("by"%in%names(value)) by <- value["by"]    
  if("value"%in%names(value)) value <- value["value"]      
  item <- gtkHScaleNewWithRange(min=from, max=to, step=by)
  if(length(value) > 1) value <- value[1]  
  item$setValue(value)
  gSignalConnect(item, "value-changed", function(obj){
    signal(item, "default")
    return(value)
  })
    
  item$setData("name", name)
  item$setData("get.value", function(x) {
     x$getValue()
    })
  item$setData("set.value", function(item, value, propagate=TRUE, ...){
    item$setValue(value)
    if(propagate) signal(item, "default")
  })
  return(item)
}

fileItem <- function(value, name="", multi = FALSE, filters=NULL, extension = "*", type = "open", ...){

  if(identical(.Platform$OS.type, "windows") && is.null(filters)) filters <- get("Filters", envir=.GlobalEnv)
  
  stopifnot(type%in%c("open", "save", "selectdir"))
  item <- gtkHBoxNew()
  #label <- gtkLabelNew(value)
  #label$setEllipsize(PangoEllipsizeMode['start'])
  label <- gtkEntryNew()
  label$setText(value)
  label$setAlignment(0)  
  button <- gtkButtonNewWithLabel(paste(c("Select", ifelse(identical(type,"selectdir"), " Folder", " File"), ifelse(multi, "(s)", "")), collapse=""))
  item$add(label)
  item$packEnd(button, FALSE, FALSE, 0)
  item$setData("label", label)
  gSignalConnect(button, "clicked", data=label, function(obj, label){
    #tryCatch({
     #fn <- file.path(my.getwd(), paste("*", extension, sep=""))
     fn <- extension
     if(identical(.Platform$OS.type, "windows")){
       if(identical(type, "selectdir")) 
         f <- choose.dir(fn)
       else
         f <- choose.files(fn, multi = multi, filters=filters)
     } else {
       f <- my.gfile(fn, multi = multi, type = type)
     }
     if(!is.na(f) && length(f) > 0 && nchar(f) > 0){
        set.value(item, f)
        #tryCatch({
        setwd(dirname(f))
        #}, error = function(e) setwd(dirname(f)))
      }         
    #}, error = function(e){print(e)})
  })
    
  item$setData("name", name)

  item$setData("get.value", function(x){
     val <-  x$getData("label")$getText()
       # Check we're evaluating either a string or a call to c(...)
     rv <- tryCatch({
       pt <- parse(text=val)       
       if( length(pt) > 0 && (identical(call("c"), pt[[1]][1]) ||  identical(class(pt[[1]]), "character"))){
         eval(pt)
       } else {
         val
       }
     }, error = function(e) val)
     return(rv)
   })

  item$setData("set.value", function(x, value, propagate=TRUE, ...) {
    x$getData("label")$setText(deparse(value))
    if(propagate) signal(item, "default")
  })

  return(item)
}

objectItem <- function(value=NULL, name="", tooltip="", data.types=NULL, parent.window=NULL, ...){
  item <- gtkHBoxNew()
  label <- gtkEntryNew()
  if(!length(value)){
    label$setText("")
  } else {
    label$setText(value)
  }  
  #label$setText(deparse(substitute(value)))
    
  label$setAlignment(0)    
  #label$setSizeRequest(150, -1)
  button <- gtkButtonNewWithLabel("Select Data...")  
  item$add(label)
  bb1 <- gtkHBoxNew()
  item$packEnd(button, FALSE, FALSE, 0)
  item$setData("label", label)
  
  gSignalConnect(button, "clicked", data=label, function(obj, data){
    ow = options("warn")[[1]]
    options(warn = -1)  # suppress warnings from vb
    label=data      
      # Changed to "select" from "gtk-close"
    dialog <- gtkDialog("Select...", parent.window, "modal", "Select", 1,show = T)
    dialog$setSizeRequest(400, 300)
    
    dialog$setPosition(GtkWindowPosition["center-on-parent"])              
    dialog$setTransientFor(parent.window)          
    dialog$showAll()

    # library(gWidgetsRGtk2)
    tryCatch({   # this warns about missing .tags on windows
print("Got here")
      vb <- my.gvarbrowser(action = list(item=item, label=label), 
        handler= function(obj, tv, path, column, ...){
          action <- obj$action
          label <- action$label
          item <- action$item 
          if(identical(data.types, "function"))
             my.getToolkitWidget((vb@widget)@filter)$setActive(4)
             
          if(my.getToolkitWidget(vb)$getSelection()$countSelectedRows() > 0){
            choice <- my.svalue(vb)
            get.choice <- safe.eval(choice, envir=.GlobalEnv)
            if(is.null(data.types) || class(get.choice)%in%data.types){          
              set.value(item, choice)
              dialog$destroy()        
            }
          }
          options(warn = ow)
        })
     }, error= function(e) {}, warning= function(w) print(w))

     dialog[["vbox"]]$add(my.getToolkitWidget(vb)$getParent()$getParent()$getParent())

     if( dialog$run() == 1){
       if(my.getToolkitWidget(vb)$getSelection()$countSelectedRows() > 0){     
         choice <- my.svalue(vb)
         get.choice <- safe.eval(choice, envir=.GlobalEnv)
         if(is.null(data.types) || class(get.choice)%in%data.types) set.value(item, choice)
       }
       options(warn = ow)
       dialog$destroy()       
     }
  })

  item$setData("button", button)   
  item$setData("name", name)

  item$setData("get.value", function(x){
    #txt = x$getData("label")$getText()
    #return(as.call(parse(text=txt))[[1]])  
    return(x$getData("label")$getText())
  })
  item$setData("set.value", function(x, value, propagate=TRUE, ...) {
    #x$getData("label")$setText(deparse(substitute(value)))
    label <- x$getData("label")
    if(!length(value)){
      label$setText("")
    } else {
      label$setText(value)
    }
    if(propagate) signal(item, "default")
  })    

  gSignalConnect(label, "activate", function(...) {
    signal(item)
  }) 
  
  #set.value(item, value)
  return(item)
}

# An item to allow any data frame to be selected
dataframeItem <- function(data.types = .DATA_OBJECTS, ...)
  objectItem(data.types=data.types, ...)

quick_message <- function(message="", win=NULL, caption="Warning") {
  dialog <- gtkDialog(caption, NULL, "destroy-with-parent", "gtk-ok", 1, show = FALSE)
  dialog$setPosition(GtkWindowPosition["center"])                      
  label <- gtkLabel(message)     
  label$setLineWrap(TRUE)
  gSignalConnect(dialog, "response", gtkWidgetDestroy)  
  hbox <- gtkHBoxNew()
  hbox$packStart(label, TRUE, TRUE, 20)
  dialog[["vbox"]]$packStart(hbox, TRUE, TRUE, 10)
  dialog$setResizable(FALSE)  
  dialog$showAll()
  dialog$setKeepAbove(TRUE)
  dialog$run()
}

# source("/home/larman/research/RGtk/inst/demo/treeStore.R")
# From Rgtk demo

# Choice item
create.model <- function(row.items, col.items, initial=TRUE){

  theFrame <- data.frame("rows" = row.items, array(initial, c(length(row.items), length(col.items))))
  colnames(theFrame)[-1] <- col.items  
  model <- rGtkDataFrame(theFrame)
  return(model)
}

item.toggled <- function(cell, path.str, data){
  treeview <- data$treeview
  vbox <- data$vbox
  checkPtrType(treeview, "GtkTreeView")
  model <- gtkTreeViewGetModel(treeview)
#  theFrame <- as.data.frame(model)
  
  path <- gtkTreePathNewFromString(path.str)
  column <- cell$getData("column")
#  theFrame[as.numeric(path$toString())+1, column+1] <- !theFrame[as.numeric(path$toString())+1, column+1]
#  treeview$setModel(rGtkDataFrame(theFrame))
  i <- as.numeric(path$toString())+1
  j <- column + 1
  model[i, j] <- !model[i, j]
  return()
}

 add.item <- function(treeview, item) {
     model <- treeview$getModel()
     iter <- model$append()$iter     
     model$set(iter, 0, item)
 }

  list.selected.items <- function(treeview)
  {
     checkPtrType(treeview, "GtkTreeView")
     model <- treeview$getModel()
     selection <- treeview$getSelection()
   
     selected <- selection$getSelectedRows()$retval
     rv <- vector("character", length(selected))
     if (length(selected)) for(ii in 1:length(selected)){
        path <- selected[[ii]]
        iter <- model$getIterFromString(path$toString())$iter
        rv[ii] <- model$get(iter, 0)[[1]]
      }         
     return(rv)
  }

column.clicked <- function(column, data){
	col.idx <- column$getData("column.index")
	#return(FALSE)
	treeview <- data
  columns <- treeview$getColumns()
  ncols <- length(columns)
	model <- treeview$getModel()
	new.state <- column$getData("next.state")
	column$setData("next.state", !new.state)
#	path <- gtkTreePathNewFromString("0")
#	iter <- model$getIter(path)$iter
#  model$foreach(user.data=col.idx, function(model, path, iter, col.idx){
#    model$set(iter, col.idx, new.state)
#    return(FALSE)
#  })	
#  
  theFrame <- as.data.frame(model)
  if(dim(theFrame)[1] > 0) 
    theFrame[, col.idx+1] <- !new.state
  treeview$setModel(rGtkDataFrame(theFrame))  
  
  return(FALSE)
}

add.columns <- function(treeview, col.items, initial=integer(0), vbox) {

  # column for holiday names
  renderer <- gtkCellRendererTextNew()
  renderer['ellipsize'] <- PangoEllipsizeMode['end']
  renderer['ellipsize-set'] <- TRUE  
  renderer$set(xalign = 0.0)

  col.offset <- treeview$insertColumnWithAttributes(-1, "Item", renderer, 
  								text = 0)
								
  column <- treeview$getColumn(col.offset - 1)
  column$setClickable(TRUE)
	gtkTreeViewColumnSetFixedWidth(column, 150)
	gtkTreeViewColumnSetSizing(column, GtkTreeViewColumnSizing['fixed'])
  
  if(length(col.items)) column$setResizable(TRUE)
  
  if(length(col.items)) for(jj in 1:length(col.items)){
		renderer <- gtkCellRendererToggleNew()		
		renderer$set(xalign = 0.0)
		renderer$setData("column", jj)                           
		gSignalConnect(renderer, "toggled", item.toggled, data=list(treeview=treeview, vbox=vbox))
		col.offset <- treeview$insertColumnWithAttributes(-1, col.items[jj], renderer,
									  active = jj)
		column <- treeview$getColumn(col.offset - 1)
		column$setData("column.index", jj)		
		column$setSizing("fixed")
		column$setFixedWidth(50)
    column$setClickable(TRUE)		
		if(jj != length(col.items)) column$setResizable(TRUE)
    column$setData("next.state", !jj%in%initial)
    gSignalConnect(column, "clicked", column.clicked, data=treeview)    
  }
 
}

choice.grid <- function(row.items, col.items=character(0), initial=integer(0), xor=TRUE, 
  headers.visible=T, select.mode="multiple", handler.selected=NULL, user.data=NULL){

  GetItems <- function(treeview){
  	model <- treeview$getModel()
  	ncols <- length(treeview$getColumns())-1
    dim <- treeview$getData("dim")
    retval <- array(F, dim) 
    the.dimnames <- treeview$getData("dimnames")
  
    if(ncols > -1){      
#      .env <- new.env()
#      .env$retval <- retval
#      
      theFrame <- as.data.frame(model)
#
##      model$foreach( user.data=.env, function(model, path, iter, .env){
#     		iter <- model$getIter(path)$iter  
#        ii <- path$getIndices()+1   
#        if(ncols > 0) for(jj in 1:ncols) .env$retval[ii, jj] <- model$get(iter, jj)[[1]]     
#        return(FALSE)
#      })
#      retval <- .env$retval
                     
      if(length(dim(theFrame)) > 1 && dim(theFrame)[2] > 0){
        rn <- theFrame[,1]
        cn <- colnames(theFrame)
        retval <- array(as.logical(as.matrix(theFrame[,-1,drop=F])), dim(theFrame[,-1,drop=F]))
        rownames(retval) <- rn  
        colnames(retval) <- cn[-1]
      } else {
        #print("Frame NOT seen")
      }
      
  	} else {  # just return selected values for 0 columns
  	    retval <- integer(0)
        tryCatch({    
          retval <- gtkTreeSelectionGetSelectedRows(gtkTreeViewGetSelection(treeview))$retval
          if(length(retval)){
            retval <- sapply(retval, gtkTreePathGetIndices)+1
          }        
        },
        error=function(e) {
          warning(e)
        })
        if(!length(retval)) retval <- integer(0)      
        retval <- the.dimnames$rows[retval]
  	}
    retval
  }
  
	# create tree view
	make.treeview <- function(vbox, row.items, col.items, initial=integer(0)){
	
  	sw <- gtkScrolledWindowNew(NULL, NULL)
  	sw$setShadowType("etched-in")
  	sw$setPolicy("automatic", "automatic")
  	vbox$packStart(sw, TRUE, TRUE, 0)	
    vbox$setData("sw", sw)
  	
  	 ## this returns an RGtkDataF
  	model <- create.model(row.items, col.items, initial)
  	treeview <- gtkTreeViewNewWithModel(model)
  	treeview$setRulesHint(TRUE)
  	treeview$setHeadersVisible(headers.visible)
    treeview$getSelection()$setMode(GtkSelectionMode[select.mode])
    vbox$setData("model", model)  	        
    sw$setData("treeview", treeview)  	    
  
#    selectedColor <- as.GdkColor(c(49, 106, 197)*256) # Linux
#    treeview$modifyBase(GtkStateType['selected'],  selectedColor)
#    treeview$modifyBase(GtkStateType['active'], selectedColor)
#    treeview$modifyText(GtkStateType['selected'], as.GdkColor("white"))
#    treeview$modifyText(GtkStateType['active'], as.GdkColor("white"))	
  	
    treeview$setData("xor", xor)
    treeview$setData("dim", c(length(row.items), length(col.items)))
    treeview$setData("dimnames", list(rows=row.items, columns=col.items))
    treeview$setData("GetItems", GetItems)  
    
  	add.columns(treeview, col.items, initial, vbox)		
  
  	sw$add(treeview)
    gSignalConnect(treeview$getSelection(), "changed", 
      data=list(vbox=vbox, treeview=treeview),    
      function(treeselection, data){
        handler.selected <- vbox$getData("handler.selected") 
        user.data=vbox$getData("user.data")                 
        if(is.function(handler.selected)){
          treeview <- data$treeview
          retval <- GetItems(treeview)
          handler.selected(retval, user.data)
        }
        return(FALSE)
      })    

    # watch out for mapping more than once          	
    #gSignalConnect(sw, "map", map.treeview, data=list(treeview=treeview, col.items=col.items))
  	return(treeview)
 	}
 	
 
 	replace.treeview <- function(vbox, row.items, col.items=character(0), initial=TRUE){

    treeview <- vbox$getData("treeview")
    theFrame <- data.frame("rows" = row.items, array(initial, c(length(row.items), length(col.items))))
    colnames(theFrame)[-1] <- col.items    
 	  treeview$setModel(rGtkDataFrame(theFrame))
 	
    vbox$setData("treeview", treeview)	 	  
 	}

	vbox <- gtkVBoxNew(FALSE, 8)
	vbox$setBorderWidth(8)
  vbox$setData("handler.selected", handler.selected)    	
  vbox$setData("user.data", user.data)    	  
  vbox$setData("replace.treeview", replace.treeview)   
	treeview <- make.treeview(vbox, row.items, col.items, initial)	
  vbox$setData("treeview", treeview)   
  return(vbox)
}

# This creates a widget that allows selection of multiple ranges of values.
variableSelectorItem <- function(value, values=integer(0), name="", tooltip="", ...){
  if(!length(values)) values <- "Select All" # our default column name 

  vbox <- choice.grid(character(0), col.items="Select All", initial =TRUE,
    headers.visible=T, select.mode="none")
    
  vbox$setBorderWidth(0)        
  vbox$setSizeRequest(200, 150)  
  vbox$setData("name", name)
  vbox$setData("values", values) 
  vbox$setData("get.value", function(vbox){
    treeview <- vbox$getData("treeview")
    GetItems <- treeview$getData("GetItems")
    stopifnot(is.function(GetItems))
    values <- vbox$getData("values")  
    retval <- GetItems(treeview)
    return(row.names(retval)[retval])
  })

  vbox$setData("set.value", function(vsi, values, initial=TRUE, propagate=TRUE, ...){
		replace.treeview <- vsi$getData("replace.treeview")
		#if(length(dim(values)) < 2){
  	replace.treeview(vsi, values, col.items="Select All", initial)
		#} else if (length(dim(values)) == 2){
		#}
    vsi$setData("value", values)
    if(propagate) signal(vsi, "default")
  })
  
  tvc <- vbox$getData("treeview")$getColumns()
  if( length(tvc) > 1){
    gSignalConnect(tvc[[2]], "clicked", function(...) signal(vbox, "default"))
    renderer <- tvc[[2]]$getCellRenderers()[[1]]		
		gSignalConnect(renderer, "toggled", function(...) signal(vbox, "default"), vbox$getData("treeview"))
  }
  vbox$setData("expand", TRUE)  # expand in column
  set.value(vbox, value)
  return(vbox)  
}

# This just makes a new name in .global guaranteed to be unique.
unique.global.name <- function(nam){
	ll <- ls(envir=.GlobalEnv)  							
	new.name <- make.unique(c(ll, nam))[length(ll)+1]
	return(new.name)
}

# This creates a widget that allows selection of multiple ranges of values.

listItem <- function(value, 
  return.selection = TRUE,
  name="", 
  show.arrows= TRUE, 
  max.items = NULL,
  remove.items = TRUE, 
  select.mode = "multiple", 
  ...){
                                        
  hbox0 <- gtkHBoxNew(FALSE, 0)

  stopifnot(select.mode%in%c("none", "single", "multiple"))
  vbox <- choice.grid(character(0), col.items=character(0), initial =character(0),
    headers.visible=F, select.mode = select.mode)                    
  hbox0$packEnd(vbox, TRUE, TRUE, 0)

  box3 <- gtkVBoxNew(FALSE, 0)
  box2 <- gtkVBoxNew(FALSE, 0)
  hbox0$packStart(box3, FALSE, FALSE, 0)
  if(show.arrows) box3$packStart(box2, FALSE, FALSE, 10)

  #button1 <- gtkButtonNewWithLabel("/\\")
  #box2$packStart(button1, FALSE, FALSE, 0)

  #button2 <- gtkButtonNewWithLabel("\\/")
  #box2$packStart(button2, FALSE, FALSE, 0)

  button.add <- gtkButtonNewWithLabel(">>")
  box2$packStart(button.add, FALSE, FALSE, 0)

  button.remove <- gtkButtonNewWithLabel("<<")
  box2$packStart(button.remove, FALSE, FALSE, 0)
  
  gSignalConnect(button.add, "clicked", function(...) signal(hbox0, "add"))
  gSignalConnect(button.remove, "clicked", function(...) signal(hbox0, "subtract"))

  hbox0$setBorderWidth(0)        
  hbox0$setSizeRequest(200, 150)  
  if(show.arrows) hbox0$setSizeRequest(220, 150)  
  hbox0$setData("name", name)
  hbox0$setData("remove.items", remove.items)
  hbox0$setData("max.items", max.items)

  hbox0$setData("get.value", function(hbox0, selected=FALSE){
    treeview <- hbox0$getData("vbox")$getData("treeview")
    rv <- NULL
    if(selected) {
      rv <- list.selected.items(treeview)
    } else {
      rv <- rownames(treeview$getData("GetItems")(treeview))
    }
    return(rv)
  })  

  hbox0$setData("vbox", vbox)
  hbox0$setData("expand", TRUE)    
  hbox0$setData("set.value", function(item, values=NULL, ...){
    if(is.numeric(max.items) && max.items >= 1 && length(values) > max.items) 
      values <- values[1:max.items] 
	  checkPtrType(item, "GtkBox")
    vbox <- item$getData("vbox")
	  replace.treeview <- vbox$getData("replace.treeview")
	  replace.treeview(item$getData("vbox"), values, col.items=character(0), initial=1)
	  item$setData("value", values)
  }) 

  selection <- vbox$getData("treeview")$getSelection()
  gSignalConnect(selection, "changed", function(...) {
    signal(hbox0, "default")
  })  
  #gSignalConnect(selection, "changed", function(...) {
  #  signal(hbox0, "default")
  #  print(get.value(hbox0))
  #})

  set.value(hbox0, value)    
  return(hbox0)
}

# this function turns a list of markup items into a list that has 
# one element per markup item
process.markup <- function(dlg.list, func=NULL){
                                                                                   
  if(any(!nchar(names(dlg.list)))) 
    stop(paste("Dialog markup contains unnamed elements at", paste(which(!nchar(names(dlg.list))), collapse=", ")))

  items.regex <- paste("[.]", .ITEMS, "$", collapse="|", sep="")
  pos.of.items <- regexpr(items.regex, names(dlg.list)) # where are the Items?
  
  items.idx <- which(pos.of.items > -1)
  possible.items.idx <- which(regexpr(paste("[iI]tem$", collapse="|", sep=""), names(dlg.list)) > -1)
  if(!all(possible.items.idx%in%items.idx)) {
    cat(paste('Found', paste(names(dlg.list)[possible.items.idx[!possible.items.idx%in%items.idx]], collapse=", "), "in dialog list but could not identify as markup item(s).\n"))
    flush.console()
    }
  dlg.item.names <- substr(names(dlg.list), pos.of.items+1, pos.of.items+attr(pos.of.items, "match.length"))[items.idx]      
  retval <- list()
      
  markup <- list()
  if(length(items.idx) > 0 && items.idx[1] > 1)
    markup <- (dlg.list[1:(items.idx[1]-1)])
  retval[[1]] <- markup
  names(retval)[1] <- "main"
  
  arg.names <- substr(names(dlg.list), 1, pos.of.items-1)[items.idx]
  if(any(duplicated(arg.names))) stop("Duplicate arguments in function dialog")  
#    # Warn for potential mismatches between function and dialog arguments
#  if(!is.null(func) && is.function(func)){
#    nff <- names(formals(func))
#    if(length(nff) && !"..."%in%arg.names && any(!arg.names%in%nff))
#       cat(paste("Found argument(s) called:", paste(arg.names[!arg.names%in%nff], collapse = ", "), "in function, but not in dialog.\n"))
#    if(length(arg.names) && !"..."%in%nff && any(!nff%in%arg.names))
#      cat(paste("Found argument(s) called:", paste(nff[!nff%in%arg.names], collapse = ", "), "in dialog, but not in function.\n"))
#  }
#  
  if(length(items.idx)) for(jj in 1:length(items.idx)){
    markup <- list()    
    pos <- items.idx[jj] # position in values
    value <- dlg.list[[pos]] # the actual value that the dialog item has
    markup['value'] <- list(value) # to deal with NULL
    markup$signals <- list()
    kk <- pos+1 # index we're at in the dialog list
    while(!kk %in% items.idx && kk <= length(dlg.list)){
      if(names(dlg.list)[kk] == "signal"){
        markup$signals[[length(markup$signals)+1]] <- as.list(dlg.list[[kk]])
      } else {
        markup[[length(markup)+1]] <- dlg.list[[kk]]
        if(is.null(dlg.list[[kk]])) # for nulls
           markup[[length(markup)+1]] <- list(dlg.list[[kk]])
        names(markup)[[length(markup)]] <- names(dlg.list)[kk]
      }
      kk <- kk +1
     }
   markup$name <- arg.names[jj]
   markup$dlg.item.name <- dlg.item.names[jj]

   retval[[length(retval)+1]] <- markup
   names(retval)[length(retval)] <- names(dlg.list)[items.idx[jj]]
   }
  retval
}    

  # this function interprets the "signals" markup and performs signal.connect
  # appropriately
wire.up.items <- function(a, dlg.items, envir){

  locate.arg <- function(a, arg) {
    idx <- which(unlist(sapply(a, function(x) x$name==arg), recursive=F))
    if(length(idx) != 1) 
      stop(paste("locate.arg: argument named", 
         arg, "either not found or duplicated"))
    return(idx)
  }
  # Connect signals
  for(jj in 1:length(a)){
    for(sig in a[[jj]]$signals){
      if(length(sig) < 2) stop("Signal argument needs to contain at least 2 items")
      signl <- sig[[1]]      
      if(!is.character(signl)) stop("Signal argument 1 needs to be a string containing the signal")      
      user.data <- sig$user.data
      ww <- which(names(sig)%in%c("user.data"))
      if(length(ww) > 0) sig <- sig[-ww]
      sf <- sig[[2]]
      if(is.character(sf)) {
         # Look for signal function first in calling envir, then in local, then die
        if (exists(sf, envir=envir, mode="function")) {
          func <- get(sf, envir=envir)                   
        } else if (exists(sf, mode="function")) {
          func <- get(sf)                           
        } else {
          stop(paste("Can't find signaling function", sf))
        }
      } else if(is.function(sf)) {
        func <- sf
      } else {
        stop("Signal argument #2 needs to be a function name or function")
      }
      obj <- dlg.items[[names(a)[jj]]]
      data.args <- NULL
      if(length(sig) >= 3)
        data.args <- sapply(sig[(3:length(sig))], 
          function(x) dlg.items[[locate.arg(a, x)]])
      if(!is.null(user.data))
        data.args <- append(data.args, list(user.data=user.data))
      connect.args = list(signal = signl, obj=obj, f=func, data=data.args)
      do.call(signal.connect, connect.args)
    }
  }
    # Send the default signal from all items marked sensitive
      # This way if an item's marked insensitive, all its subordinates won't be 
       # toggled to sensitive
  for(item in dlg.items){
    if(item$isSensitive() && !identical(item$getData("signal.on.starting"), FALSE))
      tryCatch(signal(item), error = function(e){
        #cat(paste("Couldn't initialize dialog element", item$getData("name"), "to", deparse(get.value(item)), "\n")); flush.console()
        })
    }    
}
  
create.dialog.items <- function(a){
  dlg.items <- list()
  for(jj in 1:length(a)){
    obj.name <- names(a)[jj]
    if(obj.name == "main") 
      next;
    markup <- a[[jj]]
    
    if(is.call(markup$value)) 
      markup['value'] <- list(tryCatch({
        eval(markup$value)
        }, error = function(e) {                                                      
          #deparse(markup$value)
          #cat(paste("Couldn't set dialog element", markup$name, "to", deparse(markup$value), "\n"))
          #flush.console()
          return(NULL)
        }))
    dlg.item <- do.call(get(markup$dlg.item.name), markup)
    dlg.items[[obj.name]] <- dlg.item
    
    if(is.logical(markup$sensitive)) 
       dlg.items[[obj.name]]$setSensitive(markup$sensitive)    
    if(is.logical(markup$signal.on.starting)) # option to signal from initialize dialog item on starting
       dlg.items[[obj.name]]$setData("signal.on.starting", markup$signal.on.starting)

  }
  return(dlg.items)
}

create.panel <- function(a, dlg.items){
  # this function creates a box with appropriate elements from the 
  # processed markup list
  vbox.main <- gtkVBoxNew()

  if(!is.null(a$main$label)) 
    vbox.main$packStart(gtkLabelNew(a$main$label), FALSE, FALSE, 5)

  hbox0 <- gtkHBoxNew() # hbox0 is where we pack our vboxes
  curr.vbox <- gtkVBoxNew()
  #tooltips <- gtkTooltipsNew()


  # Assign items
  for(jj in 1:length(a)){
    obj.name <- names(a)[jj]
    if(obj.name == "main") 
      next;

    markup <- a[[jj]]

    if(identical(markup$visible, FALSE))
      next;
    
    dlg.item <- dlg.items[[obj.name]]
   
    nam <- markup$name
    lab.str <- nam
    if(!is.null(markup$label)) lab.str <- markup$label

    
    if(!is.null(markup$tooltip) && nchar(markup$tooltip)) 
      dlg.item$setTooltipText(markup$tooltip) 
        
    expand <- !is.null(dlg.item$getData("expand"))
    bin <- myFrame(dlg.item, lab.str)  
    # add indentation if requested
    if(!is.null(markup$indent))  {
      hbox <- gtkHBoxNew()
      vbox <- gtkVBoxNew()      
      hbox$packStart(vbox, FALSE, FALSE, markup$indent)
      hbox$packEnd(bin, TRUE, TRUE, 0)      
      curr.vbox$packStart(hbox, expand, TRUE, 0)    
    } else {
      curr.vbox$packStart(bin, expand, TRUE, 0)    
    }
      # go to the next column if you see a BREAK
    if("BREAK"%in%names(markup)) {
      hbox0$packStart(curr.vbox, TRUE, TRUE, 0)
      curr.vbox <- gtkVBoxNew()
    }
  }
    # To pack the last vbox
  if(!"BREAK"%in%names(markup)) hbox0$packStart(curr.vbox, TRUE, TRUE, 0)
  vbox.main$packStart(hbox0, TRUE, TRUE, 0) 
  return(vbox.main)
}

run.dialog <- function(func, 
   pos = -1,
   envir = as.environment(pos),
   dlg.list = tryCatch(get(paste(func.name, "dialog", sep="."), envir=envir), error = function(e) NULL),
   var.browser=NULL,
   parent.window=NULL, 
   auto.assign=TRUE, 
   do.long.running=FALSE,
   OK_handler = default.handler,
   output.name.rule = "append",
   output.name = NULL,
   do.logging = TRUE,
   log.handler = NULL,
   user.args = NULL, # the user can specify which args to use
   ...)
{  

   func.name <- NULL
   if(missing(func)) # return NULL if no function
     func <- function(...) {}
   if(is.character(func)) {
     func.name <- func   
     func <- get(func, envir=envir)
   }
   stopifnot(is.function(func))  
   if(is.null(func.name)) func.name <- deparse(substitute(func))
   if(is.null(dlg.list))  dlg.list=tryCatch(get(paste("", func.name, "dialog", sep="."), envir=envir),
      error = function(e) get(paste(func.name , "dialog", sep="."), envir=envir))
      
   stopifnot(is.list(dlg.list))
#return()    
      # turn the flattened markup list into a list with $main and one item per dialog item
    dlist <- process.markup(dlg.list, func)
    if(!missing(user.args)){
#      uas <- substitute(user.args)
        # turn any symbols into their deparsed characters, for objectItems
#      for(ii in 2:length(uas)) if(is.symbol(uas[[ii]])) uas[[ii]] <- deparse(uas[[ii]])
#      symbol2names.user.args <- uas
         # Make sure that you pass a string to this for objectItems, not a symbol
      for(ii in seq(length=length(dlist))){
        item <- dlist[[ii]] 
        if(isTRUE(item$name%in%names(user.args))){
          user.value <- user.args[[item$name]]
          if(item$dlg.item.name%in%c("radiobuttonItem", "choiceItem")){ # named value arg
            value.names <- names(item$value)
            if(any("value"%in%value.names)) names(value.names)["value"%in%value.names] <- ""
            if(any(item$value%in%user.value)) names(dlist[[ii]]$value)[which(item$value%in%user.value)] <- "value"
          } else if (item$dlg.item.name%in%"integerItem"){
            dlist[[ii]]$value["value"] <- user.value
          } else if (item$dlg.item.name%in%c("dataframeItem", "objectItem")){
            dlist[[ii]]$value <- user.value
          } else {
            dlist[[ii]]$value <- user.value
          }
        }
      }
    }
    #print(dlist)
      
    input.name <- NULL
    tryCatch({
		  # set hints if varbrowser is set
		if(!is.null(var.browser)){
		  item.hint <- names(which(sapply(dlist, function(x) 
		     !is.null(x$dlg.item.name) && x$dlg.item.name%in%c("dataframeItem", "objectItem"))))
		  nsr <- my.getToolkitWidget(var.browser)$getSelection()$countSelectedRows()
		  if(nsr > 0 && length(item.hint) > 0) {
		    vb.choice <- my.svalue(var.browser)
		    if(!identical(dlist[[item.hint[1]]]$take.hint, FALSE)) # take.hint markup cancels this
		      dlist[[item.hint[1]]]$value <- vb.choice
		  }}
    }, error = function(e){
      print("Can't find var browser to set")
      print(e)
     })      
      # create the list of dialog items
    dlg.items <- create.dialog.items(dlist)
      
      # create a box containing the dialog elements
    box <- create.panel(dlist, dlg.items)
      # set the signaling properties between dialog elements
    wire.up.items(dlist, dlg.items, envir)
  
    title <- func.name
    if(!is.null(dlist$main$title)) title <- dlist$main$title
    
    close.str <- "gtk-cancel"
    if(identical(PLATFORM_OS_TYPE, "windows") && identical(dlist$main$keep.open, TRUE)) close.str <- "gtk-close"    

    if(identical(dlist$main$ok.button, FALSE)) 
      dialog <- gtkDialog(title, parent.window, "modal", "gtk-close", 0, show = F)              
    else 
      dialog <- gtkDialog(title, parent.window, "modal", "gtk-ok", 1, close.str, 0, show = F)
  
    dialog$setPosition(GtkWindowPosition["center-on-parent"])
    dialog[["vbox"]]$packStart(box, TRUE, TRUE, 0)
    dialog$showAll()
    dialog_button_box <-  dialog$getChildren()[[1]]$getChildren()[[2]]
    dgc <- dialog$getChildren()[[1]]$getChildren()
    dialog_button_box <- dgc[[length(dgc)]]
    
    runMe <- function(){ 
      dialog_button_box$setSensitive(FALSE)      
      the.args <- lapply(dlg.items, get.value)
      names(the.args) <- sapply(dlg.items, get.name)
  
      dataset.name <- NULL
      if(length(the.args)) for(jj in 1:length(the.args)){
          # this is to get the name or call as required    
        if(dlist[[names(dlg.items)[jj]]]$dlg.item.name%in%c("dataframeItem", "objectItem") 
           && !identical(dlist[[names(dlg.items)[jj]]]$as.character, TRUE)){
          the.args[jj] <- list(tryCatch(as.call(parse(text=the.args[[jj]]))[[1]], error = function(e) NULL))        
          if(is.null(dataset.name) && !is.null(the.args[[jj]])){
            dataset.name <- deparse(the.args[[jj]])
          }
        }
      }
        # if the "suppress" markup is used, then suppress this value from being passed to the function
      suppress.idx <- which(unlist(sapply(dlist, function(x) !is.null(x$suppress) && identical(x$suppress, TRUE))))
      for(nam in names(suppress.idx)) the.args[[dlist[[nam]]$name]] <- NULL    
      
      if("..."%in%names(the.args)){
        ss <- paste("list(",the.args$...,")", sep="")
        s2 <- eval(parse(text=ss))
        for(nn in names(s2)) the.args[[nn]] <- s2[[nn]]
        the.args$... <- NULL
      }
      if(is.null(output.name)){
        if(output.name.rule == "append"){
      		if(is.null(dataset.name)){
      	    output.name <- paste(func.name, "output", sep="_")
      	  } else {
      			output.name <- paste(func.name, dataset.name, sep="_")
      	  }
     	  } else if (output.name.rule == "replace") {
    	    output.name <- dataset.name
  	    } else {
  	      stop("Unknown output.name.rule")
  	    }
    	  output.name <- gsub(".", "_", make.names(output.name), fixed=T)  	    
	    }
  	  #the.args <<- the.args
        # Execute the function
      retval <- NULL               
      if(.Platform$OS.type == "windows" && do.long.running || identical(dlist$main$long.running, TRUE)) {
        do.long.running.task(func=func, values=the.args, func.name=func.name, input.name=input.name, output.name=output.name)
      } else { 
      
        dcan <- NULL
        if(identical(dlist$main$show.progress, TRUE)){
          dcan = gtkDialog("Progress", parent.window, "modal", "gtk-cancel", 1, show=F)
          pb <- gtkProgressBar()
          dcan[["vbox"]]$packStart(pb, TRUE, FALSE, 10)        
          dcan$setPosition(GtkWindowPosition["center-on-parent"])
          dcanbb <- rev(dcan$getChildren()[[1]]$getChildren())[[1]]
          b = dcanbb$getChildren()[[1]]
          gSignalConnect(b, "clicked", function(button, data=list(dcan=dcan, dlg=dialog)){
            data$dcan$destroy()
            dialog_button_box$setSensitive(TRUE)
            #dlg$destroy()
            .C("rgtk2extras_interrupt")
          })
            # Pass the p.b. to the function if it's looking for one 
          pw_contains_bar <- FALSE
          pw_contains_label <- FALSE        
          if("progressbar" %in% names(formals(func)) && !"progressbar" %in% names(the.args)) { 
            pw_contains_bar <- TRUE
            the.args <- append(the.args, list(progressbar = pb))
          } 

          if("progresslabel" %in% names(formals(func)) && !"progresslabel" %in% names(the.args)) { 
            pw_contains_label <- TRUE
            pl <- gtkLabelNew()
            dcan[["vbox"]]$packStart(pl, TRUE, FALSE, 10)                  
            the.args <- append(the.args, list(progresslabel = pl))
          }                  
          if(!pw_contains_bar && !pw_contains_label){
            pb$setText("Function is running...")
          }
          dcan$show()                 
          while(gtkEventsPending()) gtkMainIteration()        
        }

   	    retval <- OK_handler(func=func, values=the.args, func.name=func.name, input.name=input.name, output.name=output.name, cancel.dialog=dcan)
#  	    while(TRUE) Sys.sleep(0.1)
        if(!is.null(dcan) && !inherits(dcan, "<invalid>")) dcan$destroy()    
  	    
      }
        # Assign the output. Can do it within hierarchical list or else create a new name.
  		if(!is.null(retval) && auto.assign){
        assign(output.name[1], retval, envir=.GlobalEnv)
      }
      # Log the action to a file
    if(!do.long.running && do.logging && !is.null(log.handler) && is.function(log.handler)){
      the.args$progresslabel <- NULL
      the.args$progressbar <- NULL
      sub.values <- as.character(substitute(the.args))
      cidx <- sapply(the.args, function(x) is.character(x) && length(x) == 1)
      sub.values[cidx] <- sapply(as.character(the.args[cidx]), deparse)
      cmd.string <- paste(func.name, "(",
        paste(paste(names(the.args), sub.values, sep="="), collapse=", "),
        ")\n", sep="")
      cmd.string <- paste(output.name, "<-", cmd.string)
      log.handler(cmd.string)
    }
    dialog_button_box$setSensitive(TRUE)    
    return(list(retval=retval, args = the.args))
  } # end runMe
  
  retval <- NULL
  if(identical(PLATFORM_OS_TYPE, "windows") && 
     identical(dlist$main$keep.open, TRUE)) {
    while(dialog$run() == 1){
#      tryCatch({
        retval <- runMe()    
#      }, interrupt = function() {
#        cat("Function interrupted.\n")
#      })
    }
    dialog$destroy()    
  } else {
    if(dialog$run() == 1){
      retval <- runMe()  
      dialog$destroy()        
    } else {
      dialog$destroy()
    }  
  }
  return(invisible(retval))

}

# This is our handler for computations and returning things
# Try to call the output "Dataset.Function" or "Function.Output"
# input name is the first dataset name we've got, if any
default.handler = function(func, values, func.name=NULL, input.name=NULL, output.name=NULL, data=NULL, 
    confirm.overwrite=F, cancel.dialog=NULL) {
  retval <- NULL
  stopifnot(is.function(func))
  tryCatch({
    retval <- do.call(func, values)        
	},
    interrupt = function(ex)
  {
    cat("Function was interrupted.\n");
    #print(ex);
  },
  error = function(ex)
  {
    #gmessage(paste("An error occurred in", func.name, "\n\n", as.character(ex)))
    if(!is.null(cancel.dialog) && !inherits(cancel.dialog, "<invalid>")) cancel.dialog$destroy()
    quick_message(paste("An error occurred in", func.name, "\n\n", as.character(ex)))
  },
  #warning = function(w)
  #{
  #  gmessage(paste("A warning occurred in", func.name, "\n\n", as.character(w)))
  #},
  finally =
  {
    if(!is.null(cancel.dialog) && !inherits(cancel.dialog, "<invalid>")) cancel.dialog$destroy()
  }) # tryCatch()
  return(retval)
}

##
##
# These are some functions for common dialog signaling needs

 # This function should trigger the same result as if the rightmost button were clicked
run.it <- function(item){
  while(!inherits(item$getParent(), "GtkWindow")) item <- item$getParent()
  rev(rev(item$getChildren())[[1]]$getChildren())[[1]]$clicked()
}


# Get names from the item and set them to list
get.names <- function(item, b, user.data=NULL) {
  #obj <- get(get.value(item))
  obj <- safe.eval(get.value(item))
  set.value(b, names(obj))
}

# user.data (list(extra.names = "row.names")) puts row.names at the start
get.colnames <- function(item, ..., user.data=NULL) {
  #obj <- get(get.value(item))
  if(object.exists(get.value(item))){
    obj <- safe.eval(get.value(item))
    cn <- colnames(obj)
    if(is.null(cn)) cn <- 1:dim(obj)[2]
    if(!is.null(user.data$extra.names)) cn <- c(user.data$extra.names, cn)
    sapply(list(...), function(x) set.value(x, cn))
  }
}

get.rownames <- function(item, ..., user.data=NULL) {
  #obj <- get(get.value(item))
  obj <- safe.eval(get.value(item))  
  rn <- colnames(obj)
  if(is.null(rn)) rn <- 1:dim(obj)[1]  
  if(!is.null(user.data$extra.names)) rn <- c(user.data$extra.names, rn)  
  sapply(list(...), function(x) set.value(x, rn))
}

# Get the data object name from b and the column from item and report 
# the column objects
get.column <- function(item, b, c, user.data=NULL) {
  colname <- get.value(item)[1]
  #obj <- get(get.value(b))
  obj <- safe.eval(get.value(b))  
  set.value(c, obj[,colname])
}

# Set all arguments to user.data
set.to <- function(item, ..., user.data=NULL){
    sapply(list(...), function(x) set.value(x, user.data))
  }

toggle.sensitive = function(item, ...){
    #x$setSensitive(get.value(item))
    sapply(list(...), gtkWidgetSetSensitive, get.value(item))    
      # Send on the signal from the item only if it's sensitive
        # This way you can have hierarchies of toggles for example
    sapply(list(...), function(x) if(x$isSensitive()) signal(x))
  }
      

  # Transfer selected list rows from the second argument to the item
push.selection = function(item, b) {    
  xx <- union(get.value(item), get.value(b, selected=T))
#  print("push.selection1")
#  print(xx)
  set.value(item, xx)
  xx <- intersect(get.value(item), get.value(b, selected=T))
  xx <- setdiff(get.value(b), xx)
#  print("push.selection2")  
#  print(xx)
  set.value(b, xx)
}

  # Transfer selected list rows from the item to the second argument
pop.selection = function(item, b) {    
  push.selection(b, item)
}

  # Copy selected list rows from the second argument to the item
  # Assumes all items are unique.
copy.selection = function(item, b) {    
  xx <- union(get.value(item), get.value(b, selected=T))
  set.value(item, xx)
}
  # Remove selected list rows
remove.selection = function(item) {
  xx <- setdiff(get.value(item), get.value(item, selected=T))
  set.value(item, xx)  
}

get.choice.names <- function(choice){
  get.choice <- safe.eval(choice, envir=.GlobalEnv)
#	if(exists(choice, envir=.GlobalEnv)){ # Selected a data frame
#		get.choice <- get(choice, envir=.GlobalEnv)
#	} else { # we've specified a data frame inside a list
#		parse.choice <- parse(text=choice)
#		get.choice <- eval(parse.choice, envir=.GlobalEnv)
#	}      
	stopifnot(!is.null(names(get.choice)))
	values <- names(get.choice)
	return(values)
}	

  # Utility function to return all tables in global
get_all_tables <- function(){
  .e <- new.env()
  .e$ll <- c()
  list_all_tables <- function(item, nam){ 
    if(is.data.frame(item)||is.matrix(item)){
      if(nam%in%.e$ll) warning(paste("Name shadowed", nam))
      .e$ll <- c(.e$ll, nam)
    }
    if(is.list(item) && length(item)>0)
      for(ii in 1:length(item)){
       if(!is.null(names(item))) {         
         list_all_tables(item[[ii]], paste(nam, names(item)[ii], sep="$"))
       }# else {
       #  list_all_tables(item[[ii]], paste(nam, "[[", ii, "]]", sep=""))
       #}
      }
  }
    for(item in ls(envir=.GlobalEnv)){
    list_all_tables(get(item), item)
  }
  return(.e$ll)
}


setClass("gMenuRGtk",
         contains="gComponentRGtk",
         prototype=prototype(new("gComponentRGtk"))
         )


## menulist is a list of lists with named components. Each named sub
## is a submenu.  a leaf consistis of handler= (required), lab


## override the defaults of Michael
gtkMenuPopupHack <-gtkMenuPopup
## gtkMenuPopupHack <-
## function (object, parent.menu.shell = NULL, parent.menu.item = NULL,
##     func = NULL, data = NULL, button, activate.time)
## {
##     checkPtrType(object, "GtkMenu")
##     if (!is.null(parent.menu.shell))
##         checkPtrType(parent.menu.shell, "GtkWidget")
##     if (!is.null(parent.menu.item))
##         checkPtrType(parent.menu.item, "GtkWidget")
##     if(!is.null(func))
##       func <- as.function(func)
##     button <- as.numeric(button)
##     activate.time <- as.numeric(activate.time)
##     w <- RGtk2:::.RGtkCall("S_gtk_menu_popup", object, parent.menu.shell,
##         parent.menu.item, func, data, button, activate.time,
##         PACKAGE = "RGtk2")
##     return(invisible(w))
## }

## gtkMenuPopupHack = function (object, parent.menu.shell=NULL,
##   parent.menu.item=NULL, func=NULL, 
##   data = NULL, button, activate.time) 
## {
##     checkPtrType(object, "GtkMenu")
## #    checkPtrType(parent.menu.shell, "GtkWidget")
## #    checkPtrType(parent.menu.item, "GtkWidget")
## #    func <- as.function(func)
##     button <- as.numeric(button)
##     activate.time <- as.numeric(activate.time)
##     w <- RGtk2:::.RGtkCall("S_gtk_menu_popup", object, parent.menu.shell, 
##         parent.menu.item, func, data, button, activate.time, 
##         PACKAGE = "RGtk2")
##     return(invisible(w))
## }


## put menu in group,
## a menubar is a map from a list into a menubar
## constructor
setMethod(".gmenu",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   menulist, 
                   popup = FALSE,
                   action = NULL,
                   container=NULL, ...) {
            
            force(toolkit)
            
            if(popup) {
              mb = gtkMenuNew()
            } else {
              mb = gtkMenuBarNew()
            }
            
            group = ggroup(spacing=0); svalue(group) <- 0
            mbgroup = ggroup(spacing=0); svalue(mbgroup) <- 0
            if(popup) {
              obj = new("gMenuRGtk", block=mb, widget=mb, toolkit=toolkit)
            } else {
              add(mbgroup, mb, expand=TRUE)
              add(group, mbgroup, expand=TRUE)

              obj = new("gMenuRGtk", block=group, widget=mb, toolkit=toolkit)
              
            }
  
            
            tag(obj, "menulist") <- menulist
            tag(obj, "action") <- action
            tag(obj,"popup") <- popup
            tag(obj, "mbgroup") <- mbgroup
            tag(obj, "mb") <- mb                  # the real menubar
            
  
            .addSubMenu(mb,menulist, action=action)
            
            
            if(!is.null(container)) {
              if(is.logical(container) && container == TRUE) {
                add(gwindow(visible=TRUE), obj)
              } else {
                add(container, obj,...)
              }
            }
            
            invisible(obj)
          })


### methods
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gMenuRGtk"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {
            tag(obj, "menulist")
          })

## three cases for value: list, gMenuRGtk, guiWidget push down
## make a menubar, then replace current -- isn't working for popup case
setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gMenuRGtk",
                           value="list"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   popup = tag(obj, "popup")
                   if(!is.null(popup) && popup == TRUE)
                     mb = gtkMenuNew()
                   else
                     mb = gtkMenuBarNew()
                   mbgroup = ggroup(spacing=0); svalue(mbgroup) <- 0
                   add(mbgroup, mb, expand=TRUE)
                   
                   menulist = value           # value is a list
                   if(!is.list(menulist))
                     stop("value is not a menubar or a list")
                   
                   .addSubMenu(mb,menulist, action=tag(obj,"action"))
                   
                   delete(obj@block, tag(obj,"mbgroup")) # delete from group()
                   add(obj@block, mbgroup, expand=TRUE) # add to block
                   tag(obj,"mbgroup") <- mbgroup
                   tag(obj,"menulist") <- menulist

                   return(obj)
                 })

## get list, and then call previous
setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gMenuRGtk",
                           value="gMenuRGtk"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   .svalue(obj,toolkit, index, ...) <- svalue(value)
                   return(obj)
                 })
## call previous after getting list
setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gMenuRGtk",
                           value="guiWidget"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   .svalue(obj,toolkit,index, ...) <- svalue(value@widget)
                   return(obj)
                 })

setMethod(".add",
          signature(toolkit="guiWidgetsToolkitRGtk2", obj="gMenuRGtk", value="guiWidget"),
          function(obj, toolkit,  value, ...) {
            .add(obj, toolkit, value@widget)
          })

setMethod(".add",
          signature(toolkit="guiWidgetsToolkitRGtk2", obj="gMenuRGtk", value="gMenuRGtk"),
          function(obj, toolkit,  value, ...) {
            orig.list = svalue(obj)
            add.list = svalue(value)
            new.list = c(orig.list, add.list)
            svalue(obj) <- new.list
          })


setMethod(".add",
          signature(toolkit="guiWidgetsToolkitRGtk2",
                    obj="gMenuRGtk", value="list"),
          function(obj, toolkit,  value, ...) {
            orig.list = svalue(obj)
            new.list = c(orig.list, value)
            svalue(obj) <- new.list
          })

## "wdget" is either a gMenu, list or just names to delete
setMethod(".delete",
          signature(toolkit="guiWidgetsToolkitRGtk2", obj="gMenuRGtk",
                    widget="guiWidget"),
          function(obj, toolkit, widget, ...) {
            .delete(obj,toolkit,widget@widget,...)
          })
setMethod(".delete",
          signature(toolkit="guiWidgetsToolkitRGtk2", obj="gMenuRGtk",
                    widget="gWidgetRGtk"),
          function(obj, toolkit, widget, ...) {
            .delete(obj,toolkit,widget@widget, ...)
          })
setMethod(".delete",
          signature(toolkit="guiWidgetsToolkitRGtk2", obj="gMenuRGtk",
                    widget="gMenuRGtk"),
          function(obj, toolkit, widget, ...) {
            .delete(obj,toolkit,svalue(widget), ...)
          })

setMethod(".delete",
          signature(toolkit="guiWidgetsToolkitRGtk2", obj="gMenuRGtk",
                    widget="list"),
          function(obj, toolkit, widget, ...) {
            lst = widget                    # else assume its a character
            
            cur.list = svalue(obj)
            for(i in lst) {
              ## we delete *last* entry with this name, hence this awkwardness
              theNames = names(cur.list)
              if(i %in% theNames) {
                j = max(which(i == theNames))
                if(!is.null(cur.list[[j]])) cur.list[[j]] <- NULL
              }
            }
            ## now update menubar
            svalue(obj) <- cur.list
          })

## give vector notation
setMethod("[",
          signature(x="gMenuRGtk"),
          function(x, i, j, ..., drop=TRUE) {
            .leftBracket(x, x@toolkit, i, j, ..., drop=drop)
          })
setMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gMenuRGtk"),
          function(x, toolkit, i, j, ..., drop=TRUE) {
            lst = svalue(x)
            if(missing(i))
              return(lst)
            else
              return(lst[i])
          })

setReplaceMethod("[",
                 signature(x="gMenuRGtk"),
                 function(x, i, j,..., value) {
                   .leftBracket(x, x@toolkit, i, j, ...) <- value
                   return(x)
                 })

setReplaceMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gMenuRGtk"),
          function(x, toolkit, i, j, ..., value) {
            lst = svalue(x)
            theNames = names(lst)
            if(is.character(i))
              i = max(which(i %in% theNames))
            lst[[i]] <- value[[1]]
            theNames[i] = names(value)
            names(lst) = theNames
            svalue(x) <- lst
            return(x)
          })

###

## some helper functions for this
.isLeaf = function(lst) {
  if(.isgAction(lst) ||
     (is.list(lst) & (!is.null(lst$handler) | !is.null(lst$separator)))
     ) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

## workhorse for this
.addSubMenu = function(subMenu, menu.list, action=NULL, ...) {
  
  for(i in names(menu.list)) {
    data = menu.list[[i]]

    if(.isgSeparator(data)) {
      data <- list(separator = TRUE)
    }
      
    
    if(.isgAction(data)) {
      action <- getWidget(data)
      item <- gtkImageMenuItem("")
      if("always-show-image" %in% names(item))
        item['always-show-image'] <- TRUE
      subMenu$Append(item)
      ##action$connectProxy(item)
      item$setRelatedAction(action)
    }  else if(!.isLeaf(data)) {
      ## do submenu
      item = gtkMenuItem(i)
      subMenu$Append(item)
      newSubMenu =  gtkMenu()
      .addSubMenu(newSubMenu, data)
      item$SetSubmenu(newSubMenu)
    } else if(!is.null(data$separator)) {
      ## add a separator
      item = gtkSeparatorMenuItem()
      subMenu$Append(item)
    } else {
      ## what name
      if(!is.null(data$label))
        theName = data$label
      else
        theName = i
      ## make a menuitem
      item <- gtkImageMenuItemNewWithLabel(theName)
      if(!is.null(data$icon)) {
        icon = data$icon
        if(file.exists(icon)) {
          ## a file on system
          image = gtkImageNewFromFile(icon)
        } else {
          ## assume a stock icon file
          icon = getstockiconname(icon)
          image = gtkImageNew()
          image$SetFromStock(icon,size=GtkIconSize["menu"])
        }
        item$SetImage(image)
        if("always-show-image" %in% names(item)) # newer GTK
          item['always-show-image'] <- TRUE
      }  else {
        item = gtkMenuItem(theName)
      }
      subMenu$Append(item)
      item$AddCallback("activate",data$handler, data=list(action=action))
    } 
  }
}


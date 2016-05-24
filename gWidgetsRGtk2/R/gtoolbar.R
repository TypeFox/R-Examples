## gtoolbar, similar to gmenu

setClass("gToolbarRGtk",
         contains="gComponentRGtk",
         prototype=prototype(new("gComponentRGtk"))
         )

## turn a list into a uimgr object
setMethod(".gtoolbar",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   toolbarlist,
                   style = c("both","icons","text","both-horiz"),
                   action=NULL,
                   container=NULL, ...) {

            force(toolkit)
            
            style = match.arg(style)
            toolbar = .mapListToToolBar(toolbarlist, style)


            group = ggroup(spacing=0)
            svalue(group) <- 0                    # border
            add(group, toolbar, expand=TRUE)

            obj = new("gToolbarRGtk",block=group, widget=toolbar, toolkit=toolkit)

            tag(obj,"toolbarlist") <- toolbarlist
            tag(obj,"toolbar") <- toolbar
            tag(obj,"style") <- style
            tag(obj,"group") <- group
  
            ## attach to container
            if (!is.null(container)) {
              if(is.logical(container) && container == TRUE)
                container = gwindow(visible=TRUE)
              add(container, obj, ...)  # was expand=TRUE
            }

            ## warn
            if(!is.null(action))
              warning("The action argument is not yet defined for gtoolbar.")
            
            invisible(obj)
  
          })


.mapListToToolBar = function(lst, style, ...) {
  
  ## some helper functions for this
  is.leaf = function(lst) {
    if(.isgAction(lst) ||
       .isgSeparator(lst) ||
       (is.list(lst) & (!is.null(lst$handler) | !is.null(lst$separator)))
       ) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }

  ### XXX ### CHECK THIS
  toolbar <- gtkToolbar()
  ## toolbar has no recurse
  for(i in names(lst)) {
    itemData <- lst[[i]]
    if(.isgSeparator(itemData))
      itemData <- list(separator = TRUE)
    if(.isgAction(itemData)) {
      action <- getToolkitWidget(itemData)
      item <- action$createToolItem()
      toolbar$insert(item,pos = -1)
    } else if(is.list(itemData)) {
      if(!is.null(itemData$separator)) {
        ## XX put in separator into toolbar
        item <- gtkSeparatorToolItemNew()
        toolbar$insert(item,pos = -1)
      } else if(!is.null(itemData$handler)) {
        ## XXX put in value from itemData (label handler, ...)
        ## JV: need to call this constructor with label, otherwise gtk error pops up
        item <- gtkToolButtonNew(label=itemData$icon)
        if(!is.null(itemData$icon))
          item$setStockId(getstockiconname(itemData$icon))
        

        gSignalConnect(item,signal="clicked",
                       f = function(a,...) {
                         handler <- a$handler
                         action <- a$action
                         h <- list(action = action)
                         handler(h,...)
                       },
                       data = itemData,
                       user.data.first=TRUE
                       )
        toolbar$insert(item,pos = -1)
      }
    } else if(is(itemData,"guiWidget")) {
      ## can add in a button or popup et
      widget <- getBlock(itemData)
      item <- gtkToolItemNew()
      item$Add(widget)
      toolbar$insert(item,pos = -1)
    }
  }
    
  return(toolbar)
}
### main function returns toolbar from list
## replaced by above -- simpler, uses older API
.mapListToToolBar.old = function(lst, style, ...) {
  
  ## some helper functions for this
  is.leaf = function(lst) {
    if(is.list(lst) & (!is.null(lst$handler) | !is.null(lst$separator))) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  
  quiet.cb = function(h,...) {}
  Cat = function(..., file="", append=FALSE) {
    cat(Paste(...),file=file,append=append)
  }
  
  
  ## This function is called recursively to make both the actions and
  ## the UI.  The structure of lst is fairly basic -- lists of
  ## lists. The terminal nodes contain components name, icon, label, and
  ## handler. Null values are okay for all but name, handler.
  
  ## assigns to str (a string), and actions (a list) which must be defined
  ## previously
  make.ui = function(lst,name="root",no=1) {
    for(i in names(lst)) {
      if(!is.leaf(lst[[i]])) {
        ## make item
        Cat(Paste(rep("\t",no-1)),"<toolitem action=\"",i,"\">\n",file=filename,append=TRUE)
        ## add action, tooltip
        actions[[length(actions)+1]] <<- c(i,"gtk-null",i,"","",quiet.cb)
        
        ## call recursively
        make.ui(lst[[i]],i,no+1)
        ## close tag
        Cat(Paste(rep("\t",no-1)),"</toolbar>\n",file=filename,append=TRUE)
      } else {
        ## UImgr
        if(!is.null(lst[[i]]$separator)) {
          ## add a separator here
        Cat("<separator/>\n",file=filename, append=TRUE)
      } else {        
        Cat(Paste(rep("\t",no)),"<toolitem action=\"",i,"\"/>\n",file=filename,append=TRUE)
        
        ## actions
        tmp = lst[[i]]
        ## add action
        ## fill in missing values. Only a handler is needed
          if(is.null(tmp$icon)) tmp$icon="null"
          if(is.null(tmp$label)) tmp$label=i
          if(is.null(tmp$handler)) {
            tmp$handler = function(h,...) gwCat(i,"in toolbar list needs a handler")
          }
        ## tooltip
        if( is.null(tmp$tooltip))
          tmp$tooltip = tmp$label
        
          if(!is.function(tmp$handler)) {
            ## call using ggenericwidget
            if(is.character(tmp$handler)) {
              dalst = get(tmp$handler)        # called within closure
              tmp$handler = function(...)
                do.call("ggenericwidget",list(dalst))
            } else if(is.list(tmp$handler)) {
              dalst = tmp$handler             # call within closure
              tmp$handler = function(...)
                do.call("ggenericwidget",dalst)
            }
          }
        }

        ## now add to  actions
        actions[[length(actions)+1]] <<-
          list(name=i,
               "stock_id" = getstockiconname(tmp$icon),
               label = tmp$label,
               accelerator = "",
               tooltip = tmp$label,
               callback = tmp$handler) 
      }
    }
  }
  
  

  ## to create the UI and actions
  ## initialize
  filename = tempfile()
  
  actions = list()
  Cat("<ui>\n<toolbar name=\"","toolbar","\">\n", file=filename,append=FALSE)
  ## call function
  make.ui(lst)
  ## finish
  Cat("</toolbar>\n</ui>\n", file=filename,append=TRUE)

  acgrp = gtkActionGroupNew(name="ActionGroup")
  acgrp$AddActions(entries=actions) 
  uimgr = gtkUIManagerNew()
  uimgr$InsertActionGroup(acgrp,0)
  uimgr$AddUiFromFile(filename)

  unlink(filename)

  ## menubar
  toolbar = uimgr$GetWidget(Paste('/',"toolbar"))
  toolbar$setStyle(style)

  return(toolbar)

  
}

### methods
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gToolbarRGtk"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {
            tag(obj, "toolbarlist")
          })

setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gToolbarRGtk"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   if(!is.list(value)) 
                     stop("A toolbar requires a list to define it.")

                   toolbar = .mapListToToolBar(value, tag(obj,"style"))
                   ## swap out
                   delete(tag(obj,"group"), tag(obj,"toolbar") )
                   add(tag(obj,"group"), toolbar, expand=TRUE)
                   ## replace
                   tag(obj,"toolbar") <- toolbar
                   tag(obj,"toolbarlist") <- value
                   
                   ##  all done
                   return(obj)
                 })

## returns list, or part of list
setMethod("[",
          signature(x="gToolbarRGtk"),
          function(x, i, j, ..., drop=TRUE) {
            .leftBracket(x, x@toolkit, i, j, ..., drop=drop)
          })
setMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gToolbarRGtk"),
          function(x, toolkit, i, j, ..., drop=TRUE) {
            lst = tag(x,"toolbarlist")
            if(missing(i))
              return(lst)
            else
              return(lst[[i]])
          })

setReplaceMethod("[",
                 signature(x="gToolbarRGtk"),
                 function(x, i, j,..., value) {
                   .leftBracket(x, x@toolkit, i, j, ...) <- value
                   return(x)
                 })

setReplaceMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gToolbarRGtk"),
          function(x, toolkit, i, j, ..., value) {
            if(!is.list(value))
              stop("assignment must be a list defining a (part) of a toolbar.")
            lst = tag(x,"toolbarlist")
            if(missing(i))
              lst = value
            else
              lst[[i]] = value
            
            svalue(x) <- lst
            
            return(x)
          })


setMethod(".add",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gToolbarRGtk", value="list"),
          function(obj, toolkit, value,  ...) {
            svalue(obj) <- c(svalue(obj), value)
          })

## (from gmenu)
setMethod(".delete",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gToolbarRGtk"),
          function(obj, toolkit, widget,  ...) {
            ## widget can be gToolBar or a list
            if(is.character(widget)) {
              lst = widget                    # else assume its a character
            } else if(is(widget,"gComponentRGtk")) {
              lst = svalue(widget)
              lst = names(lst)
            } else if(is.list(widget)) {
              lst = names(widget)
            } else {
              warning("Must be either a vector of names, a list, or a gToolbar instance")
              return()
            }
            
            cur.list = svalue(obj)             
            for(i in lst) {
              ## we delete *last* entry with this name, hence this awkwardness
              theNames = names(cur.list)
              if(i %in% theNames) {
                j = max(which(i == theNames))
                if(!is.null(cur.list[[j]])) cur.list[[j]] <- NULL
              }
            }
            ## now update toolbar
            svalue(obj) <- cur.list
          })

### no method to set style, use tag(obj,"style")<-"theStyle" instead

## gtoolbar, similar to gmenu
## need to incorporate delete/add methods as in imenu


setClass("gToolbartcltk",
         representation = representation("gComponenttcltk",
           style="character"),
         contains="gComponenttcltk",
         prototype=prototype(new("gComponenttcltk"))
         )

## turn a list into a uimgr object
setMethod(".gtoolbar",
          signature(toolkit="guiWidgetsToolkittcltk"),
          function(toolkit,
                   toolbarlist,
                   style = c("both","icons","text","both-horiz"),
                   action=NULL,
                   container=NULL, ...) {

            force(toolkit)

            if(is(container,"logical") && container)
              container <- gwindow("Toolbar")
            if(!is(container,"guiWidget")) {
              warning("Container is not correct. No NULL containers possible\n" )
              return()
            }


            style = match.arg(style)
            tkstyle <- c("both"="top",
                         "icons"="image",
                         "text"="text",
                         "both-horiz"="left")


            ## container must be a gwindow
            toplevel <- getTopLevel(container)
            
            if(!(is(toplevel,"gWindowtcltk") || is(toplevel@widget,"gWindowtcltk"))) {
              message(gettext("gtoolbar: container must be gwindow instance\n"))
            }
##            tt <- getBlock(container)
            tt <- tag(toplevel, "tb")
            gp <- ttkframe(tt)
            
            tb <- ttkframe(gp)
            tkpack(tb, side="left",anchor="w", expand=TRUE, fill="x")

            
            obj <- new("gToolbartcltk",block=gp, widget=tb,
              toolkit=toolkit, ID=getNewID(),e = new.env(),
              style=style)

            tag(obj,"toolbarlist") <- toolbarlist

            add(container, obj, ...)

            .mapListToToolBar(tb, toolbarlist, tkstyle[style])

            invisible(obj)
  
          })


## helpers
.addToolbarButton <- function(tb, style, label=NULL, icon=NULL,handler=NULL, action=NULL) {
  ## get icon
  if(!is.null(icon)) {
    file <- findTkIcon(icon)
    icon <- tcl("image","create","photo",file=file)
    ## make a button with icon
    b <- ttkbutton(tb, image=icon, text=label, compound=style)
  } else {
    b <- ttkbutton(tb, text=label)
  }
  ## add in handler
  handler = force(handler)              # need to force so scoping works in this call
  if(!is.null(handler)) {
    tkbind(b,"<Button-1>", function(...) {
      h = list(obj=b, action=action)
      handler(h,...)
    })
  }

##   slaves <- tclvalue(tcl("grid","slaves",tb))
##   slaves <- unlist(strsplit(slaves," "))
##   n <- length(slaves)
##   tkgrid(b, row=0, column=n, sticky="ns")
  
  tkpack(b, side="left",anchor="w",expand=TRUE,fill="y")
  return(b)
}


.mapListToToolBar = function(tb, lst, style) {
  ## list is simple compared to menubar
  for(i in names(lst)) {
    tmp <- lst[[i]]
    label <- i

    if(.isgSeparator(tmp))
      tmp <- list(separator=TRUE)
    ## is it a gaction?
    if(.isgAction(tmp)) {
      tmp <- getToolkitWidget(tmp)
      label <- tmp$label
    }
    
    
    if(!is.null(tmp$separator)) {
      ## add separator
      gseparator(horizontal=FALSE, container=tb)
    } else if(!is.null(tmp$handler)) {
      ## how to decide there are no text parts?
      b <- .addToolbarButton(tb, style, label, tmp$icon, tmp$handler, tmp$action)
      if(.isgAction(lst[[i]])) {
        if(is(lst[[i]],"gActiontcltk"))
          e <- lst[[i]]@e
        else
          e <- lst[[i]]@widget@e
        l <- e$toolbaritems
        l[[length(l) + 1]] <- b
        e$toolbaritems <- l
      }
      
    }
  }
}


### methods
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gToolbartcltk"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {
            tag(obj, "toolbarlist")
          })

setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gToolbartcltk"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   if(!is.list(value)) 
                     stop("A toolbar requires a list to define it.")

                   toolbar = obj@widget
                   ## delete from toolbar
                   n = length(tag(obj,"toolbarlist"))

                   ## how to delete from group
                   gwCat(gettext("No method to delete toolbar components\n"))
                   
                   .mapListToToolBar(toolbar, value, obj@style)

                   tag(obj,"toolbarlist") <- value
                   
                   ##  all done
                   return(obj)
                 })

## returns list, or part of list
setMethod("[",
          signature(x="gToolbartcltk"),
          function(x, i, j, ..., drop=TRUE) {
            .leftBracket(x, x@toolkit, i, j, ..., drop=drop)
          })
setMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkittcltk",x="gToolbartcltk"),
          function(x, toolkit, i, j, ..., drop=TRUE) {
            lst = tag(x,"toolbarlist")
            if(missing(i))
              return(lst)
            else
              return(lst[[i]])
          })

setReplaceMethod("[",
                 signature(x="gToolbartcltk"),
                 function(x, i, j,..., value) {
                   .leftBracket(x, x@toolkit, i, j, ...) <- value
                   return(x)
                 })

setReplaceMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkittcltk",x="gToolbartcltk"),
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
          signature(toolkit="guiWidgetsToolkittcltk",obj="gToolbartcltk", value="list"),
          function(obj, toolkit, value,  ...) {
            svalue(obj) <- c(svalue(obj), value)
          })

## (from gmenu)
setMethod(".delete",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gToolbartcltk"),
          function(obj, toolkit, widget,  ...) {
            ## widget can be gToolBar or a list
            if(is.character(widget)) {
              lst = widget                    # else assume its a character
            } else if(is(widget,"gComponenttcltk")) {
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

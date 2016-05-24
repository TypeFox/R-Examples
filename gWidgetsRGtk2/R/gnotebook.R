## class previously defined
setMethod(".gnotebook",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   tab.pos = 3,                          # same as pos= in text
                   closebuttons = FALSE,
                   dontCloseThese = NULL,                 # integer of tabs not to close
                   container=NULL,                           # add to this container
                   ...) {
            
            force(toolkit)

            ## beginvv
            notebook = gtkNotebookNew()
            notebook$SetScrollable(TRUE)

            
            ## tab placement: 1,2,3,4 -> 3,0,2,1
            types = c("bottom","left","top","right")
            tabposition = GtkPositionType[types]
            notebook$SetTabPos(tabposition[tab.pos])

            ## add close button, in same level as tab.pos
            if(tab.pos == 1 || tab.pos == 3) 
              group = ggroup(container=container, ...)
            else
              group = ggroup(horizontal=FALSE,container=container, ...)
            add(group,notebook,expand=TRUE)
            
            ## create gnotebook object
            obj = new("gNotebookRGtk", block=group, widget=notebook, toolkit=toolkit)

            tag(obj,"closebuttons") <- closebuttons
            tag(obj,"dontCloseThese") <- dontCloseThese
            


            invisible(obj)
          })

as.gWidgetsRGtk2.GtkNotebook <- function(widget,...) {
  ## no group here
  obj <- new("gNotebookRGtk", block=widget, widget=widget,
    toolkit=guiToolkit("RGtk2"))

  tag(obj,"closebuttons") <- FALSE
  tag(obj,"dontCloseThese") <- FALSE

  return(obj)
}

### methods

## different, set notebook, not group
setReplaceMethod(".size",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gNotebookRGtk"),
                 function(obj, toolkit, ..., value) {
                   width = value[1]; height = value[2]
                   obj@widget$SetSizeRequest(width,height)
                   return(obj)
                 })

## ## close buttons? Call this by default
## defaultCloseButtonHandler = function(h,...) {
##   notebook = h$action@notebook             # gtk notebook, not gnotebook
##   currentPage = notebook$GetCurrentPage()
##   notebook$RemovePage(currentPage)
##   svalue(h$action)  <- currentPage
## }



## return the current tab
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gNotebookRGtk"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {
            if(!is.null(index)) {
              warning("No index argument for a gnotebook instance")
            }
            notebook = obj@widget
            return(notebook$GetCurrentPage() + 1)
          })

## set the current tab to value
setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gNotebookRGtk"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   notebook = obj@widget
                   nPages = notebook$GetNPages()
                   notebook$SetCurrentPage(min(nPages,as.numeric(value)-1))

                   return(obj)
                 })

## ## set label or page values
## ## avlues is list with components page.no, label, and/or page
## ## setting page is not implemented (add, remove?)
## set.values.gNotebook = function(obj,values,...) {
##   page.no = values$page.no
##   if(is.null(page.no))
##     page.no = svalue(obj)

##   if(!is.null(values$label)) {
##     labelgroup = ggroup()

##     label = values$label
##     if(is.character(label)) {
##       label = glabel(label)
##     }
##     add(labelgroup, label)
##     ## label should be glabel instance
##     ## this code could be consolidated. It is taken from add.gNotebook
##     if(obj$closebuttons) {
##       add(labelgroup,label)
##       closeImage = gimage("gtk-close",dirname="stock",
##         handler = function(h,...) {
##           curPage = notebook$notebook$GetCurrentPage()
##           if(!is.null(notebook$dontCloseThese) &&
##              !((curPage + 1)  %in% notebook$dontCloseThese)) {
##             dispose(notebook)
##           }
##         })
##       add(labelgroup,closeImage, expand=FALSE)
##     }
##     obj$notebook$SetTabLabel(obj$notebook[[page.no]],labelgroup$ref)
##   }
##   ## now update page widget
##   if(!is.null(values$page)) {
##     warning("Setting a page in gNotebook is note implemented. Try adding, disposing")
##   }
## }


## remove the current tab
## this should be called delete -- which is used to remove objects
setMethod(".dispose",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gNotebookRGtk"),
          function(obj, toolkit,  ...) {
            theArgs = list(...)
            to.right=ifelse(!is.null(theArgs$to.right), theArgs$to.right,FALSE)
            dontCloseThese = tag(obj,"dontCloseThese")
            cur.page = svalue(obj)
  
            if(to.right) {
              no.pages = length(obj)
              no.right = no.pages - cur.page
    
              if(no.right > 0) {
                ## clear out, must work from last backwards
                for(i in no.right:1) {
                  if(!is.null(dontCloseThese) &&
                     !((cur.page - 1 + i + 1) %in% dontCloseThese)) {

                    ## destroy widget
                    theWidget = obj@widget$getNthPage(cur.page - 1 + 1)
                    obj@widget$RemovePage(cur.page - 1 +i) # cur.page 1-based
                    gtktry(theWidget$destroy())
                    
                    svalue(obj) <- cur.page
                  }
                }
              }
            } else {
              ## just this page
              if(!is.null(dontCloseThese)) {
                if(!((cur.page - 1 + 1) %in% dontCloseThese)) {
                  theWidget = obj@widget$getNthPage(cur.page - 1)                  
                  obj@widget$RemovePage(cur.page - 1) # value is 1 based, not 0
                  gtktry(theWidget$destroy())
                  svalue(obj) <- cur.page
                }
              } else {
      ## no restriction of closing page
                theWidget = obj@widget$getNthPage(cur.page - 1)                                  
                obj@widget$RemovePage(cur.page - 1) # value is 1 based, not 0
                gtktry(theWidget$destroy())
                svalue(obj) <-cur.page
              }
            }
          })

## remove the widget form the notebook
setMethod(".delete",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gNotebookRGtk"),
          function(obj, toolkit, widget,  ...) {
            obj@widget$remove(getWidget(widget))
          })


### add() is a workhorse method here. Several args available in ...
#add.gNotebook = functionf(obj, value,
#  label="", markup = FALSE,
#  index = NULL, override.closebutton = FALSE, ...) {
setMethod(".add",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gNotebookRGtk",
                    value="guiWidget"),
          function(obj, toolkit, value,  ...) {
            .add(obj, toolkit, value@widget, ...)
          })

setMethod(".add",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gNotebookRGtk",
                    value="gWidgetRGtk"),
          function(obj, toolkit, value,  ...) {
            ## in ... we have many possibilies
            ## label -- for setting label  (also look for name)
            ## index for setting the index of page to add
            ## markup -- markup label
            ## override.closebutton -- to not put closebutton even if set in constructor
            
            ## process ...
            theArgs = list(...)                      # for making generic
            if(!is.null(theArgs$label)) {
              label = theArgs$label
            } else if(!is.null(theArgs$name)) {
              label = theArgs$name
            } else {
              label = id(obj)
              if(is.null(label))
                label = "unnamed"
            }
            
            index = if (is.null(theArgs$index)) NULL else theArgs$index
            if(!is.null(theArgs$pageno)) index = theArgs$pageno # also called paegno
            markup = if (is.null(theArgs$markup)) FALSE  else theArgs$markup
            override.closebutton =
              if (is.null(theArgs$override.closebutton))
                FALSE
              else
                as.logical(theArgs$override.closebutton)

            
            ## let's go
            gtknotebook = obj@widget

            ## label
            labelgroup = ggroup()
            if(is.character(label))
              label = glabel(label, markup = markup)
            add(labelgroup,label)

            ## closebutton
            if(!is.null(tag(obj,"closebuttons")) &&
               as.logical(tag(obj,"closebuttons")) &&
               !override.closebutton) {
              closeImage = gimage("gtk-close",dirname="stock",
                handler = function(h,...) {
                  dispose(obj)
                  return(TRUE)
#                  cat("DEBUG: all done\n")
#                  curPage = gtknotebook$GetCurrentPage()

#                 if(is.null(tag(obj,"dontCloseThese")) || # dont close if in dontCloseThes
#                     (!is.null(tag(obj,"dontCloseThese")) &&
#                      !((curPage + 1)  %in% tag(obj,"dontCloseThese"))
#                      )) {
#                    gtknotebook$RemovePage(gtknotebook$GetCurrentPage())
#                    svalue(obj) <- curPage+1
                                        #                }
                })
              add(labelgroup,closeImage, expand=FALSE)
            }
            
            ## store widget into label for dnd
            tag(labelgroup,"widget") <- value

  
            group =ggroup()
            add(group, value, expand=TRUE)        # get to expand
            ## what to add
            gtkpage = getWidget(group)
            labelWidget = getBlock(labelgroup)
            ## where to add
            if(is.null(index) | !is.numeric(index)) {
              thePage = gtknotebook$AppendPage(gtkpage, labelWidget) 
            } else if(index < 0) {
              thePage = gtknotebook$PrependPage(gtkpage, labelWidget)
            } else {
              thePage = gtknotebook$InsertPage(gtkpage,
                labelWidget, position=(index - 1))
            }

            ## Add DND actions for labels
            theLabel = gtknotebook$GetTabLabel(gtknotebook$GetNthPage(thePage))
            adddropsource(theLabel,
                          targetType = "object",
                          handler = function(h,...) {
                            ## dump the widget attached to label
                          },
                          action = value
                          )
            
 
  ## add drop motion for labels
            adddroptarget(theLabel)
            adddropmotion(theLabel,handler = function(h,...) gtknotebook$SetCurrentPage(thePage))
            
            ## uncomment below, and comment above, to change drop motion
            ## add drop motion to label if no close button
            ##   if(!obj$closebuttons || override.closebutton) {
            ##     theLabel = obj$notebook$GetTabLabel(obj$notebook$GetNthPage(thePage))
            ##     adddroptarget(theLabel)
            ##     adddropmotion(theLabel,handler = function(h,...) obj$notebook$SetCurrentPage(thePage))
            ##   }
            
            ## move to newpage
            svalue(obj) <- thePage + 1
          })
            
## Regular R methods treat gnotebook like a vector

## find out number of pages
setMethod(".length",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gNotebookRGtk"),
          function(x, toolkit) {
            x@widget$GetNPages()
          })

## return tabnames
setMethod(".names",signature(toolkit="guiWidgetsToolkitRGtk2",x="gNotebookRGtk"),
          function(x, toolkit) {
            notebook = x@widget
            NPages = notebook$GetNPages()
            if(NPages == 0) {
              return(c())
            } else {
              theNames = sapply(1:NPages, function(i) {
                notebook$GetTabLabel(notebook$GetNthPage(i-1))[[1]][[1]]$GetText()
              })
              return(theNames)
            }
          })

## can assigne with names(x) <-x or even names(x)[i] <- "single name"
setReplaceMethod(".names",
                 signature(toolkit="guiWidgetsToolkitRGtk2",x = "gNotebookRGtk"),
                 function(x,toolkit, value) {
                   n = length(x)
                   if(length(value) != n)
                     stop("New names for notebook must have proper length")
                   
                   notebook = x@widget
                   
                   NPages = notebook$GetNPages()
                   if(NPages == 0) {
                     return(c())
                   } else {
                     for(i in 1:NPages)
                       notebook$GetTabLabel(notebook$GetNthPage(i-1))[[1]][[1]]$SetText(value[i])
                   }
                   invisible(x)
                 })


## return widget contained in notebook page i as a  list or single widget
setMethod("[",
          signature(x="gNotebookRGtk"),
          function(x, i, j, ..., drop=TRUE) {
            .leftBracket(x, x@toolkit, i, j, ..., drop=drop)
          })
setMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gNotebookRGtk"),
          function(x, toolkit, i, j, ..., drop=TRUE) {
            if(missing(i))
              i = seq_along(x)

            i <- i[i <= length(x)]
            i <- i[i > 0]

            if(length(i) == 0) {
              warning("No widget for that index")
              return(NULL)
            }
           

            if(length(i) > 1) {
              lst = lapply(i,function(j)
                getNotebookPageWidget(x,pageno = j-1)
                )
              return(lst)
            } else {
              return(getNotebookPageWidget(x, pageno = i-1))
            }
          })


## Puts widget into a position
setReplaceMethod("[",
                 signature(x="gNotebookRGtk"),
                 function(x, i, j,..., value) {
                   .leftBracket(x, x@toolkit, i, j, ...) <- value
                   return(x)
                 })

setReplaceMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gNotebookRGtk"),
          function(x, toolkit, i, j, ..., value) {
            n = length(x)
            if(missing(i)) {
              add(x,value)                        # append!
            } else {
              if(length(i) == 1) {
                add(x, value, index = i)
              } else {
                warning("Can't use '[' method for more than 1 element")
              }
            }
          })


### handlers
setMethod(".addhandlerchanged",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gNotebookRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            ## put page number into h$pageno
            widget <- getWidget(obj)
            ID <-
              gSignalConnect(widget,signal = "switch-page",
                             f = function(d,widget,page, pageno,...) {
                               h <- list(obj=d$obj,action=d$action, pageno=pageno+1)
                               if(!is.null(d$handler) &&
                                  is.function(d$handler))
                                 d$handler(h,...)
                               return(FALSE) # propogate
                             },
                             user.data.first = TRUE,
                             data = list(obj=obj,handler=handler, action=action)
                             )
            invisible(ID)
#            addhandler(obj,"switch-page",handler,action)
          })


setMethod(".addhandlerexpose",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gNotebookRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            addhandlerchanged(obj,handler, action,...)
          })

### helpers
## used in [. method
getNotebookPageWidget = function(obj, pageno =obj@widget$GetCurrentPage() + 1) {
  theLabel = obj@widget$GetTabLabel(obj@widget$GetNthPage(pageno))
  widget = tag(theLabel,"widget")
  return(widget)
}


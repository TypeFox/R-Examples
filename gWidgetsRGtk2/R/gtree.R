setClass("gTreeRGtk",
         contains="gComponentRGtk",
         prototype=prototype(new("gComponentRGtk"))
         )

## offspring takes two argument

##' toolkit constructor for gtree
setMethod(".gtree",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   offspring = NULL,
                   hasOffspring = NULL,                 # for defining offspring. FUN
                                        # of children only.
                   offspring.data = NULL,
                   col.types = NULL, # data frame with logical
                   icon.FUN = NULL,                      # return stock name --called
                                        # on offspring, returns a
                                        # vector of length nrow
                   chosencol = 1,
                   multiple = FALSE,
                   handler = NULL,
                   action=NULL,
                   container=NULL,
                   ...
                   ) {
            
            force(toolkit)

            ## do we have first col. for icons?
            iconFudge <- ifelse(is.null(icon.FUN), 0, 1)
            ## is second column for offspring?

            
            ## get base offspring
            children <- offspring(c(), offspring.data)

            ## we have some hacks here. First we place icon info into data frame if icon.FUN non NULL
            ## as well, we also strip out hasOffspring info into doExpand variable. This might be found
            ## from a function or from the second column -- if logical, or is just FALSE.
            
            ## we can have icons, if so we place in column 1
            ## column 2 can have offspring data!
            ## put in icons if needed
            lst <- getOffSpringIcons(children, hasOffspring, icon.FUN)
            children <- lst$children
            doExpand <- lst$doExpand
            
            ## ask before we put in icon info if asked
            if(is.null(col.types)) {
              col.types <- children[1,]
              if(iconFudge)
                col.types <- col.types[, -1] # shift out icon info
            }
            
            ## get GTK types -- force first to be character
            if(length(col.types) > 1) {
              types = c("gchararray", sapply(col.types[ ,-1],RtoGObjectConversion))
            } else {
              types <- "gchararray"
            }
            if(iconFudge == 1)
              types <- c("gchararray", types)       # stores filename of image

 

            
            
            ## define treestore, sorted, view
            treestore <- gtkTreeStoreNew(types)
            treestoreModel <- gtkTreeModelSortNewWithModel(treestore)
            view <- gtkTreeViewNewWithModel(treestoreModel)

            ##  if(nrow(children) > 15)
            ##    view$SetFixedHeightMode(TRUE)       # speeds up this. FAILED?
            view$SetSearchColumn(iconFudge)         # for CTRL-f
            
            ## define cellrender
            colHeaders <- names(children)
            
            for(i in (1+iconFudge):ncol(children)) {
              cellrenderer = gtkCellRendererTextNew()
              view.col = gtkTreeViewColumnNew()
              ## properties
              view.col$SetResizable(TRUE)
              ## title
              if(!is.na(colHeaders[i]) && !is.null(colHeaders[i]))
                view.col$SetTitle(colHeaders[i])
              view.col$SetSortColumnId(i-1) # allow sorting
              view.col$PackStart(cellrenderer, TRUE)
              view.col$AddAttribute(cellrenderer, "text", i-1)
              view$InsertColumn(view.col,i-1)
            }
            
            if(iconFudge == 1) {
              cellrenderer = gtkCellRendererPixbufNew()
              view.col = gtkTreeViewColumnNew()
              ## properties
#              view.col$SetMaxWidth(20) # 20 pixel icons
              view.col$PackStart(cellrenderer, TRUE)
              view.col$AddAttribute(cellrenderer, "stock-id", 0)
              view$InsertColumn(view.col,0)
            }  
            
            ## pack into scrolled window
            group = ggroup()
            sw <- gtkScrolledWindowNew()
            sw$SetPolicy("GTK_POLICY_AUTOMATIC","GTK_POLICY_AUTOMATIC")
            sw$Add(view)
            add(group, sw, expand=TRUE)
            
            ## allow multiple if asked
            if(multiple) {
              treeselection = view$GetSelection()
              treeselection$SetMode(GtkSelectionMode["multiple"])
            }
            
            
            ## turn on alternating shading if more than 1 column
            if(ncol(children) > 1)
              view$SetRulesHint(TRUE)
            

            obj = new("gTreeRGtk", block=group, widget=view, toolkit=toolkit)

            tag(obj,"store") <- treestore
            tag(obj,"SortedStore") <- treestoreModel
            tag(obj,"view") <- view
            tag(obj,"offspring") =offspring
            tag(obj,"hasOffspring") = hasOffspring
            tag(obj,"offspring.data") = offspring.data
            tag(obj,"icon.FUN") = icon.FUN
            tag(obj,"iconFudge") = iconFudge
            tag(obj,"chosencol") = chosencol
            tag(obj,"multiple") = multiple
            tag(obj,"ncols") = length(types)
          
            ## put in children, handler for expand-row
            addChildren(treestore, children, doExpand, iconFudge, parent.iter=NULL)
            
            ## now add a handler to row-exapnd
            addhandler(obj,"row-expanded",
                       handler = function(h,view, iter, path,...) {
                         ## get unsorted iter from path
                         uspath <- treestoreModel$ConvertPathToChildPath(path)
                         iter <- treestore$GetIter(uspath)$iter
                         path <- .getValuesFromIter(h$obj,iter)
                         children <- offspring(path,tag(obj, "offspring.data"))


                         lst <- getOffSpringIcons(children, hasOffspring, icon.FUN)
                         children <- lst$children
                         doExpand <- lst$doExpand
                         
                         addChildren(treestore, children, doExpand,
                                     tag(h$obj,"iconFudge"), iter)
                         ## remove errant offspring
                         child.iter <- treestore$IterChildren(iter)
                         if(child.iter$retval)
                           treestore$Remove(child.iter$iter)
                       })
            
            
            addhandler(obj,"row-collapsed",
                       handler = function(h, view, iter, path, ...) {

                         ## get unsorted iter from path
                         uspath = treestoreModel$ConvertPathToChildPath(path)
                         iter = treestore$GetIter(uspath)$iter

                         ## get children, remove
                         n = treestore$IterNChildren(iter)
                         if(n > 1) { ## n=1 gets removed when expanded
                           for(i in 1:(n-1)) {
                             child.iter = treestore$IterChildren(iter)
                             if(child.iter$retval)
                               treestore$Remove(child.iter$iter)
                           }
                         }
                       })
            
            if(!is.null(handler)) {
              id = addhandlerdoubleclick(obj,handler,action)
              tag(obj, "handler.id") <- id
            }
            
            ## attach to container
            if (!is.null(container)) {
              if(is.logical(container) && container == TRUE)
                container = gwindow(visible=TRUE)
              add(container, obj,...)
            }
            
            return(obj)
          })

## Take the data frame and massage it to return
## icons if asked, and figure out offspring
getOffSpringIcons = function(children, hasOffspring, icon.FUN) {

  ## do we expand?
  ## how to determine if offspring are needed?
  ## default to hasOffspring, then second column, then default to FALSE
  if(!is.null(hasOffspring)) {
    doExpand = hasOffspring(children)
  } else {
    ## if second column is logical, we use that
    if(is.logical(children[,2])) {
      doExpand = children[,2]
      children = children[,-2, drop=FALSE]
    } else {
      doExpand = rep(FALSE, nrow(children))
    }
  }
  
  ## make icons first column if there
  ## icon.FUN is called on data.frame, returns vector to cbind to children.
  if(!is.null(icon.FUN)) {
    if(nrow(children) > 0) {
      icons = getstockiconname(icon.FUN(children))
      children = data.frame(icons=I(icons), children)
    } else {
      children = data.frame( icons = character(0), children)
    }
  }
    
  return(list(children=children, doExpand=doExpand))
}


## children has label, logical, ...
## used to update tree
addChildren = function(treestore, children, doExpand, iconFudge, parent.iter=NULL) {
  if(nrow(children) == 0)
    return(NULL)

  ## load row by row, column by column
  for(i in 1:nrow(children)) {
    iter <- treestore$Append(parent=parent.iter)$iter
    ## no write values for each column
    for(j in 1:ncol(children)) {
      treestore$SetValue(iter,column=j-1, children[i,j])
    }
    ## add branch?
    if(!is.na(doExpand[i]) && doExpand[i]) {
      treestore$Append(parent=iter)
    }
  }
}

## has different arguments, but we mask this with ...
## this has offspringdata as first argument

setMethod("update",
          signature(object="gTreeRGtk"),
          function(object,...) {
            .update(object, object@toolkit, ...)
          })
setMethod(".update",
          signature(toolkit="guiWidgetsToolkitRGtk2",object="gTreeRGtk"),
          function(object, toolkit, ...) {
            
            theArgs = list(...)
            offspring.data = theArgs$offspring.data
            if(is.null(offspring.data) && length(theArgs))
              offspring.data = theArgs[[1]]
            if(!is.null(offspring.data))
              tag(object, "offspring.data") <- offspring.data
            

            obj = object                          # rename, object from update generic
            ## what should now be in this part of the tree
            newchildren <- tag(object,"offspring")(c(), tag(object, "offspring.data"))
            newvalues <- newchildren
            stillThere <- c()

            ## allow override by passing in function isStillThere into object via tag
            ## you may want to use get and digets here
            ## val is c(name, type) of item from tree;
            ## allVals is df with nameType of the newvalues to add to tree
            isStillThere <- function(val, allVals) {
              if(length(val) && length(allVals))
                val[1] %in% allVals[,1,drop=TRUE]
              else
                FALSE
            }
            isStillThere <- getWithDefault(tag(obj, "isStillThere"), isStillThere)


            
            ## loop over values in the treestore, if not in newchildren, remove
            i <- 0
            remove.these <- c()
            iter = tag(obj,"store")$GetIterFromString(i)
            while(iter$retval) {
              n <- ncol(newchildren) - is.null(tag(obj, "hasOffspring"))
              old <- sapply(1:n, function(i) {
                tag(obj,"store")$GetValue(iter$iter, i - 1 + tag(obj,"iconFudge"))$value
              })
#              treeValue <- tag(obj,"store")$GetValue(iter$iter,0  +tag(obj,"iconFudge"))$value
#              treeValueType <- tag(obj,"store")$GetValue(iter$iter,0+ 1 +tag(obj,"iconFudge"))$value

#              if(isStillThere(c(treeValue, treeValueType), newvalues)) {
              if(isStillThere(old, newvalues)) {
                stillThere <- c(stillThere, old[1])
              } else {
                ## need to delete
                remove.these = c(remove.these, i)
              }
              i = i + 1
              iter = tag(obj,"store")$GetIterFromString(i)
            }


            
            if(length(remove.these)>0) {
              for(i in rev(sort(remove.these))) {
                iter = tag(obj,"store")$GetIterFromString(i)
                tag(obj,"store")$Remove(iter$iter)
              }
            }
            
            didThese = newvalues[,1,drop=TRUE] %in% stillThere
            newchildren = newchildren[!didThese, , drop=FALSE] # don't drop dimension
            ## add these to end
            if(nrow(newchildren) > 0) {
              
              lst = getOffSpringIcons(newchildren, tag(obj,"hasOffspring"),
                tag(obj,"icon.FUN"))
              newchildren = lst$children
              doExpand = lst$doExpand
              ## add the children
              addChildren(tag(obj,"store"), newchildren, doExpand, tag(obj,"iconFudge"))
            }
          })

## XXX OLDuse index for the column to override the column returned
##' XXX Index should be for index of selected, e.g 1:2:3 type thing -- aka the path
##' @param index if TRUE, then return either a numeric vector or list of numeric vectors (if multiple selection)
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gTreeRGtk"),
          function(obj, toolkit, index=NULL, drop=NULL,...) {
            theArgs = list(...)

            index <- getWithDefault(index, FALSE)
            if(index) {
              ## make change -- return tree index
              twidget <- obj@widget
              sel <- twidget$getSelection()
              selectedRows <- sel$getSelectedRows()
              selList <- selectedRows$retval # list of GtkTreePaths
              if(length(selList) == 0) {
                ## no selection
                return(NULL)
              }
              out <- lapply(selList, function(i) {
                tmp <- i$toString()
                vals <- as.numeric(unlist(strsplit(tmp, ":"))) + 1
              })
              if(length(out) == 1)
                out <- out[[1]]         # return a list only if 2 or more
              
              return(out)
            } 
            
            ## we had case for both multiple or not, but we can use the same code for each
            ##            if(tag(obj,"multiple")) {
            treeselection = obj@widget$GetSelection()
            out = treeselection$GetSelectedRows() # 2 parts, paths, model
            if(length(out$retval) == 0) {
              return(NULL)
            } else {
              model = out$model
              tmp = c()
              for(i in out$retval) {
                iter = model$GetIter(i)$iter
                value = model$GetValue(iter,
                  tag(obj, "chosencol") - 1 + tag(obj,"iconFudge"))$value
                tmp = c(tmp,value)
              }
              return(tmp)
            }
          ## } else {
          ##     ## single selection
          ##     iter = obj@widget$GetSelection()$GetSelected()
          ##     if(iter$retval) 
          ##       return(obj@widget$GetModel()$GetValue(iter$iter,whichCol-1 +
          ##                                              tag(obj,"iconFudge"))$value)
          ##     else
          ##       return(NULL)              # nothing selected
          ##   }
          })


##' svalue<-
##'
##' Set selection by index. A path looks like c(a,b,c) 1-based
##' @param value indices. Either a vector for single selection or list of vectors for multiple selection.
##' @param index must be TRUE
setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gTreeRGtk"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   index <- getWithDefault(index, TRUE)
                   if(!index) {
                     gwCat(gettext("Need to have index=TRUE (or NULL)"))
                     return(obj)
                   }
                   ## value is character vector of paths. 
                   tr <- getWidget(obj)
                   sel <- tr$getSelection()
                   sel$unselectAll()    # clear selection

                   if(is.atomic(value))
                     value <- list(value)

                   value <- lapply(value, function(i) i-1)
                   
                   lapply(value, function(tmp) {
                     for(j in 1:length(tmp)) {
                       tpath <- gtkTreePathNewFromString(paste(tmp[1:j], collapse=":"))
                       tr$expandRow(tpath, open.all=FALSE)
                     }
                     sel$selectPath(tpath) ## adds if selection is multiple
                   })
                   return(obj)
                 })

### need to figure this out
## return the path in values. i,j ignored
setMethod("[",
          signature(x="gTreeRGtk"),
          function(x, i, j, ..., drop=TRUE) {
            .leftBracket(x, guiToolkit("RGtk2"), i, j, ..., drop=drop)
          })
setMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gTreeRGtk"),
          function(x, toolkit, i, j, ..., drop=TRUE) {
            obj = x
            
            ## XXX We had different cases for multiple and not, but this isn't necessary
            ##            if(tag(obj,"multiple"))  {
            twidget <- obj@widget
            sel <- twidget$getSelection()
            selectedRows <- sel$getSelectedRows()
            selList <- selectedRows$retval # list of GtkTreePaths
            model <- selectedRows$model
            if(length(selList) == 0) {
              ## no selection
              return(NULL)
            }
            
            out <- lapply(selList, function(path) {
              string <- path$toString()
              indices <- unlist(strsplit(string,":"))
              thePath <- c()
              for(j in 1:length(indices)) {
                npath <- paste(indices[1:j],collapse=":")
                iter <- tag(obj,"SortedStore")$GetIterFromString(npath)
                thePath[j] <- tag(obj,"SortedStore")$GetValue(iter$iter,0+
                                                        tag(obj,"iconFudge"))$value
              }
              thePath
            })
            if(length(out) == 1)      # if only 1 selection return it, o/w give as list
                out <- out[[1]]
            return(out)
          ## } else {
          ##   sel <- obj@widget$GetSelection()$GetSelected()
          ##     if(!sel$retval) {
          ##       ## no selection
          ##       return(character(0))
          ##     }
          ##     iter <- sel$iter
          ##     ## need to convert to unsorted
          ##     iter = tag(obj,"SortedStore")$ConvertIterToChildIter(iter)$child.iter
              
          ##     string = tag(obj,"store")$GetPath(iter)$ToString()
          ##     indices = unlist(strsplit(string,":"))
          ##     thePath = c()
          ##     for(j in 1:length(indices)) {
          ##       path = paste(indices[1:j],collapse=":")
          ##       iter = tag(obj,"store")$GetIterFromString(path)
          ##       thePath[j] = tag(obj,"store")$GetValue(iter$iter,0+
          ##                tag(obj,"iconFudge"))$value
          ##     }
          ##     if(missing(i))
          ##       return(thePath)
          ##     else
          ##       return(thePath[i])
          ##   }
          })
          
          
          
### methods
## row-activated in gtable gives double click
setMethod(".addhandlerdoubleclick",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gTreeRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
             addhandler(obj, "row-activated",handler,action,...)
           })

setMethod(".addhandlerchanged",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gTreeRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            .addhandlerdoubleclick(obj, toolkit, handler, action, ...)
          })


## clicked is on selection
setMethod(".addhandlerclicked",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gTreeRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            widget <- getWidget(obj)
            widget <- widget$getSelection()
             addhandler(widget, "changed",handler,action,actualobj=obj,...)
           })

## used internally
.getValuesFromIter = function(obj, iter) {
  string = tag(obj,"store")$GetPath(iter)$ToString()
  indices = unlist(strsplit(string,":"))
  thePath = c()
  for(i in 1:length(indices)) {
    path = paste(indices[1:i],collapse=":")
    iter = tag(obj,"store")$GetIterFromString(path)
    ## need to fudge here if necessary
    thePath[i] = tag(obj,"store")$GetValue(iter$iter,0+tag(obj,"iconFudge"))$value
  }
  return(thePath)
}



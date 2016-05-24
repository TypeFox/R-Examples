### file for making gtable and gdf.
## TODO:
## * always have rownames, switch for showing
## * ggrid is called by gtable, gdf. So ggrid (not exported) can haeva argument such as doIcons, doRownNames, doFilter, doSort
## ad drin colors for rfg, rbg
## svalue.gtkTreeView and svalue gtkTreeViewColumn
## for column index=TREE returns col.no. o/w vector
## function for showing a vector or data frame

setClass("gGridRGtk",
         contains="gComponentRGtk",
         prototype=prototype(new("gComponentRGtk"))
         )


## constructor for selecting values from a data set -- not meant for editing
setMethod(".gtable",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   items,
                   multiple = FALSE,
                   chosencol = 1,                        # for drag and drop, value
                   icon.FUN = NULL,
                   filter.column = NULL,
                   filter.labels = NULL,
                   filter.FUN = NULL,   # two args gtable instance, filter.labels element
                   handler = NULL,
                   action = NULL,
                   container = NULL,
                   ...) {
            

            force(toolkit)
            
            obj = .ggrid(
              toolkit,
              items=items,
              multiple = multiple,
              chosencol = chosencol,
              editable = FALSE,
              icon.FUN = icon.FUN,
              filter.column = filter.column,
              filter.labels = filter.labels,
              filter.FUN = filter.FUN,
#              doSort = FALSE,           # makes visible work
              doSort = TRUE,           # makes visible work
              doRownames = FALSE,
              handler=handler,
              action=action,
              container = container,
              ...)

            tag(obj,"type") <- "gtable" ## should be a class -- ughh
            return(obj)
          })
          
## constructor for editing a data frame
setMethod(".gdf",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   items = NULL,
                   name = deparse(substitute(items)),
                   do.subset = FALSE,
                   container=NULL,...)  {

            force(toolkit)
            ## the colors
            theArgs = list(...)
            .colors = theArgs$colors
            if(is.null(colors))
              .colors = c(
                bg = "navajo white",fg = "black",
                rbg = "white smoke",rfg="red"
                )
            
            obj = NULL
            if(is.null(items)) {
              ## popup dialog to get first column, then ring back
              ## need to use gbasicdialog to make modal
              group = ggroup(horizontal=FALSE)
              tbl = glayout(); add(group, tbl, expand=TRUE)
              theName = gedit("X1")
              theType = gdroplist(c("numeric","character"))##,"factor"))
              theNoRows = gspinbutton(from=1,to=100,by=1,value=1)
              tbl[1,1] = glabel("First variable name:");tbl[1,2] = theName
              tbl[2,1] = glabel("Its type:");tbl[2,2] = theType
              tbl[3,1] = glabel("No. rows:");tbl[3,2] = theNoRows
              visible(tbl) <- TRUE
              gbasicdialog(title="Describe first variable",
                           widget=group,
                           handler = function(h,...) {
                             tmp = cbind(do.call(paste("as.",svalue(theType),sep=""),
                               list(rep(NA, length=svalue(theNoRows)))))
                             colnames(tmp)[1] = svalue(theName)
                             items <<- tmp
                           })
  }
            obj <- .ggrid(
                          toolkit,
                         items = items,
                         multiple=FALSE,
                         chosencol = 1,
                         editable=TRUE,
                         doFilter=FALSE,
                         doIcons=FALSE,
                         doSort = FALSE,
                         doRownames=TRUE,
                         doSubsetBy = do.subset,
                         handler=NULL,
                         action=NULL,
                         container=container,
                         colors=.colors,
                          diy = theArgs$diy, # suppress ke or popup
                         ...)
            tag(obj,"type") <- "gdf" ## should be a class -- ughh

            
##             ## add 3rd mouse handler for the view
##             lst = list()
##             lst$"Apply function to column"$handler = function(h,...) {
##               col.no = h$action
              
##               win = gwindow("Apply function to column",visible=TRUE)
##               group = ggroup(horizontal = FALSE, container=win)
##               glabel("<b>Apply function to column</b>", markup=TRUE, container=group)
##               tmpGroup = ggroup(container=group)
##               glabel("<b>function(x) = {</b>", markup=TRUE,container=tmpGroup)
##               addSpring(tmpGroup)
##               FUN = gtext(container=group)
##               tmpGroup = ggroup(container=group)
##               glabel("}", container=tmpGroup)
##               addSpring(tmpGroup)
##               buttonGroup = ggroup(container=group)
##               addSpring(buttonGroup)
##               gbutton("ok",container=buttonGroup,handler = function(h,...) {
##                 FUN = Paste("function(x) {",svalue(FUN),"}")
##                 f = eval(parse(text=FUN))
##                 theNewVals = f(obj[,col.no, drop=FALSE])
##                 obj[,col.no] = theNewVals
##                 dispose(win)
##               })
##               gbutton("cancel",container=buttonGroup, handler = function(h,...)
##                       dispose(win))
##             }
##             lst$"Sort by column (increasing)"$handler = function(h,...) {
##               col.no = h$action
##               newOrder = order(obj[,col.no], decreasing = FALSE)
##               obj[,] = obj[newOrder,]
##               ## signal? -- is killing R
##               ##      cr = view.col$GetCellRenderers()[[1]] 
##               ##      try(cr$SignalEmit("edited"), silent=TRUE) # notify
##             }
##             lst$"Sort by column (decreasing)"$handler = function(h,...) {
##               col.no = h$action
##               newOrder = order(obj[,col.no], decreasing = TRUE)
##               obj[,] = obj[newOrder,]
##               ## signal?
##               ##      cr = view.col$GetCellRenderers()[[1]] 
##               ##      try(cr$SignalEmit("edited"), silent=TRUE) # notify
##             }
##             ## can't easily do this, as obj[,] wants to keep the same types
## ##             lst$"Coerce column type"$handler = function(h,...) {
## ##               colNum = h$action
## ##               theData = obj[,colNum,drop=TRUE]
## ##               theClass = class(theData)
## ##               allClasses = c("numeric","integer","character","factor","logical")
## ##               win = gwindow("Coerce column data")
## ##               g = ggroup(horizontal=FALSE, cont=win)
## ##               add(g,glabel("Select the new column type"))
## ##               gdroplist(allClasses,cont=g,selected = which(theClass == allClasses),
## ##                         handler = function(h,...) {
## ##                           newClass = svalue(h$obj)
## ##                           theData = do.call(paste("as.",newClass,sep="",collapse=""),list(theData))
## ##                           df = obj[,]; df[,colNum] <- theData
## ##                           obj[,] <- df
## ##                           dispose(win)
## ##                         })
## ##               add(g, gbutton("close",handler = function(...) dispose(win)))
## ##             }
##             ## rename -- tedious. Was better when label was editable
##             lst$"Rename column"$handler = function(h,...) {
##               col.no = h$action
##               view.col = tag(obj,"view")$GetColumn(
##                 col.no-1+tag(obj,"doRownames") + tag(obj,"doIcons"))
              
##               win = gwindow("Change name", visible=TRUE)
##               group = ggroup(horizontal=FALSE, container=win)
##               ok.handler = function(h,...) {
##                 names(obj)[col.no] <- svalue(h$action)
##                 dispose(win)
##                 if(tag(obj,"doSubsetBy")) {
##                   subsetBy = tag(obj,"subsetBy")
##                   update(subsetBy)
##                 }
##                 return(FALSE)
##               }
##               newName = gedit(id(view.col),container=group)
##               addhandlerchanged(newName, handler=ok.handler, action=newName)
##               buttonGroup = ggroup(container=group);addSpring(buttonGroup)
##               add(buttonGroup,gbutton("ok", handler = ok.handler, action=newName))
##               add(buttonGroup,gbutton("cancel",handler=function(h,...) dispose(win)))
##               return(TRUE)
##             }
            

##             f = function(h, widget, event,...) {
##               if(event$GetButton() != 3) {
##                 return(FALSE)                     # propogate signal
##               } else {
##                 cursor = widget$GetCursor()                
##                 view.col = cursor[["focus_column"]]
##                 if(is.null(view.col)) {
##                   view.col = cursor[['focus.column']] # view.col is the column
##                 }
##                 column.number = tag(view.col,"column.number")
##                 if(is.null(column.number)) {
##                   cat("Select a cell first by clicking once\n")
##                   return()
##                 }
##                 column.number = column.number - 1 + tag(obj,"doRownames") + tag(obj,"doIcons")
                
                
                
##                 mb = gmenu(h$action, popup = TRUE, action=column.number) #  action argument?
##                 mb = tag(mb,"mb")                 # actual gtkwidget
##                 print(class(mb))
##                 gtkMenuPopupHack(mb,
##                              button = event$GetButton(),
##                              activate.time=event$GetTime(),
##                              func = NULL
##                              )
##                 return(TRUE)
                
##               }
##             }

##             ## This isn't working!
## ##            addhandler(tag(obj,"view"),signal = "button-press-event",
## ##                       handler=f, action=lst)
## ##            
            return(obj)
            
          })
  

## make this generic, its not part of gWidgets API
setGeneric(".ggrid",function(toolkit,
                   items,                                # items to show: vector, matrix or df
                   multiple = FALSE,                     # allow multiple selection
                   chosencol = 1,                        # for drag and drop, svalue
                   editable = FALSE,                     # T -> gDF, F gtable
                   icon.FUN = NULL,                      # make icons?
                   filter.column = NULL,                 # do we filter easily?
                   filter.labels = NULL,                 # if we filter harder
                   filter.FUN = NULL,   # two args gtable instance, filter.labels element
                   doIcons = ifelse(is.null(icon.FUN),FALSE, TRUE),
                   doFilter = FALSE,    # see belo
                   doSort = TRUE,
                   doRownames = FALSE,
                   doSubsetBy = FALSE,
                   handler = NULL,                       # double click handler
                   action = NULL,                        # passed to handler
                   container = NULL,                     # optional container
                   ...) standardGeneric(".ggrid"))

           
setMethod(".ggrid",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   items,                                # items to show: vector, matrix or df
                   multiple = FALSE,                     # allow multiple selection
                   chosencol = 1,                        # for drag and drop, svalue
                   editable = FALSE,                     # T -> gDF, F gtable
                   icon.FUN = NULL,                      # make icons?
                   filter.column = NULL,                 # do we filter easily?
                   filter.labels = NULL,                 # if we filter harder
                   filter.FUN = NULL,   # two args gtable instance, filter.labels element
                   doIcons = ifelse(is.null(icon.FUN),FALSE, TRUE),
                   doFilter = FALSE,
                   doSort = TRUE,
                   doRownames = FALSE,
                   doSubsetBy = FALSE,
                   handler = NULL,                       # double click handler
                   action = NULL,                        # passed to handler
                   container = NULL,                     # optional container
                   ...) {
            
            force(toolkit)
            
            theArgs = list(...)                   # for colors
            if(is.null(theArgs$colors)) {
              theColors = c(
                bg = "navajo white",fg  = "black",
                rbg = "white smoke",rfg ="red"
#                fg="black",bg="white",            # default for fg and bg
#                rfg="white",rbg="black"          # for rows
                )
            } else {
              theColors = theArgs$colors          # must have names
            } 
            
            ## define the object
            group = ggroup(horizontal=FALSE, container = container, ...)
            ## make widget=group, later set to view
            obj = new("gGridRGtk", block=group, widget=group, toolkit=toolkit)

            tag(obj,"chosencol") <- chosencol
            tag(obj,"filter.column") <- filter.column # 1:n based
            tag(obj,"theColors") <- theColors
            tag(obj,"do.it.yourself") <- theArgs$diy # c("suppress.key","suppress.popup") || NULL


            ## what are we doing?
            iconFudge = ifelse(as.logical(doIcons), 1, 0)
            tag(obj,"doIcons") <- doIcons
            tag(obj,"icon.FUN") <- icon.FUN
            tag(obj,"doRownames") <- doRownames

            ## sort offilter?
            if(doFilter || !is.null(filter.column) ||
               !is.null(filter.FUN) ) {
              doFilter <- TRUE
            }
            if(doFilter) doSort <- FALSE          # can't sort and filter
            tag(obj,"doSort") <- doSort
            tag(obj,"doFilter") <- doFilter
            tag(obj,"doSubsetBy") <- doSubsetBy
            
            
            
            items = hack.as.data.frame(items)
            if(class(items)[1] != "data.frame") {
              warning("The items can not be coerced into a data frame")
              return(NA)                            # error message?
            }
            m = nrow(items); n = ncol(items)
            
            
            itemsPadded = makePaddedDataFrame(
              obj, items,
              visible = rep(TRUE, length=n)
              )
            
            
            store = rGtkDataFrame(itemsPadded)
            tag(obj,"store") <- store
            ## figure out whether we filter or sort
            ## do we filter? edit? neither?
            tag(obj,"editable") <- editable
            if(tag(obj,"doFilter")) {
              filter.popup = gdroplist(c("")) # replace with values if defined
              ## we filter *if$ filter.column is set or if filter.FUN is non null
              if(!is.null(filter.column)) {
                ## we filter based on value in this column. Define filter.labels 
                filter.labels = c("",sort(unique(as.character(store[,3*(filter.column+1)]))))
                filter.FUN = function(obj, filter.by) {
                  if(length(filter.by) == 0 || filter.by == "") {
                    vals = rep(TRUE, dim(obj)[1])
                  } else {
                    vals = as.character(obj[,filter.column,drop=TRUE]) == as.character(filter.by)
                  }
                  return(vals)
                }
                filterGroup = ggroup(container = group)
                glabel("Filter by:", container=filterGroup)
                filter.popup = gdroplist(filter.labels, container=filterGroup)
              } else {
                if(is.function(filter.FUN)) {
                  filterGroup = ggroup(container = group)
                  glabel("Filter by:", container=filterGroup)
                  filter.popup = gdroplist(filter.labels, container=filterGroup)
                }
              }
              tag(obj,"filter.FUN") <- filter.FUN

              ## if filter.FUN = is non null and *not* a function, no sorting
              if(is.function(filter.FUN)) {
                tag(obj,"filter.popup") <- filter.popup
                addhandlerchanged(filter.popup, action=obj,handler = function(h,...) {
                  vals = tag(h$action,"filter.FUN")(h$action, svalue(h$obj))
                  visible(h$action) <- vals
                })
              }
            }


            if(doSort) {
              ##model = gtkTreeModelSort(store)
              model = store
            } else {
              model = store$FilterNew()
              model$SetVisibleColumn(0)
            }
            
            view <- gtkTreeViewNew(TRUE)
            view$SetModel(model)
            
            tag(obj,"view") <- view
            tag(view, "gridObj") <- obj # no toolkit inside view
            obj@widget <- view          # replace widget
            
            ## add scroll window for tree view
            sw <- gtkScrolledWindowNew()
            sw$SetPolicy("GTK_POLICY_AUTOMATIC","GTK_POLICY_AUTOMATIC")
            sw$Add(view)
            add(group,sw, expand=TRUE, fill="both")

            ## properties
            if(multiple) {
              treeselection = view$GetSelection()
              treeselection$SetMode(GtkSelectionMode["multiple"])
            }
            ## turn on alternating shading if more than 1 column
            if(ncol(items) > 1)
              view$SetRulesHint(TRUE)
            ## search
            view$SetEnableSearch(TRUE)
            if(doRownames)
              view$SetSearchColumn(2)
            else
              view$SetSearchColumn(3*(chosencol+1)-1) # -1 to put in GTK 0 base

            ## Now to display the data 
            if(!is.null(icon.FUN))
              addIcons(view)
            for(j in 1:n) {
              if(tag(obj,"editable")) {
                view.col = addTreeViewColumnWithEdit(obj,j, colnames(items)[j])
              } else {
                view.col = addTreeViewColumnNoEdit(obj, j, colnames(items)[j])
              }
            }
            
            if(tag(obj,"doRownames")) {
              view.col = addTreeViewColumnWithEdit(obj, 0,"Row.names")
            }
            
            ## do we add subsetBy
            if(doSubsetBy) {
              ## now add subset by to group.cycling
              subsetByGroup = gexpandgroup("subset=", container=group)
              subsetBy = gsubsetby(obj, container=subsetByGroup,
                handler = function(h,...) {
                  visValues = h$value
                  return()
                  if(is.na(visValues[1]))
                    visible(h$action) <- rep(TRUE, dim(obj)[2])
                  else
                    visible(h$action) <- visValues
                })
              tag(obj,"subsetBy")<-subsetBy

              ## add subsetby to each view.col for a variable
              for(i in view$GetColumns()) tag(i,"subsetBy") <- subsetBy
#              sapply(view$GetColumns(), function(view.col)
#                     tag(view.col,"subsetBy") <- subsetBy)
              
            }

            if(tag(obj,"editable") &&
               ((is.null(tag(obj,"do.it.yourself")) ||
                !("suppress.key" %in% tag(obj,"do.it.yourself"))))
               ) {
              ## handler for moving around
              addKeyMotionHandler(obj)
            }
            

            ## add handler for double click
            if(!is.null(handler)) {
              id = addhandlerdoubleclick(obj, handler, action)
            }
            
            return(obj)
          })



### methods
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="GtkTreeView"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {
            svalue(tag(obj,"gridObj"),index, drop, ...)
          })
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gGridRGtk"),
          function(obj, toolkit, index=NULL, drop=NULL,...) {
            theArgs = list(...)
            view = tag(obj,"view")
            indices = .getSelectedIndices(obj,view)

            ## careful, we may be filtering!
            vals = visible(obj)
            
            if(!is.null(index) && index == TRUE) {
              ## we may be filtering. Need to undo
              return(which(vals)[indices])
            }

            ## Now a value. Works if filtering
            df = obj[vals,,drop=FALSE]
            
            if(!is.null(drop) && drop == FALSE)
              df[indices,,drop=FALSE]
            else
              df[indices, tag(obj,"chosencol"), drop=TRUE] # no drop=FALSE here
          })
          
## return indices for the original store, not filtered or sorted
.getSelectedIndices = function(obj, view, ...) {
  selection <- view$GetSelection()$GetSelectedRows()$retval

  if(length(selection) == 0)
    return(NULL)

  store <- view$GetModel()

  if(is.null(tag(obj,"type")) || tag(obj,"type") == "gdf") {
  
    indices <- sapply(selection,function(i) {
      ind <- store$ConvertPathToChildPath(i)$ToString()
      as.numeric(ind) + 1                 # shift to 1:m base
    })
  } else {
    indices <- sapply(selection,function(i) {
      ind <- as.numeric(i$ToString())
      as.numeric(ind) + 1                 # shift to 1:m base
    })
  }
            
  return(indices)
  
}
          
## set by index value selected value
setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gGridRGtk"),
                 function(obj, toolkit, index=NULL, ..., value) {

                   view <- tag(obj,"view")
                   selection <- view$GetSelection()

                   ## get indices, then select
                   if((!is.null(index) &&index == TRUE) ||
                      (is.null(index) && is.numeric(value))) {
                     ind <- as.integer(value) - 1L
                     
                   } else {
                     ## set by value -- not by index
                     curVals <- obj[,tag(obj,"chosencol")]
                     ind <- match(value, curVals)
                     if(length(ind) == 1 && is.na(ind)) {
                       ## exit if no match
                       ind <- NULL
                       selection$unselectAll()
                       return(obj)
                     }
                     ind <- ind -1L
                   }
                   ind <- ind[ind >= 0] # only non-negative indices
##                    if((!is.null(index) && index == FALSE) || !is.integer(value)) {
##                      ## get indices, then select
##                      (is.null(index) && is.integer(value))) {
##                      ind <-  as.character(as.integer(value) - 1L)
##                    } else {
##                      ind <-  as.character(as.integer(value) - 1L)
##                    }

                   
                   ## block handlers to quiet down change signal
                   if(length(ind)) {
                     blockhandler(selection)
                     selection$unselectAll()
                     unblockhandler(selection)
                   } else {
                     ## we want to call handler when 0 or negative index
                     selection$unselectAll()
                   }

                   lapply(ind, function(i) {
                     path <- gtkTreePathNewFromString(i)
                     selection$SelectPath(path)
                   })

                   ## move to cell unless none selected
                   if(is.null(ind) || length(ind) == 0 || (length(ind) ==1 && is.na(ind))) 
                     return(obj)
                   
                   i <- min(ind)
                   path <- gtkTreePathNewFromString(i)
                   view$scrollToCell(path)
                   
                   return(obj)
                 })

## helper function here
## unlike make.names this doesn't put "X" prefix
make.row.names <- function(x, unique=TRUE) {
  dups = duplicated(x)
  if(any(dups))
    x[dups] <- make.names(x,unique=unique)[dups]
  return(x)
}


setMethod("[",
          signature(x="GtkTreeView"),
          function(x, i, j, ..., drop=TRUE) {
            .leftBracket(x, x@toolkit, i, j, ..., drop=drop)
          })
setMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="GtkTreeView"),
          function(x, toolkit, i, j, ..., drop=TRUE) {
            gwCat("DEBUG: call leftBracket on gtkTreeView: deprecate?\n")

            gridObj = tag(x,"gridObj")
            if(missing(i) && missing(j))
              tmp = gridObj[,,...,drop=drop]
            else if(missing(i))
              tmp = gridObj[,j,...,drop=drop]
            if(missing(j))
              tmp = gridObj[i,,...,drop=drop]
            else
              tmp = gridObj[i,j,...,drop=drop]
            
            return(tmp)
          })

## refers to the entire data frame
## index returned by svalue(index=T) works here
setMethod("[",
          signature(x="gGridRGtk"),
          function(x, i, j, ..., drop=TRUE) {
            .leftBracket(x, x@toolkit, i, j, ..., drop=drop)
          })
setMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gGridRGtk"),
          function(x, toolkit, i, j, ..., drop=TRUE) {
            theArgs = list(...)         
            ## look for visible=TRUE to show only visible items
            showVisible = ifelse(is.null(theArgs$visible),FALSE,theArgs$visible)
            
            ## can't increase size of data frame *or* change the class
            store <- .getRGtkDataFrame(x)
            n <- (dim(store)[2] - 2)/3 -1
            
            frame <- store[ , 3*((1:n)+1), drop=FALSE]
            
            rownames(frame) <- make.row.names(store[,3])
            names(frame) <- names(x)

            ## handle missing values
            if(missing(i) && missing(j)) {
              i = if(showVisible) which(visible(x)) else seq_len(nrow(x))
              j = seq_len(n)
            } else if (missing(i)) {
              i = if(showVisible) which(visible(x)) else seq_len(nrow(x))
            } else if (missing(j)) {
              j = seq_len(n)
            }
            if(showVisible) 
              i = intersect(i,which(visible(x)))

            ## return
            return(frame[i,j,drop=drop])
          })

## [<-
setReplaceMethod("[",
                 signature(x="GtkTreeView"),
                 function(x, i, j,..., value) {
                   .leftBracket(x, x@toolkit, i, j, ...) <- value
                   return(x)
                 })
setReplaceMethod(".leftBracket",
                 signature(toolkit="guiWidgetsToolkitRGtk2",x="GtkTreeView"),
                 function(x, toolkit, i, j, ..., value) {
                   gridObj = tag(x,"gridObj")
                   if(missing(i) && missing(j))
                     gridObj[,,...] <- value
                   if(missing(i))
                     gridObj[,j,...] <- value
                   if(missing(j))
                     gridObj[i,,...] <- value
                   else
                     gridObj[i,j,...] <- value
                   return(x)
                 })

## Refers to the entire data frame, unsorted.
## This is kind of a brutal hack. I'd like to add in
## a way to colorize the NA values that are added, but gave up
setReplaceMethod("[",
                 signature(x="gGridRGtk"),
                 function(x, i, j,..., value) {
                   .leftBracket(x, x@toolkit, i, j, ...) <- value
                   return(x)
                 })

setReplaceMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gGridRGtk"),
          function(x, toolkit, i, j, ..., value) {
            d = dim(x); m= d[1]; n=d[2]
            dv = dim(as.data.frame(value))
            theColors = tag(x,"theColors")
            store = .getRGtkDataFrame(x)
            view = tag(x,"view")
            
            frame = as.data.frame(store)
            

            ## we have to be careful if we are *replacing*. If the size isn't
            ## the same, then we make a new store.
            if(missing(i)  && missing(j)) {
              if(dv[2] == 0) dv[2] <- 1
              if(dv[2] < n) {
                warning("Can't replace with fewer columns")
                return(x)
              } else if(dv[2] == n) {
                ## same number of columns. Now rows?
                ## what size is dv[1]?
                if(dv[1] == m) {
                  if(m > 0) {
                    ## straight replace -- same no. cols, rows
                    store[,3*((1:n)+1)] <- value
                    if(!is.null(rownames(value)))
                      store[,3] <- rownames(value)
                    return(x)
                  }
                } else {
                  ## fewer or more rows
                  ## make a new padded rGtkDataFrame, then replace model
                  if(dv[1] == 0) {
                    ## zero out
                    frame = frame[0,,drop=FALSE]
                  } else if(dv[1] < m) {
                    ## fewer rows, same columns
                    frame = frame[1:dv[1],,drop=FALSE]
                    if(dv[1] > 0)       # something to replace?
                      frame[,3*((1:n)+1)] <- value
                    ## leave row and column names out of this
                    ## user can replace with dimnames
                  } else {
                    ## more rows, same columns
                    ## need to lengthen data frame
                    ## strategy -- replace first rows, then add one at atime
                    ## for k=1 case
                    value = as.data.frame(value,stringsAsFactors=FALSE)

                    frame[1:m, 3*((1:n)+1)] <- value[1:m,]
                    
                    for(i in (m+1):dv[1]) {
                      replaceList = list(TRUE,"",i,frame[1,4],frame[1,5])
                      for(k in 1:n) {
                        replaceList[[3*(k+1)]] <- value[i,k] # value
                        replaceList[[3*(k+1)+1]] <- frame[1,3*(k+1)+1] #fg
                        replaceList[[3*(k+1)+2]] <- frame[1,3*(k+1)+2] #bg
                      }
                      frame[i,] <- replaceList
                    }
                  }
                  ## now swap out model in tree view
                  newstore = rGtkDataFrame(frame)
                  if(tag(x,"doSort")) {
                    ##model = gtkTreeModelSort(newstore)
                    model = newstore
                  } else {
                    model = newstore$FilterNew()
                    model$SetVisibleColumn(0)
                  }
                  view$SetModel(model)
                }
              } else if(dv[2] > n) {
                ## more columns, need to extend.
                ## first get right number of rows
                ## add /replace rows
                ## then add columns
                if(dv[1] <= m) {
                  ## fewer rows, truncate
                  frame = frame[1:dv[1],] 
                } else {
                  ## more rows and more columns, first add rows
                  ## lengthen rows
                  for(i in (m+1):dv[1]) {
                    newRowName = rownames(value)[i]
                    newRowName = make.row.names(c(newRowName,rownames(frame)), unique=TRUE)
                    replaceList = list(TRUE,"",newRowName,frame[1,4],frame[1,5])
                    for(j in 1:n) {
                      replaceList[[3*(j+1)]] <- value[i,j] # value
                      replaceList[[3*(j+1)+1]] <- frame[1,3*(j+1)+1] #fg
                      replaceList[[3*(j+1)+2]] <- frame[1,3*(j+1)+2] #bg
                    }
                    frame[i,] <- replaceList
                  }
                }
                ## finished with rows,
                ## now we need to add columns. We do so one column at a time
                for(j in (n+1):dv[2]) {
                  newVals = value[,j,drop=TRUE]
                  newPart = data.frame(
                    a=newVals,
                    b=rep(theColors['fg'],length=dv[1]),
                    c=rep(theColors['bg'],length=dv[1]),
                    stringsAsFactors = FALSE
                    )
                  names(newPart)[1] <- colnames(value)[j]
                  frame[,(3*(j+1)):(3*(j+1)+2)] = newPart
                }
                                        #      frame = adjustNA(value, frame)
                ## now swap out frame
                newstore = rGtkDataFrame(frame)
                if(tag(x,"doFilter")) {
                  model = newstore$FilterNew()
                  model$SetVisibleColumn(0)
                } else {
                  model = gtkTreeModelSort(newstore)
                }
                view$SetModel(model)
                ## now extend view -- 
                for(j in (n+1):dv[2]) {
                  if(tag(x,"editable"))
                    view.col = addTreeViewColumnWithEdit(x,j, colnames(value)[j])
                  else
                    view.col = addTreeViewColumnNoEdit(x, j,  colnames(value)[j])
                }
              }
              ## fix up the rownames
              store = .getRGtkDataFrame(x)
              if(dv[2] > 1 && dv[1] > 0)
                store[,3] <- rownames(value)
            } else {
              if(missing(i)) {
                ## no j is missing, just i
                if(dv[1] != m) {
                  warning("Sorry, you can't shorten or lengthen number of rows
        without replacing data frame. Try x[,] <- value instead.")
                  return(x)
                } else {
                  i = 1:m
                }
              } else if(missing(j)) {
                if(length(value) < n) {
                  warning("can't shorten number of columns")
                  return(x)
                } else if(length(value) > n) {
                  warning("To lengthen data frame, you must replace: x[,]<-value")
                  return(x)
                } else {
                  j = 1:n
                }
              }
              ## now we can assign with i and j
              store[i, 3*(j+1) ] <- value
            }
            
            ## fix icons if there`
             if(tag(x,"doIcons")) {
               store = .getRGtkDataFrame(x)
               d <- dim(store)
               if(d[1] > 0) {
                 n = (d[2]-2)/3 - 1
                 frame = store[,3*((1:n)+1)]
                 store[,2] = getstockiconname(tag(x,"icon.FUN")(frame))
               }
             }

            ## update filter
            if(tag(x,"doFilter") && !is.null(tag(x,"filter.column"))) {
              popup = tag(x,"filter.popup")
              vals = frame[,3*(tag(x,"filter.column")+1), drop=TRUE]
              popup[] <- c("",as.character(sort(unique(vals))))
              svalue(popup, index=TRUE) <- 1
            }
            
            
            return(x)
          })
                 
## first column is the visible row

setMethod(".visible",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gGridRGtk"),
          function(obj,toolkit,set=NULL, ...) {
            ## visible is only of value if sorting is not taking place
            if(tag(obj, "doSort"))
              return(rep(TRUE, dim(obj)[1]))

            ## not sorting, so first column holds visibility info
            frame = .getRGtkDataFrame(obj)
            return(frame[,1, drop=TRUE])
          })

## sets the first column
setReplaceMethod(".visible",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gGridRGtk"),
                 function(obj, toolkit, ..., value) {
                   if(tag(obj, "doSort")) {
                     #gwcat(gettext("Can't use visible<- method unless filtering is being used\n"))
                     return(obj)        # no means to set
                   }
                   
                   frame = .getRGtkDataFrame(obj)
                   m = nrow(frame)
                   frame[,1] <- rep(value, length=m)
                   
                   try(tag(obj,"view")$GetModel()$Refilter(), silent=TRUE)            # show
                   return(obj)
                 })


## data frame like
setMethod(".dim", 
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gGridRGtk"),
          function(x,toolkit) {
            store = .getRGtkDataFrame(x)
            tmp = dim(store)
            return(c(tmp[1], (tmp[2]-2)/3 - 1))
          })

setMethod(".dimnames",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gGridRGtk"),
          function(x,toolkit) {
            store = .getRGtkDataFrame(x)
            rownames = make.row.names(store[,3])
            colnames = names(x)
            return(list(rownames, colnames))
          })

setReplaceMethod(".dimnames",
                 signature(toolkit="guiWidgetsToolkitRGtk2",x="gGridRGtk"),
                 function(x, toolkit,  value) {
                   if(!is.list(value))
                     stop("value is a list with first element the row names, and second the column names")
                   rnames = value[[1]]
                   cnames = value[[2]]
                   d = dim(x)
                   if(is.null(rnames) || length(rnames) != d[1])
                     stop("Row names are the wrong size")
                   if(is.null(cnames) || length(cnames) != d[2])
                     stop("Column names are the wrong size")
                   
                   ## set column names
                   names(x) <- cnames
                   ## set row names
                   store = tag(x,"store")
                   if(is.null(store))
                     store = tag(x, "store")
                   if(length(rnames) > 0)
                     store[,3] <- make.row.names(rnames)
                   
                   return(x)
                 })

setMethod(".length",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gGridRGtk"),
          function(x,toolkit) return(dim(x)[2]))


setMethod(".names",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gGridRGtk"),
          function(x, toolkit) {
            view.cols = tag(x,"view")$GetColumns()
            theNames = character(length(view.cols))
            for(i in 1:length(view.cols))
              theNames[i] = id(view.cols[[i]])
            
            if(tag(x,"doRownames"))
              theNames = theNames[-1]
            if(tag(x,"doIcons"))
              theNames = theNames[-1]
            return(theNames)
          })

setReplaceMethod(".size", 
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gGridRGtk"),
                 function(obj, toolkit, ..., value) {

                   w <- getWidget(obj)

                   if(is.list(value) && !is.null(value$columnWidths)) {
                     colWidths <- value$columnWidths
                     colWidths <- rep(colWidths, length.out=dim(obj)[2])
                     sapply(seq_len(colWidths), function(i) {
                       col <- w$getColumn(i-1)
                       col$setMinWidth(colWidths[i])
                     })
                   }

                   if(is.list(value) && !is.null(value$rowHeights)){
                     ## no height method
                   }

                   ## width/height now
                   if(is.list(value)) {
                     width <- value$width # possibly NULL
                     height <- value$height
                   } else {
                     width <- value[1]
                     height <- ifelse(length(value) > 1, value[2], -1)
                   }
                   if(is.null(width))
                     width <- -1
                   if(is.null(height))
                     height <- -1

                   if(!is.null(width))
                     w$SetSizeRequest(width,height)

                  return(obj)
                })

setReplaceMethod(".names",
                 signature(toolkit="guiWidgetsToolkitRGtk2",x="gGridRGtk"),
                 function(x, toolkit, value) {
                   ## check that dimensions are correct
                   n = length(x)
                   if (length(value) != n)
                     stop("vector of names is the wrong length")
                   ## fix up names as needed
                   value = make.names(value, unique=TRUE)
                   view.cols = tag(x,"view")$GetColumns()
                   ## which cols depends on rownames
                   sapply(1:n, function(i) { # adding logical
                     id(view.cols[[i  + tag(x,"doRownames") + tag(x,"doIcons") ]]) <- value[i] 
                   })
                   
                   return(x)
                 })

## method to place edit
setCursorAtCell = function(obj, i,j,start.editing=TRUE) {
  if(is(obj,"gGridRGtk") || is(obj,"guiWidget"))
    view = tag(obj,"view")
  else
    view = obj

      
  path = gtkTreePathNewFromString(i-1)  # offset
  view.col = view$GetColumn(j-1 + tag(obj,"doRownames")+tag(obj,"doIcons"))        # offset
  view$SetCursor(path=path, focus.column=view.col,start.editing=start.editing)
  return(TRUE)
}


## Functions  for colors -- were methods, but didn't use

fgcolors = function(obj,i,j, ...) {
  d = dim(obj)
  ## return grid of colors
  if(missing(i)) i = 1:d[1]
  if(missing(j)) j = 1:d[2]

  frame = as.data.frame(.getRGtkDataFrame(obj))
  frame[i, 3*(j+1)+1, drop=FALSE]
}

"fgcolors<-" = function(obj, i,j, ..., value) {
  ## if both missing, assume value repeats down columns
  d = dim(obj)
  store = .getRGtkDataFrame(obj)
  
  if(missing(i) && missing(j)) {
    value = rep(value, length=d[2])     #recycle
    sapply(1:d[2], function(k)
           store[1:d[1],3*(k+1)+1] <- rep(value[k],length=d[1]))
  } else if(missing(i)) {
    ## only assign to j values
    rep(value, length=length(j))
    sapply(1:length(j), function(k)
           store[1:d[1],3*(j[k]+1)+1] <- rep(value[k],length=d[1]))
  } else if(missing(j)) {
    ## assume it runs down rows
    rep(value, length=length(i))
    sapply(1:length(i), function(k)
           store[i[k],3*((1:d[2])+1)+1] <- rep(value[k],length=d[2]))
   } else {
     store[i,3*(j+1)+1] <- value
   }
  return(obj)
}

bgcolors = function(obj, i,j,...) {
  d = dim(obj)
  ## return grid of colors
  if(missing(i)) i = 1:d[1]
  if(missing(j)) j = 1:d[2]

  frame = as.data.frame(.getRGtkDataFrame(obj))
  frame[i, 3*(j+1)+2, drop=FALSE]
}

"bgcolors<-" = function(obj, i,j, ..., value) {
  ## if both missing, assume value repeats down columns
  d = dim(obj)
  store = .getRGtkDataFrame(obj)
  
  if(missing(i) && missing(j)) {
    value = rep(value, length=d[2])     #recycle
    sapply(1:d[2], function(k)
           store[1:d[1],3*(k+1)+2] <- rep(value[k],length=d[1]))
  } else if(missing(i)) {
    ## only assign to j values
    rep(value, length=length(j))
    sapply(1:length(j), function(k)
           store[1:d[1],3*(j[k]+1)+2] <- rep(value[k],length=d[1]))
  } else if(missing(j)) {
    ## assume it runs down rows
    rep(value, length=length(i))
    sapply(1:length(i), function(k)
           store[i[k],3*((1:d[2]) +1)+2] <- rep(value[k],length=d[2]))
   } else {
     store[i,3*j+2] <- value
   }
  return(obj)
}



## handlers
setMethod(".addhandlerdoubleclick",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gGridRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            ## need to put onto view -- not group
##            id = addhandler(tag(obj,"view"), "row-activated",handler,action)
            id = addhandler(obj, "row-activated",handler,action,...)
            invisible(id)
          })

## gdf: click on headers -- passed on to each treeview
## gtable: click on row
##         note: to block gtable, need blockHandler(getToolkitWidget(b)$getSelection(), id)
setMethod(".addhandlerclicked",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gGridRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            if(tag(obj,"type") == "gdf") { # hack, should have a class here
              sapply(tag(obj,"view")$GetColumns(), function(object) {
                addhandlerclicked(tag(object,"widget"), handler, action,...)
              })
            } else {
              ## gtable -- put onto selection
              sel <-  obj@widget$GetSelection()
              ID <- gtktry(connectSignal(sel,
                signal = "changed",
                f = function(h,...) {
                  h$handler(h,...)
                },
                data = list(obj=obj, action=action, handler=handler),
                user.data.first = TRUE,
                after = FALSE
                ),
                silent=TRUE)

              ## add to selection
              l <- tag(sel,"handler.id")
              if(is.null(l))
                l <- list()
              l <- c(l, ID)
              tag(sel,"handler.id", replace=TRUE) <- l
              
              invisible(ID)
            }
          })

## for gdf -- change value
## for gtable -- double click
setMethod(".addhandlerchanged",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gGridRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            if(tag(obj,"type") == "gdf") {
              ## apply handler to change of each treeviewcolumn
              if(!missing(handler)) {     # only if handler is not missing
                view = tag(obj,"view")
                for(i in view$GetColumns())
                  addhandlerchanged(i, handler, action,...)
              }
            } else {
              ## gtable -- double click
              addhandler(obj, "row-activated",handler,action,...)
            }
          })


## Header handlers -- column clicked and columnrightclick
setMethod(".addhandlercolumnclicked",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gGridRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            ## apply handler to change of each treeviewcolumn
            if(!missing(handler)) {     # only if handler is not missing
              view = tag(obj,"view")
              viewCols <- view$getColumns()
              IDs <- sapply(seq_along(viewCols), function(i) {
                vc <- viewCols[[i]]
                widget <- vc$getWidget()
                widget <- widget$getParent()$getParent()$getParent()
                addhandlerclicked(widget, handler, action,
                                  column=i,
                                  ...)
              })
              invisible(IDs)
            }
          })
setMethod(".addhandlercolumnrightclick",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gGridRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            ## apply handler to change of each treeviewcolumn
            if(!missing(handler)) {     # only if handler is not missing
              view = tag(obj,"view")
              viewCols <- view$getColumns()
              IDs <- sapply(seq_along(viewCols), function(i) {
                vc <- viewCols[[i]]
                widget <- vc$getWidget()
                widget <- widget$getParent()$getParent()$getParent()
                addhandlerrightclick(widget, handler, action,
                                     column=i,
                                     ...)
              })
              invisible(IDs)
            }
          })

setMethod(".addhandlercolumndoubleclick",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gGridRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            ## apply handler to change of each treeviewcolumn
            if(!missing(handler)) {     # only if handler is not missing
              view = tag(obj,"view")
              viewCols <- view$getColumns()
              IDs <- sapply(seq_along(viewCols), function(i) {
                vc <- viewCols[[i]]
                widget <- vc$getWidget()
                widget <- widget$getParent()$getParent()$getParent()
                addhandlerdoubleclick(widget, handler, action,
                                      column=1,
                                      ...)
              })
              invisible(IDs)
              }
          })


### helpers
###############################
##################
## internal function
  ## return the rGtk data frame stored in obj
.getRGtkDataFrame = function(obj, ...) {
  view = tag(obj,"view")
  model <- view$getModel()
  if(!is(model, "RGtkDataFrame"))
    model <- model$getModel()

  return(model)
}


## XXX This is a mess -- replace me one day XXX
## the data frame has columns
## 1 -- visibilit
## 2 -- either rownames or icon names
## (3i, 3i+1, 3i+2) - df[,i],colors[1],colors[2]
makePaddedDataFrame <- function(obj,
                                items,
                                visible = rep(TRUE,length=ncol(items))
                                ) {
  
  
  
  m = nrow(items)
  n = ncol(items)
  cnames = colnames(items)
  reducedClass <- function(x) {
    out <- class(x)
    out[length(out)]
  }
  theClass = sapply(items, reducedClass)
  
  
  firstCol = rep(visible, length=m)
  
  ## we always add 5 cols at first: visible, icons, rownames, and rfg, rbg
  ## then for each col, we add 3 values, fg, bg
  if(is.null(tag(obj,"icon.FUN")))
    secondCol = rep("", length=m)
  else
    secondCol = getstockiconname(tag(obj,"icon.FUN")(items))
  
  ## go get em -- ugly code, how to make data frame *without* factor class?
  lst = list()
  lst[[1]] = firstCol;
  lst[[2]] = secondCol
  lst[[3]] = rownames(items)
  theColors = tag(obj,"theColors")
  lst[[4]] = rep(theColors['rfg'], length=m)
  lst[[5]] = rep(theColors['rbg'], length=m)

  bgColors = rep(theColors["bg"], length=m)
  fgColors = rep(theColors["fg"], length=m)
  sapply(1:n, function(j) lst[[3*(j+1)]] <<- items[,j])
  sapply(1:n, function(j) lst[[3*(j+1) + 1]] <<- fgColors) # foreground first
  sapply(1:n, function(j) lst[[3*(j+1) + 2 ]] <<- bgColors)

#  frame = do.call("data.frame", lst)
  frame <- data.frame(lst, stringsAsFactors=FALSE)
  
  ## coerce to proper class
  frame[,1] <- as.logical(frame[,1])       # visible
##   for(i in 2:5) frame[,i] = as.character(frame[,i])

##   for(j in seq_along(n)) {
##     if(theClass[j] != "AsIs")
##       frame[,3*(j+1)] = do.call(Paste("as.",theClass[j]),
##              list(x=frame[,3*(j+1)]))
##     frame[,3*(j+1) + 1] = as.character(frame[,(3*(j+1)+1)])
##     frame[,3*(j+1) + 2] = as.character(frame[,(3*(j+1)+2)])
##   }

  ## if we want to trap NA values, this works 
#  ## for NA values, we change colors
#  areNA = which(is.na(items), arr.ind=TRUE)
#  if(is.matrix(areNA)) {
#    lst = split(areNA, areNA[,2])
#    for(j in names(lst)) {
#      i = lst[[j]]; if(is.matrix(i)) i = i[,1]
#      frame[i,3*(as.numeric(j)+1) + 1] = theColors["nafg"]
#      frame[i,3*(as.numeric(j)+1) + 2] = theColors["nabg"]
#    }
#  }

  ## colnames
  colNames = rep("",3*n+2)
  colNames[1] = "visible"; colNames[2]="icons"
  colNames[3] = "Row.names"
  colNames[4] = "rfg"; colNames[5] = "rbg"
  colNames[3*((1:n)+1)] = cnames
  colNames[3*((1:n)+1) + 1]  = paste("fgCol",1:n, sep="")
  colNames[3*((1:n)+1) + 2]  = paste("bgCol",1:n, sep="")

  colnames(frame) <- colNames
  
  
  return(frame)
}

addIcons = function(view) {
  cellrenderer = gtkCellRendererPixbufNew()
  view.col = gtkTreeViewColumnNew()
  view.col$PackStart(cellrenderer, TRUE)
  view.col$AddAttribute(cellrenderer, "stock-id", 1)
  view$InsertColumn(view.col,0)
}
  
## j is in 1:n *or* 0 for rownames
addTreeViewColumnNoEdit <- function(obj, j,label) {
  view <- tag(obj, "view")
  
  cellrenderer <- gtkCellRendererTextNew()
  view.col <- gtkTreeViewColumnNew()
  view.col$PackStart(cellrenderer, TRUE)

  
  id(view.col) <- label
  ## store these
  tag(view.col,"column.number") <- j# add this for later usage
  tag(view.col,"view") <- view
  tag(view.col,"gridObj") <- obj
  
  ## properties
  view.col$SetResizable(TRUE)
  view.col$SetClickable(TRUE)
  if(tag(obj,"doSort")) {
    view.col$SetSortColumnId(3*(j+1) - 1)
  }
  ##

  ## TOO SLOW and does'
##   view.col$SetCellDataFunc(cellrenderer, func=
##                            function(vc, cr, model, iter, data) {
##                              curVal <- model$getValue(iter, data)$value
##                              cr['text'] <- if(is.na(curVal)) {
##                                "NA"
##                              } else if(is.nan(curVal)) {
##                                "NaN"
##                              } else if(is.null(curVal)) {
##                                "NULL"
##                              } else {
##                                curVal
##                              }
##                              cr['foreground'] <- model$getValue(iter, data + 1)$value
##                              cr['background'] <- model$getValue(iter, data + 2)$value
##                            },
##                            func.data = 3*(j+1) - 1 + tag(obj,"doIcons") + tag(obj,"doRownames"))


   view.col$AddAttribute(cellrenderer, "text", 3*(j+1) - 1)
   if(!is.null(tag(obj,"type")) && tag(obj,"type") == "gdf") {
     view.col$AddAttribute(cellrenderer,"foreground",3 *(j+1) + 1 - 1)
     view.col$AddAttribute(cellrenderer,"background",3 *(j+1) + 2 - 1)
   }
  
  view$InsertColumn(view.col,
                    j - 1 + tag(obj,"doIcons") + tag(obj,"doRownames"))
  
  return(view.col)
  
}

##################################################
### Define some key functions
### this is the main handler for editing data
## movement after editing is wanky! This pushes down a row, then the
## handler on the view starts editing. Unfortunately, it doesn't
## save then move
edit.handler = function(h,cell,path,newtext) {
  if(is.null(path) || is.null(newtext))
    return(FALSE)                     # propogate
  ## get position
  obj = h$action                        # the gGrid object
  column.number = h$column.number
  store = .getRGtkDataFrame(obj)
  
  i = as.numeric(path) + 1           # row
  ## if visible, we need to adjust i to point to right row
  if(colnames(store)[1] == "visible")
    i <- which(store[,"visible", drop=TRUE])[i]
  j = column.number

  
  ## coerce newtext from text to proper class
  theColData = obj[,j]

  if(is.integer(theColData)) {
    newtext = as.integer(newtext)
  } else if(is.numeric(theColData)) {
    newtext = as.numeric(newtext)
  } else if(is.character(theColData)) {
    newtext = as.character(newtext)
  } else if(is.factor(theColData)) {
    if(newtext %in% levels(theColData)) {
      ## nothing
    } else {
      levels(obj[,j]) = c(levels(theColData), newtext) ## was df
        ## tried a popup window, didn't work
    }
  } else if(is.logical(theColData)) {
    newtext = as.logical(newtext)
  } else {
    newtext = newtext               # nothing p
  }


  ## update foreground color == if was NA then fg=bg
  store[i,3*(j+1)] = newtext
  if(j == 0) {
    store[i,3*(j+1)+1] = tag(obj,"theColors")['rfg']  
    store[i,3*(j+1)+2] = tag(obj,"theColors")['rbg']
  } else {
    store[i,3*(j+1)+1] = tag(obj,"theColors")['fg']  
    store[i,3*(j+1)+2] = tag(obj,"theColors")['bg']
  }
  
  ## update subsetby if there
  doSubsetBy = tag(obj,"doSubsetBy")    # a logical or noull
  if(!is.null(doSubsetBy) && doSubsetBy) update(tag(obj,"subsetBy"))

  return(TRUE)
}


addTreeViewColumnWithEdit = function(obj, j,label) {
  view = tag(obj,"view")

  cellrenderers = tag(obj,"cellrenderers")
  if(is.null(cellrenderers)) {
    cellrenderers = list()
    tag(obj,"cellrenderers") <- cellrenderers
  }
  
  cellrenderer = gtkCellRendererTextNew()
  tag(obj,"cellrenderers") <- c(cellrenderers, cellrenderer) ## store to test
  
  ## properties
  gObjectSet(cellrenderer,"editable"=TRUE)
  gObjectSet(cellrenderer,"rise"=-10)

  view.col = gtkTreeViewColumn()
  ## add these for later usage
  tag(view.col,"column.number") <- j
  tag(view.col,"view") <- view
  tag(view.col,"gridObj") <- obj
  
  view.col$SetResizable(TRUE)
#  view.col$SetClickable(TRUE)

  ## cell renderers
  view.col$PackStart(cellrenderer, TRUE)
  id(view.col) <- label
  ## Need to fix this up, if numeric then editable is wanky
  view.col$AddAttribute(cellrenderer, "text", 3*(j+1) - 1)
  view.col$AddAttribute(cellrenderer,"foreground",3 *(j+1) + 1 - 1)
  view.col$AddAttribute(cellrenderer,"background",3 *(j+1) + 2 - 1)
  
  ## edit signal
  callbackId = gtktry(connectSignal(cellrenderer,
                signal = "edited",
                f=edit.handler,
                data = list(obj=cellrenderer,action=obj,column.number = j),
                user.data.first = TRUE,
                after=FALSE),
    silent=TRUE)
  ##
  view$InsertColumn(view.col,
                    j - 1 + tag(obj,"doIcons") + tag(obj,"doRownames"))
  

  ## fix up a bit
  addDragAndDropToViewCol(view.col)

  ## add popup unless asked not to. (YN suggestion)
  if(is.null(tag(obj,"do.it.yourself")) || !("suppress.popup" %in% tag(obj,"do.it.yourself")))
    addPopupMenuToViewCol(view.col)


  ## return the column
  return(view.col)


}


addPopupMenuToViewCol = function(view.col) {
  view.col$SetClickable(TRUE)          # make clickable headers
  lst = list()
  lst$"Apply function to column"$handler = function(h,...) {
    win = gwindow("Apply function to column",visible=TRUE)
    group = ggroup(horizontal = FALSE, container=win)
    glabel("<b>Apply function to column</b>", markup=TRUE, container=group)
    tmpGroup = ggroup(container=group)
    glabel("<b>function(x) = {</b>", markup=TRUE,container=tmpGroup)
    addSpring(tmpGroup)
    FUN = gtext(container=group)
    tmpGroup = ggroup(container=group)
    glabel("}", container=tmpGroup)
    addSpring(tmpGroup)
    buttonGroup = ggroup(container=group)
    addSpring(buttonGroup)
    gbutton("ok",container=buttonGroup,handler = function(h,...) {
      obj = tag(view.col,"gridObj")
      FUN = Paste("function(x) {",svalue(FUN),"}")
      f = eval(parse(text=FUN))
      col.no = tag(view.col,"column.number") ## not correct: - 1 # rownames offset

      oldVals = obj[,col.no, drop=FALSE]
      theNewVals = f(oldVals)
      obj[,col.no] = theNewVals
      dispose(win)
    })
    gbutton("cancel",container=buttonGroup, handler = function(h,...)
            dispose(win))
  }
##   lst$"Clear column"$handler = function(h,...) {
##     col.no = tag(view.col,"column.number") ## - 1 # rownames offset
##     obj[,col.no] = rep(NA, length(view.col))
##   }
  lst$"Sort by column (decreasing)"$handler = function(h,...) {
    col.no = tag(view.col,"column.number") ## - 1 # rownames offset
    newOrder = order(obj[,col.no], decreasing = TRUE)
    obj[,] = obj[newOrder,]
#    rownames(obj) = rownames(obj)[newOrder]
    ## signal?
    ##      cr = view.col$GetCellRenderers()[[1]] 
    ##      try(cr$SignalEmit("edited"), silent=TRUE) # notify
  }
  lst$"Sort by column (increasing)"$handler = function(h,...) {
    col.no = tag(view.col,"column.number") ## - 1 # rownames offset
    newOrder = order(obj[,col.no], decreasing = FALSE)
    obj[,] = obj[newOrder,]
#    rownames(obj) = rownames(obj)[newOrder]
    ## signal? -- is killing R
    ##      cr = view.col$GetCellRenderers()[[1]] 
    ##      try(cr$SignalEmit("edited"), silent=TRUE) # notify
  }
  lst$"Rename column"$handler = function(h,...) {
    win = gwindow("Change name", visible=TRUE)
    group = ggroup(horizontal=FALSE, container=win)
    ok.handler = function(h,...) {
      newVal = make.names(svalue(h$action))
      id(view.col) <- newVal
      ## signal
      ##        cr = view.col$GetCellRenderers()[[1]] 
      ##        try(cr$SignalEmit("edited"), silent=TRUE) # notify
      
      dispose(win)
      if(!is.null(tag(obj,"doSubsetBy")) &&
         tag(obj,"doSubsetBy")
         )                    #
        update(tag(obj,"subsetBy"))                # update
      
      return(FALSE)
    }
    newName = gedit(id(view.col),container=group)
    addhandlerchanged(newName, handler=ok.handler, action=newName)
    buttonGroup = ggroup(container=group);addSpring(buttonGroup)
    add(buttonGroup,gbutton("ok", handler = ok.handler, action=newName))
    add(buttonGroup,gbutton("cancel",handler=function(h,...) dispose(win)))
    return(TRUE)
  }
  
  ## define this so that it gets picked up in popup handlers
  obj = tag(view.col,"gridObj")

  ## put popup onto this guy -- button doesn't get 3rd mouse signal
  widget = tag(view.col,"widget")
  gtkbutton = view.col$GetWidget()$GetParent()$GetParent()$GetParent()
  ## was widget
#  add3rdmousepopupmenu(widget, menulist=lst)
  add3rdmousepopupmenu(gtkbutton, menulist=lst)  
#  add3rdmousepopupmenu(tag(view.col,"header"), menulist=lst)
}


## drag and drop handler on button oflabel
addDragAndDropToViewCol = function(view.col) {
  ## the widget is set by the "id<-" method
  ## this attaches an event box to the GetWidget()
#  gtkbutton = view.col$GetWidget()$GetParent()$GetParent()$GetParent()
  force(adddropsource(view.col$GetWidget()$GetParent()$GetParent()$GetParent(),
                      targetType="object",
                      action = view.col))
}

### This one moves the cursor and then sets the state to editing
addKeyMotionHandler = function(obj, ...) {
  view = tag(obj,"view")
  addhandler(view,"key-release-event",action=obj,handler=function(h,widget,event,...) {
    obj = h$action
    d = dim(obj)
    keyval = event$GetKeyval()
    cursor = widget$GetCursor()
    ## i,j are current positions,
    i = cursor$path$ToString()
    i = as.numeric(i) + 1               # in 1:m coordinates
    view.col = cursor[['focus_column']] # view.col is the column
    if(is.null(view.col)) {
      view.col = cursor[['focus.column']] # view.col is the column
    }
    j = tag(view.col,"column.number")
    ## where to move to
    if( keyval == GDK_Down ) {
      ## do we need to add a new row?
      ## for down arrow we can, for enter we don't
      if(i == d[1]) {
        frame = obj[,,drop=FALSE]
        ## pad with NA values
        lst = list();for(tmp in 1:d[2]) lst[[tmp]]=NA
        frame[d[1]+1,] = lst
        obj[,] = frame
      }
      ## move down
      setCursorAtCell(obj,i+1,j, start.editing=TRUE)
    } else if(keyval == GDK_Return) {
      if(i !=d[1]) i  = i+ 1         # dont add unless downarrow
      setCursorAtCell(obj, i, j, start.editing=TRUE)
    } else if( keyval == GDK_Up) {
      if(i  > 1)             # can't go too long
        i = i - 1
      setCursorAtCell(obj, i, j, start.editing=TRUE)
    } else if(keyval == GDK_Tab) {
      ## move to right
      ## add new column if at d[2] already
      ## Add dialog in case we are at last column
      if(j == d[2]) {
        setCursor = FALSE
        addNewColumnDialog(obj,i, j)
      } else  {
        setCursorAtCell(obj, i, j + 1, start.editing=TRUE)
      }
    }
    return(TRUE)
  })
}

## Dialog to add a new column
addNewColumnDialog <- function(obj, i, j, ...) {
  view = tag(obj,"view")
  ## need to popup a dialog to gather name and class, set view
  win = gwindow("Add column")
  group = ggroup(horizontal=FALSE, container=win)
  tbl = glayout()
  colName = gedit(paste("X",j+1,sep=""))
  ## logical is a problem in showing (shows as TRUE even if NA)
  colClass = gdroplist(c("numeric","character"))##,"factor"))##,"logical"))
  tbl[1,1] = glabel("Column name:")
  tbl[1,2] = colName
  tbl[2,1] = glabel("Column class")
  tbl[2,2] = colClass
  visible(tbl) <- TRUE
  
  add(group, tbl, expand=TRUE)
  buttonGroup = ggroup(container=group)
  addSpring(buttonGroup)
  gbutton("ok",container=buttonGroup, handler =function(h,...) {
    frame = obj[,,drop=FALSE]
    type = svalue(colClass)
    nRows = dim(obj)[1]
    if(type == "numeric") {
      x = rep(NA,length=nRows)
      x = as.numeric(x)
      newframe = data.frame(frame,x)
    } else if(type == "character") {
      x = character(dim(obj)[1])
      newframe = data.frame(frame,x,stringsAsFactors=FALSE)
    } else if(type == "factor") {
      x = character(dim(obj)[1])
      newframe = data.frame(frame,x,stringsAsFactors=TRUE)
    }
    names(newframe)[dim(obj)[2]+1] <- svalue(colName)
    obj[,] <- newframe
    
    ## need to set cursor here, as this happens after the setCursor below
    setCursorAtCell(obj, i,j+1, start.editing=TRUE)
    ## clean up
    dispose(win)
    return(TRUE)
  })
  gbutton("cancel", container=buttonGroup, handler=function(h,...) {
          dispose(win)
        })
}
  
##################################################

addColumn = function(obj, x, name=NULL) {
  store = .getRGtkDataFrame(obj)
  d = dim(obj)
  x = rep(x, length=d[1])               # recycle
  theColors = tag(obj,"theColors")
  fgColor = rep(theColors['fg'], length=d[1])
  bgColor = rep(theColors['bg'], length=d[1])

  toAdd = data.frame(x,fg=fgColor, bg=bgColor)
  if(is.character(x)) toAdd[,1] = as.character(x)
  for(i in 2:3) toAdd[,i] = as.character(toAdd[,i])
  store$AppendColumns(toAdd)

   if(is.null(name))
     name = paste("X",d[2]+1,sep="")
   view.col = addTreeViewColumnNoEdit(obj, d[2]+1, name)
}

addRow = function(obj, x, ...) {
  store = .getRGtkDataFrame(obj)
  d = dim(obj)
  dstore = dim(store)

  if(is.null(x))
    x = rep(NA,length=d[2])
  
  if(length(x) != d[2]) {
    warning("Need to add same size row or no row")
    return()
  }
  
  theRow = list(); theColors = tag(obj,"theColors")
  theRow[[1]] = TRUE; theRow[[2]]=""
  sapply(1:d[2], function(i) theRow[[3*i]] <<- x[[i]])
  sapply(1:d[2], function(i) theRow[[3*i + 1]] <<- theColors['fg'])
  sapply(1:d[2], function(i) theRow[[3*i + 2]] <<- theColors['bg'])

  store$AppendRows(theRow)

         
}


## subsetBy part
##################################################
## subset by widget -- specific to gDF object
##################################################

## here action, after values, gives an environment to evaluate variables within
setClass("gSubsetbyRGtk",
         contains="gComponentRGtk",
         prototype=prototype(new("gComponentRGtk"))
         )

setGeneric("gsubsetby",function(gridObj,
                   handler = NULL, action = NULL,        # when changed
                   container=NULL, ...) standardGeneric("gsubsetby"))
setMethod("gsubsetby",
          signature(gridObj = "gGridRGtk"),
          function(gridObj,
                   handler = NULL, action = NULL,        # when changed
                   container=NULL, ...) {
            
            vars = names(gridObj)
            group = ggroup(container = container, ...)
            subsetVar = gdroplist(c("NA",vars),  selected=1,container = group)
            subsetHow = gdroplist(c(""),  editable=TRUE, selected=1, container=group)
            size(subsetVar) <- c(25*8, 30)
            size(subsetHow) <- c(25*8, 30)
            leftArrow = gimage("larrow",dirname="stock",container = group)
            rightArrow = gimage("rarrow",dirname="stock",container = group)

            obj = new("gSubsetbyRGtk", block=group,widget=group,toolkit=gridObj@toolkit)
            
            tag(obj, "subsetVar") <- subsetVar
            tag(obj, "subsetHow") <- subsetHow
            tag(obj, "leftArrow") <- leftArrow
            tag(obj, "rightArrow") <- rightArrow
            tag(obj, "vars") <- vars
            tag(obj, "handler") <- handler
            tag(obj, "action") <- gridObj
            
            ## add handlers
            ## changing var name resets subsetHow
            addhandlerchanged(subsetVar, handler = function(h,...) {
              varName = svalue(subsetVar)
              if(!length(varName) || varName == "NA") {
                subsetHow[] = c("")
              } else {
                theColumn = which(varName == names(gridObj))
                theValues = gridObj[,theColumn, drop=TRUE]
                theValues = sort(unique(theValues))
                if(is.factor(theValues))
                  theValues = as.character(theValues)
                if(is.character(theValues))
                  theValues = paste(paste('"',theValues,sep=""), '"', sep = "") # quote
                subsetHow[] =  c("",paste("==",theValues,sep=" "))
              }
              svalue(subsetHow,index=TRUE) <- 1
            })
            ## changing subsetHow updates gridobject
            addhandlerchanged(subsetHow,handler = function(h,...) {
              how = svalue(subsetHow)
              if(is.empty(how)) {
                visible(gridObj) <- rep(TRUE, nrow(gridObj))
              } else {
                theColumn = which(svalue(subsetVar) == names(gridObj))
                theValues =  gridObj[,theColumn] # using name to extract column
                if(is.factor(theValues))
                  theValues = as.character(theValues)
                ## subsetHow of the form '== value'
                cmd = paste("theValues",svalue(subsetHow),collapse="")
                whichRows = gtktry(eval(parse(text=cmd)),silent=TRUE)
                if(!inherits( whichRows, "try-error")) {
                  whichRows[is.na(whichRows)] <- FALSE
                  visible(gridObj) <- whichRows
                }
              }
            })
            addhandlerclicked(leftArrow, handler = function(h,...) {
              subsetHow = tag(obj,"subsetHow")
              setValues = subsetHow[]
              curIndex = svalue(subsetHow,index=TRUE)
              n = length(setValues)
              if(is.na(curIndex)) curIndex = 2 # then newi = 1
              newIndex = (curIndex-2)%%n+1       # faster than ifelse?
              svalue(subsetHow,index=TRUE) <- newIndex
              return(TRUE)
            })
            addhandlerclicked(rightArrow, handler = function(h,...) {
              subsetHow = tag(obj,"subsetHow")
              setValues = subsetHow[]
              curIndex = svalue(subsetHow,index=TRUE)
              n = length(setValues)
              if(is.na(curIndex)) curIndex = n # then newi = 1
              newIndex = curIndex %% n + 1 ## really (i-1)+1 mod n + 1
              svalue(subsetHow,index=TRUE) <- newIndex
              return(TRUE)              
            })
            return(obj)
          })
          
##################################################
## methods
## this updates the names in subsetVar
setMethod("update",
          signature(obj="gSubsetbyRGtk"),
          function(object, ...) {
            obj = object                          # subsetby guy
            gridObj = tag(obj,"action")
            tag(obj,"subsetVar")[] <- c("NA",names(gridObj))
          })

setMethod("length",
          signature(x="gSubsetbyRGtk"),
          function(x) {
            gwCat("DEBUG: length called on gSubsetbyRGtk\n")
          })

## returns a vector of TRUE or FALSE
setMethod("svalue",
          signature(obj="gSubsetbyRGtk"),
          function(obj, index=NULL, drop=NULL, ...) {
            subsetVar = tag(obj, "subsetVar")
            subsetHow = tag(obj, "subsetHow")
            varName = svalue(subsetVar)
            if(!length(varName) || varName == "NA")
              return(NA)
            ## have a variable
            values = svalue(varName)
            assign(varName,values)
            condition = svalue(subsetHow)
            ret = eval(parse(text=Paste(varName, condition)))
            return(ret)
          })

## put onto the both widgets
setMethod("addhandlerchanged",
          signature(obj="gSubsetbyRGtk"),
          function(obj, handler=NULL, action=NULL, ...) {
            subsetVar = tag(obj, "subsetVar")
            subsetHow = tag(obj, "subsetHow")
            lst = list()
            lst[["subsetVar"]] <- addhandlerchanged(subsetVar, handler, action)
### It seems that this can cause a loop: subsetHow->ggrid->tvCol->subsetBy->subsetHow
            lst[["subsetHow"]] <- addhandlerchanged(subsetHow, handler, action)
### why is this buggy?
            ##            return(lst)                 # return IDS

          })

## redirect to above
setMethod(".addhandlerchanged",
          signature(toolkit="guiWidgetsToolkitRGtk2", obj="gSubsetbyRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            addhandlerchanged(obj, handler, action, ...)
          })

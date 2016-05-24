
#' Constructor for gdfedit widget
#'
#" gdfedit is a light-weight interface to RGtk2Extra's excellent data editor widget by Thomas Taverner

gdfedit <- function(
                    items = NULL,
                    name = paste(deparse(substitute(items)), "1", sep="."),
                    container = NULL, ... ,
                    toolkit=guiToolkit()
                    ) {

  widget <- .gdfedit (toolkit, items=items, name=name, container=container ,...) 
  obj <- new( 'guiComponent', widget=widget, toolkit=toolkit) 
  return(obj)
}

#' class for the widget
setClass("gDfEditRGtk",
         contains="gComponentRGtk",
         prototype=prototype(new("gComponentRGtk"))
         )



#' generic for toolkit dispatch
setGeneric( '.gdfedit' ,
           function(toolkit,
                    items = NULL, name = paste(deparse(substitute(items)),"1", sep="."),
                    container = NULL, ... )
           standardGeneric( '.gdfedit' ))


#' gWidgetsRGtk2 interface to constructor
setMethod(".gdfedit",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   items=NULL,
                   name,
                   container=NULL,...)  {
            ## main constructor
            ## should define so that getToolkitWidget() returns the gdfedit object
            
            ## load package RGtk2DfEdit and RGtk2 but get passed R CMD check by cheating
            gtkDfEdit <- NULL
            if(!do.call("require",list("RGtk2Extras")))
              stop(sprintf("Must have %s package installed", "RGtk2Extras"))
            do.call("require", list("RGtk2"))

            
            widget <- gtkDfEdit(items, dataset.name=name) ## how to add
            obj <- new("gDfEditRGtk", block=widget, widget=widget, toolkit=toolkit)

            
            if(!is.null(container)) {
              if(is.logical(container) && container) {
                container <- gwindow()
                add(container, obj)
              } else {
                add(container, obj, ...)
              }
            }

            return(obj)
          })


##' needed when adding as child
setOldClass("GtkDfEdit")
setMethod(".tag", signature(toolkit="guiWidgetsToolkitRGtk2",obj="GtkDfEdit"),
          function(obj, toolkit, i, drop=TRUE, ...) {
            NULL
          })


## data frame methods

#' return selected value
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gDfEditRGtk"),
          function(obj, toolkit, index=NULL, drop=NULL,...) {
            object <- getWidget(obj)
            df <- object$getDataFrame()
            sel <- object$getSelection()

            if(!is.null(index) && as.logical(index)) {
              return(sel$rows)
            }
            
            if(!is.null(drop) && as.logical(drop)) {
              ## use columns
              return(df[sel$rows, sel$columns])
            } else {
              return(df[sel$rows, ])
            }
          })


## ## set by index value selected value
## ## not sure what this is to do?
## setReplaceMethod(".svalue",
##                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gDfEditRGtk"),
##                 function(obj, toolkit, index=NULL, ..., value) {
##                   ## ??? might just skip this
##                   return(obj)
##                 })

setMethod("[",
         signature(x="gDfEditRGtk"),
          function(x, i, j, ..., drop=TRUE) {
            .leftBracket(x, x@toolkit, i, j,..., drop=drop)
          })

setMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gDfEditRGtk"),
          function(x, toolkit, i, j, ..., drop=TRUE) {
            object <- getWidget(x)
            df <- object$getDataFrame()
            
            if(missing(i) && missing(j))
              df[,,drop=drop]
            else if (missing(j))
              df[i,,drop=drop]
            else if(missing(i))
              df[,j, drop=drop]
            else
              df[i,j, drop=drop]
          })


## [<-
setReplaceMethod("[",
                 signature(x="gDfEditRGtk"),
                 function(x, i, j,..., value) {
                   x <- .leftBracket(x, x@toolkit, i, j,..., value)
                  return(x)
                 })

#' Method for [<-
#' really needs to have check on column type
setReplaceMethod(".leftBracket",
                 signature(toolkit="guiWidgetsToolkitRGtk2",x="gDfEditRGtk"),
                 function(x, toolkit, i, j, ..., value) {
                   ## we could implement -- but only if we do not change the class of the object
                   df <- getWidget(x)$getModel()
                   if(missing(i))
                     i <- seq_len(dim(x)[1])
                   if(missing(j))
                     j <- seq_len(dim(x)[2])
                   ## logical to index
                   if(is.logical(i)) 
                     i <- seq_along(i)[i]
                   if(is.logical(j))
                     j <- seq_along(j)[j]
                   
                   ## value if vector
                   if(is.null(dim(value))) {
                     if(length(j) != 1) {
                       cat(sprintf("Dimension of value does not match that of j\n"))
                       return(x)
                     }
                     try(df[i, 1 + j] <- value, silent=TRUE)
                     return(x)
                   }
                   ## value is not vector
                   if(dim(value)[1] != length(i) || dim(value)[2] != length(j)) {
                     cat(sprintf("Dimension of value does not match that of i or j\n"))
                     return(x)
                   }
                   ## fill in column by column
                   sapply(seq_along(j), function(l) {
                     try(df[i, 1 + j[l]] <- value[,j], silent=TRUE)
                   })
                   return(x)
                 })


## data frame like
setMethod(".dim",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gDfEditRGtk"),
          function(x,toolkit) {
            ## basically:
            object <- getWidget(x)
            data.model <- object$getDataFrame()
            dim(data.model)
          })


setMethod(".length",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gDfEditRGtk"),
          function(x,toolkit) {
            d <- dim(x)
            return(d[2])
          })


setMethod(".dimnames",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gDfEditRGtk"),
          function(x, toolkit) {
            ##
            object <- getWidget(x)
            rnames <- object$getRowNames()
            cnames <- object$getColumnNames()
            
            list(rnames, cnames)
          })

setReplaceMethod(".dimnames",
                 signature(toolkit="guiWidgetsToolkitRGtk2",x="gDfEditRGtk"),
                 function(x, toolkit,  value) {
                   object <- getWidget(x)
                   ##
                   rnames <- value[[1]]
                   cnames <- value[[2]]
                   if(!is.null(rnames)) {
                     model <- object$getModel()
                     model[,1] <- rnames
                   }
                   if(!is.null(cnames)) {
                     sapply(seq_along(cnames), function(i) {
                       object$setColumnName(i, cnames[i])
                     })
                   }
                   ##
                   x
                 })



setMethod(".names",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gDfEditRGtk"),
          function(x, toolkit) {
            dimnames(x)[[2]]
          })


setReplaceMethod(".names",
                 signature(toolkit="guiWidgetsToolkitRGtk2",x="gDfEditRGtk"),
                 function(x, toolkit, value) {
                   object <- getWidget(x)
                   sapply(seq_along(value), function(i) {
                     object$setColumnName(i, value[i])
                   })
                   x
                 })

## handlers to add
## for gdfedit -- change a value
#' @param handler function to call when column is clicked. Along with usual obj, action, first argument
#'        (typically h) has components df for the data frame, and column.no to return the column number that
#'        was clicked on
setMethod(".addhandlercolumnclicked",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gDfEditRGtk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            object <- getWidget(obj)
            
            f <- function(df, col.idx) {
              h <- list(obj=obj, action=action, df=df, column.no=col.idx)
              handler(h)
            }
            ## doesn't give back an ID!!
            object$setColumnClickHandler(f)
            invisible()
          })

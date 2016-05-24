setClass("gLayoutRGtk",
         contains="gContainerRGtk",
         prototype=prototype(new("gContainerRGtk"))
         )

## an gWidget for tables
 

setMethod(".glayout",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   homogeneous = FALSE,
                   spacing = 10,        # amount (pixels) between row, cols, NULL=0
                   container = NULL, ...
                   ) {
            
            force(toolkit)

            tbl <- gtkTableNew(homogeneous = homogeneous)
            ## homogeneous spacing
            tbl$SetRowSpacings(spacing)
            tbl$SetColSpacings(spacing)
            
            obj <- as.gWidgetsRGtk2(tbl)
            tag(obj, "childlist") <- list()
            
            if (!is.null(container)) {
              if(is.logical(container) && container == TRUE)
                container = gwindow()
              add(container, obj,...)
            }
            
            invisible(obj)
          })
            

as.gWidgetsRGtk2.GtkTable <- function(widget, ...) {
  obj = new("gLayoutRGtk", block=widget, widget=widget,
    toolkit=guiToolkit("RGtk2"))

  return(obj)
}


            
### The add method is a stub so that this works with same
## approach as gWidgetstcltk
setMethod(".add",
          signature(toolkit="guiWidgetsToolkitRGtk2", obj="gLayoutRGtk", value="gWidgetRGtk"),
          function(obj, toolkit, value, ...) {
            ## stub
          })


## retrieve values
setMethod("[",
          signature(x="gLayoutRGtk"),
          function(x, i, j, ..., drop=TRUE) {
            .leftBracket(x, x@toolkit, i, j, ..., drop=drop) 
          })
setMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gLayoutRGtk"),
          function(x, toolkit, i, j, ..., drop=TRUE) {
            l <- tag(x, "childlist")
            ind <- sapply(l, function(comp) {
              i[1] %in% comp$x && j[1] %in% comp$y
            })
            if(any(ind))
              return(l[ind][[1]]$child) # first
            else
              NA
          })


## how we populate the table
setReplaceMethod("[",
                 signature(x="gLayoutRGtk"),
                 function(x, i, j,..., value) {
                   .leftBracket(x, x@toolkit, i, j, ...) <- value
                   return(x)
                 })

setReplaceMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gLayoutRGtk"),
          function(x, toolkit, i, j, ..., value) {

            if(missing(i))
              i <- dim(x)[1] + 1
            if(missing(j)) {
              cat(gettext("glayout: [ needs to have a column specified."))
              return(x)
            }

            ## check that all is good
            if(is.character(value)) {
              ## wrap characters into labels
              value <- glabel(value,...)
            }

            ## widgets
            tbl <- getWidget(x)
            child <- getBlock(value)
            
            theArgs <- list(...)

            ## get expand, anchor, fill
            expand <- getWithDefault(theArgs$expand, FALSE)
            if(!is.null(theArgs$align))
              theArgs$anchor <- theArgs$align
            anchor <- getWithDefault(theArgs$anchor, NULL)
            if(!is.null(anchor)) {       # put in [0,1]^2
              anchor <- (anchor+1)/2      # [0,1]
              anchor[2] <- 1 - anchor[2]     # flip yalign
            }

            default_fill <- getWithDefault(tag(value, "default_fill"), "both")
            fill <- getWithDefault(theArgs$fill, default_fill) # "", x, y or both

              ## we do things differently if there is a gtkAlignment for a block
            if(is(child, "GtkAlignment")) {
              if(expand && (fill =="both" || fill == "x")) {
                child['xscale'] <- 1
              }

              if(expand && (fill == "both" || fill == "y")) {
                child['yscale'] <- 1
              }

              if(expand && fill == "") {
                child['xscale'] <- child['yscale'] <- 1
              }
              
              if(!is.null(anchor)) {
                child['xalign'] <- anchor[1]
                child['yalign'] <- anchor[2]
              }
            } else {
              ## in gtkstuff 
              setXYalign(child, getWidget(value), anchor)
            }

            ## fix up number of columns
            d <- dim(x)
            nr <- max(i); nc <- max(j)
            if( nr > d[1] || nc > d[2])
              tbl$Resize(max(max(i), nr), max(max(j), nc))
            
            if(expand)
              opts <- c("fill","expand","shrink")
            else
              opts <- c("fill")
            
            child <- getBlock(value)
            tbl$Attach(child,
                       min(j)-1, max(j), min(i)-1, max(i),
                       xoptions=opts,yoptions=opts)

            ## store for [ method
            l <- tag(x, "childlist")
            l[[as.character(length(l) + 1)]] <- list(x=i, y=j, child=value)
            tag(x, "childlist") <- l
            

            return(x)
          })

## inherits delete method for containers

## replaced
## We like visible, return it. Unlike delete it only hides the widget
## setReplaceMethod(".visible",
##                  signature(toolkit="guiWidgetsToolkitRGtk2",obj="gLayoutRGtk"),
##                  function(obj, toolkit, ..., value) {
##                    gwCat(gettext("visible<- is now redundant for glayout in RGtk2"))
##                    return(obj)
##                  })

## get number of rows and columns
setMethod(".dim", 
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gLayoutRGtk"),
          function(x,toolkit) {
            tbl <- getWidget(x)
            return(c(nrow=tbl$GetNrows(), ncol=tbl$GetNcols()))
          })

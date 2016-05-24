##' @include guiComponent.R

##' 

##' Method to add icon to list of stock icons
##'
##' @export
addStockIcons = function(iconNames,iconFiles, ..., toolkit = guiToolkit()) {
  out =  .addStockIcons (toolkit, iconNames, iconFiles, ...)
  return(out)
}

##' generic for dispath
##' @alias addStockIcons
setGeneric( '.addStockIcons' ,
           function(toolkit, iconNames, iconFiles,... )
           standardGeneric( '.addStockIcons' ))

##' return list of available stock icons
##'
##' @export
getStockIcons = function( ..., toolkit = guiToolkit()) {
  out =  .getStockIcons (toolkit,...)
  return(out)
}

##' generic for toolkit dispatch
##' @alias getStockIcons
setGeneric( '.getStockIcons' ,
           function(toolkit,...)
           standardGeneric( '.getStockIcons' ))

##' Find a stock icon from the given class
##'
##' @export
stockIconFromClass = function(theClass, ..., toolkit = guiToolkit()) {
  out =  .stockIconFromClass (toolkit, theClass, ...)
  return(out)
}
##' generic for dispath
##' @alias stockIconFromClass
setGeneric( '.stockIconFromClass' ,
           function(toolkit, theClass,... )
           standardGeneric( '.stockIconFromClass' ))

##' Find stock icon from the given object
##'
##' @export
stockIconFromObject = function(obj, ..., toolkit = guiToolkit()) {
  out =  .stockIconFromClass (toolkit, obj, ...)
  return(out)
}

##' generic for dispath
##' @alias stockIconFromObject
setGeneric( '.stockIconFromObject' ,
           function(toolkit, obj,... )
           standardGeneric( '.stockIconFromObject' ))







## returns a list with key the icon name
## and value the filepath





## find the stock icons. This includes those added bia loadGWidgetIcons()
gWidgetsIcons = list()
assignInNamespace("gWidgetsIcons",gWidgetsIcons,ns="gWidgets")
addedStockIcons = list()
assignInNamespace("addedStockIcons",addedStockIcons, ns="gWidgets")

getgWidgetsIcons = function() {
  gWidgetsIcons = getFromNamespace("gWidgetsIcons",ns="gWidgets")
  if(length(gWidgetsIcons) == 0) {
    path = system.file("images",package="gWidgets")
    allIcons = list.files(path)
    ## create a hash with name -> location
    for(i in allIcons) {
      filename = sub("\\.xpm$|\\.gif$|\\.jpg$|\\.jpeg$|\\.png$|\\.tiff$","",i)
      gWidgetsIcons[[filename]] <- system.file("images",i,package="gWidgets")
    }
  }
  assignInNamespace("gWidgetsIcons",gWidgetsIcons,ns="gWidgets")
  return(gWidgetsIcons)
}

## incorporate both these and any additional ones added via addStockIcons
setMethod(".getStockIcons",
          signature(toolkit="ANY"),
          function(toolkit) {
            gWidgetsIcons = getgWidgetsIcons()
            addedStockIcons = getFromNamespace("addedStockIcons",ns="gWidgets")

            return(c(addedStockIcons ,gWidgetsIcons))

            
          })

setMethod(".addStockIcons",
          signature(toolkit="ANY"),
          function(toolkit, iconNames, iconFiles, ...) {
            stockIcons = getgWidgetsIcons()
            addedStockIcons = getFromNamespace("addedStockIcons", ns="gWidgets")
            
            if(length(iconNames) == length(iconFiles)) {
              for(i in 1:length(iconNames)) {
                if(!is.na(match(iconNames[i],names(stockIcons))))
                  cat("Overriding stock icon",iconNames[i],"\n")
                addedStockIcons[[iconNames[i]]] <- iconFiles[[i]]
              }
              assignInNamespace("addedStockIcons", addedStockIcons, ns="gWidgets")
            } else {
              cat("Lengths of names and file don't match\n")
            }
            invisible()
          })

## name can be a vector
## return NA, if not there
getstockiconname = function(name=NULL) {
  .stockicons = getStockIcons()         # cache?

  if(is.null(name))
    return(unlist(.stockicons))
  

  tmpfun = function(names) {
    sapply(names, function(name) {
      ## already a stock name?
      if(name %in% .stockicons)
        return(name)
      
      if(name %in% names(.stockicons)) {
        return(.stockicons[[name]])
      } else {
        return(NA)
      }
    })
  }
  
  return(tmpfun(name))
}


#################################################
## functions to deal with icons
## class to icon translation -- return stock name
## with prefix

## find the stock icons. This includes those added bia loadGWidgetIcons()
setMethod(".stockIconFromClass",
          signature(toolkit="ANY"),
          function(toolkit,theClass, ...) {

            default = "symbol_star"
            
            if(is.null(theClass) ||
               is.na(theClass) ||
               length(theClass) == 0
               )
              return(NA)
            
            if(theClass %in% .models)
              return(getstockiconname("lines"))
            if(theClass %in% .ts)
              return(getstockiconname("ts"))
            if(theClass %in% .functions)
              return(getstockiconname("function1"))
            
            ret = switch(theClass,
              "complex"="numeric",
              "character"="character",
              "date" = "date",              
              "data.frame" = "dataframe",
              "integer"= "numeric",
              "factor"="factor",
              "function"="function1",
              "list" = "dataframe",
              "logical" = "logical",
              "matrix" = "matrix",
              "numeric"= "numeric",
              "recordedplot" = "plot",
              NA)
            
            return(getstockiconname(ret))
          })


setMethod(".stockIconFromObject",
          signature(toolkit="ANY"),
          function(toolkit,obj, ...) {
            .stockIconFromClass(class(obj)[1])
          })

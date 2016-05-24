require(methods)
require(utils)

## Set up generics and classes for ANY widgets.                                 
## ANY widgets are meant to work with any toolkit. They are compound widgets
## their methods are not inherited from the toolkit ones, but rather are made explicit here.


## how to make a generic widget here?



### Methods
## We need these methods to get the dispatch correct.



## TAG method uses ID
## This is taken from gWIdgets tcltk. It uses a big hash. Ughh
## create namespace object
tags = list()
assignInNamespace("tags",list(),"gWidgets")

## clear out tags for this ID. Not exported. Is this used?
Tagsclear = function(obj) {

  id = obj@ID
  
  tags = getFromNamespace("tags",ns="gWidgets")
  allKeys = names(tags)

  inds = grep(paste("^",id,"::",sep=""),allKeys)
  if(length(inds) == 0)
    return(NA)

  ## else
  tags[[inds]] <- NULL
  assignInNamespace("tags",tags,ns="gWidgets")
}


setMethod("tag",signature(obj="gWidgetANY"),
          function(obj,i,drop=TRUE, ...) {
            if(missing(drop)) drop <- TRUE
            .tag(obj, obj@toolkit,i, drop=drop,...)
          })

setMethod(".tag", signature(toolkit="ANY",obj="guiWidget"),
          function(obj, toolkit, i, drop=TRUE, ...) {
            if(missing(i)) i = NULL
            if(missing(drop)) drop <- TRUE                        
            .tag(obj@widget,toolkit,  i, drop=drop,  ...)
          })
setMethod(".tag", signature(toolkit="ANY",obj="gWidgetANY"),
          function(obj, toolkit, i, drop=TRUE, ...) {
            if(missing(i)) i = NULL
            if(missing(drop)) drop <- TRUE                                    


            id = obj@ID

            ## get all values for this id
            tags = getFromNamespace("tags",ns="gWidgets")
            allKeys = names(tags)

            inds = grep(paste("^",id,"::",sep=""),allKeys)
            if(length(inds) == 0)
              return(NULL)

            justTheKeys = sapply(allKeys[inds],function(keyWithID) {
              sub(paste("^",id,"::",sep=""),"",keyWithID)
            })

            tagByKey = list()
            for(key in justTheKeys) 
              tagByKey[[key]] = tags[[paste(id,key,sep="::")]]
                      
            
            
            if(is.null(i)) return(tagByKey)

            if(drop) {
              if(length(i) == 1)
                return(tagByKey[[i]])
              else
                return(sapply(i, function(j) tagByKey[j]))
            } else {
              return(sapply(i, function(j) tagByKey[j]))
            }
          })

## tag <-
setReplaceMethod("tag",signature(obj="gWidgetANY"),
          function(obj, i, replace=TRUE, ..., value) {
            .tag(obj, obj@toolkit,i,replace, ...) <- value
            return(obj)
          })

## objects can be in many different flavors: guiWIdget, gWidgettcltk, tcltkObject
setReplaceMethod(".tag", signature(toolkit="ANY",obj="guiWidget"),
          function(obj, toolkit, i, replace=TRUE, ..., value) {
            if(missing(i)) i = NULL
            .tag(obj@widget,toolkit,  i, replace, ...) <- value
            return(obj)
          })

setReplaceMethod(".tag", signature(toolkit="ANY",obj="gWidgetANY"),
          function(obj, toolkit, i, replace=TRUE, ..., value) {
            if(missing(i)) i = NULL
            

            id = obj@ID
            key = paste(id,i,sep="::")
            
            ## if we append we need to work a little harder
            tags = getFromNamespace("tags",ns="gWidgets")
  
            if(replace==FALSE) {
              value = c(tags[[key]],value)
            }

            tags[[key]] <- value
            assignInNamespace("tags", tags,ns="gWidgets")

            return(obj)

          })

### svalue


## Generic svalue
setMethod(".svalue",signature(toolkit = "ANY", obj="character"),
          function(obj, toolkit, index=NULL, drop=NULL,  ...)  {
            ifelse(length(obj) == 1,
                   return(getObjectFromString(obj)),
                   return(NA)
                   )
          })

setMethod(".svalue",signature(toolkit = "ANY", obj="NULL"),
          function(obj, toolkit, index=NULL, drop=NULL,  ...)  {
            return(NULL)
          })


setMethod("svalue",signature(obj="gWidgetANY"),
          function(obj, index=NULL, drop=NULL, ... ) {
            toolkit = obj@toolkit
            .svalue(obj, toolkit, ...,index=index, drop=drop)
          })

## svalue<- -- objec specific
setReplaceMethod("svalue",signature(obj="gWidgetANY"),
          function(obj, index=NULL, ...,value) {
            .svalue(obj, obj@toolkit, index=index, ...) <- value
            return(obj)
          })

## [. [<-
## [
setMethod("[",
          signature(x="gWidgetANY"),
          function(x,i,j,...,drop=TRUE) {
            
            return(.leftBracket(x, x@toolkit,i,j,...,drop=TRUE))
          })

## [<-
setReplaceMethod("[",signature(x="gWidgetANY"),
          function(x,i,j,...,value) {
            if(missing(i) && missing(j))
              .leftBracket(x, x@toolkit,...) <- value
            else if(missing(j))
              .leftBracket(x, x@toolkit,i,...) <- value
            else 
              .leftBracket(x, x@toolkit,i,j,...) <- value
            return(x)
          })

## size<-
setReplaceMethod("size",signature(obj="gWidgetANY"),
          function(obj, ..., value) {
            .size(obj, obj@toolkit,...) <- value
            return(obj)
          })

## visible
setMethod("visible",signature(obj="gWidgetANY"),
          function(obj, set=NULL, ...) {
            .visible(obj,obj@toolkit, set=set, ...)
          })

setReplaceMethod("visible",signature(obj="gWidgetANY"),
          function(obj, ..., value) {
            .visible(obj, obj@toolkit, ...) <- value
            return(obj)
          })

## enabled -- TRUE If state is normal
setMethod("enabled",signature(obj="gWidgetANY"),
          function(obj, ...) {
            .enabled(obj, obj@toolkit,...)
          })

setReplaceMethod("enabled",signature(obj="gWidgetANY"),
          function(obj, ..., value) {
            .enabled(obj, obj@toolkit,...) <- value
            return(obj)
          })

setMethod("focus",signature(obj="gWidgetANY"),
          function(obj, ...) {
            .focus(obj, obj@toolkit,...)
          })


setReplaceMethod("focus",signature(obj="gWidgetANY"),
          function(obj, ..., value) {
            .focus(obj, obj@toolkit,...) <- value
            return(obj)
          })

## font<-
setReplaceMethod("font",signature(obj="gWidgetANY"),
          function(obj, ..., value) {
            .font(obj, obj@toolkit,...) <- value
            return(obj)
          })

## update
setMethod("update",signature(object="gWidgetANY"),
          function(object, ...) {
            .update(object, object@toolkit, ...)
          })


## dimnames
setMethod("dimnames", "gWidgetANY", function(x) .dimnames(x,x@toolkit))
setReplaceMethod("dimnames",
                 signature(x="gWidgetANY"),
                 function(x,value) {
                   .dimnames(x,x@toolkit) <- value
                   return(x)
                 })

## as of 2.5.0 this became primiive
if(as.numeric(R.Version()$major) <= 2 &
   as.numeric(R.Version()$minor) <= 4.1) {
  setGeneric("names")
  setGeneric("names<-")
}

setMethod("names", "gWidgetANY", function(x) .names(x,x@toolkit))
setReplaceMethod("names",
                 signature(x="gWidgetANY"),
                 function(x,value) {
                   .names(x,x@toolkit) <- value
                   return(x)
                 })


## add for widgets (gtext, ghelp, ...)
## Container methods are in toolkits implementations
setMethod("add",signature(obj="ANY"),
          function(obj, value, ...) {
            .add(obj, obj@toolkit,value,...)
          })

setMethod(".add",
          signature(toolkit="ANY",
                    obj="ANY", value="gWidgetANY"),
          function(obj, toolkit, value, ...) {
            .add(obj, toolkit, value@widget, ...)
          })



## undo and redo only implemented in some toolkits
setMethod("undo",signature(obj="gWidgetANY"),
          function(obj,...) {
            .undo(obj, obj@toolkit,...)
          })

setMethod(".undo", signature(toolkit="ANY",obj="guiWidget"),
          function(obj, toolkit, ...) {
            .undo(obj@widget,toolkit,...)
          })
setMethod(".undo", signature(toolkit="ANY",obj="gWidgetANY"),
          function(obj,toolkit, ...) {
            ## nothing
          })

setMethod("redo",signature(obj="gWidgetANY"),
          function(obj,...) {
            .redo(obj, obj@toolkit,...)
          })

setMethod(".redo", signature(toolkit="ANY",obj="guiWidget"),
          function(obj, toolkit, ...) {
            .redo(obj@widget,toolkit)
          })
setMethod(".redo", signature(toolkit="ANY",obj="gWidgetANY"),
          function(obj, toolkit, ...) {
            ## nothing
          })

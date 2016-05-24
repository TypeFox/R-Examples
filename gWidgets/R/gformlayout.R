## a way to layout out a form (dialog) using a list, similar to ggenericwidget
## but perhaps better

## idea is taken from extjs www.extjs.com


##' @include guiComponents.R

##' Class for widget to layout form from a specification in a list
setClass("gFormLayout",
         contains="guiComponent",
         prototype=prototype(new("guiComponent"))
         )

##' A constructor for a widget to layout a form from a specification in a list
##'
##' @export
gformlayout <- function(
                        lst, container = NULL,
                        ... ,
                        toolkit=guiToolkit()){
  widget <- .gformlayout(toolkit,
                         lst=lst, container=container ,...
                         )
  obj <- new( 'gFormLayout',widget=widget,toolkit=toolkit) 
  return(obj)
}


##' generic for toolkit dispatch
##' @alias gformlayout
setGeneric( '.gformlayout' ,
           function(toolkit,
                    lst, 
                    container = NULL, ... )
           standardGeneric( '.gformlayout' ))



##################################################
## ANY implementation
setClass("gFormLayoutANY",
         representation = representation("gComponentANY",
           lst="list"),
         contains="gComponentANY",
         prototype=prototype(new("gComponentANY"))
         )

###############
setMethod(".gformlayout",
          signature(toolkit="ANY"),
          function(toolkit,
                   lst,
                   container=NULL,
                   ...
  ) {

            ## see docs for description of the list
            e <- new.env()
            mainGroup <- ggroup(container = container, ...)
            .makeForm(lst, mainGroup, e)
            
            ## return the container now that it has all the stuff in in.
            obj = new("gFormLayoutANY",
              block=mainGroup, widget=mainGroup,
              ID = getNewID(), toolkit=toolkit,
              lst = lapply(e, function(i) i))
            
            invisible(obj)
          })

### methods
setMethod("[",
          signature(x="gFormLayoutANY"),
          function(x, i, j, ..., drop=TRUE) {
            .leftBracket(x, x@toolkit, i, j, ..., drop=drop)
          })
## return list containing widgets
setMethod(".leftBracket",
          signature(toolkit="ANY",x = "gFormLayoutANY"),
          function(x, toolkit, i, j, ..., drop = TRUE) {
            widgets <- x@lst
            if(missing(i))
              return(widgets)
            else if(length(i) == 1 && drop == TRUE)
              return(widgets[[i]])
            else
              return(widgets[i])
          })

## return a list of values
setMethod(".svalue",
          signature(toolkit="ANY",obj="gFormLayoutANY"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {
            return(lapply(obj[], svalue))
          })

setMethod(".names",
          signature(toolkit="ANY",x="gFormLayoutANY"),
          function(x, toolkit) {
            return(names(x[]))
          })
          
## changes here should be copied into gWidgetsWWW
## helper functions
.makeForm <- function(lst, parent, e, ...) {
  g <- ggroup(container = parent, expand=TRUE,...)

  ## make a local copy of lst and modify for do.call
  tmp <- lst;
  tmp$name <- tmp$type <- tmp$children <- NULL
  tmp$depends.on <- tmp$depends.FUN <- tmp$depends.signal <- NULL
  tmp$container <- g; tmp$expand <- TRUE

  ## treat fieldset differently
  if(lst$type == "fieldset") {
    .makeFieldset(lst, g, e, label = lst$label)
    return()
  } else {  
    ## make object
    newObject <- do.call(lst$type, tmp)
    ## store if a name is given
    if(!is.null(lst$name))
      e[[lst$name]] <- newObject
    ## do we enable new object
    if(!is.null(lst$depends.on)) {
      widget <- e[[lst$depends.on]]
      if(is.null(lst$depends.signal))
        lst$depends.signal <- "addHandlerChanged"
      do.call(lst$depends.signal,list(obj=widget,handler =  function(h,...) {
        value <- svalue(h$obj)
        enabled(newObject) <- lst$depends.FUN(value)
      }))
      enabled(newObject) <- lst$depends.FUN(svalue(widget))
    }
  }
   


  
  ## show children if there
  ## this recurses except on "fieldset"
  if(!is.null(lst$children)) {
    for(i in 1:length(lst$children)) {
      l <- lst$children[[i]]
      if(l$type == "fieldset") {
        if(lst$type == "gnotebook")
          .makeFieldset(l, newObject, e, label = l$label)
        else
          .makeFieldset(l, newObject, e)
      } else {
        if(lst$type == "gnotebook")
          .makeForm(l, newObject, e, label = l$label)
        else
          .makeForm(l, newObject, e)
      }
    }
  }
}

## fieldset does not recurse
.makeFieldset <- function(lst, parent, e, ...) {
  ## parent is parent container
  ## lst is list as above

  ## outer container
  if(!is.null(lst$label)) 
    g <- gframe(lst$label, container = parent, ...)
  else
    g <- ggroup(container = parent, ...)
  ## main table
  tbl <- glayout(container = g)
  
  ## do we enable new object
  if(!is.null(lst$depends.on)) {
    widget <- e[[lst$depends.on]]
    if(is.null(lst$depends.signal))
      lst$depends.signal <- "addHandlerChanged"
    do.call(lst$depends.signal, list(obj = widget,handler = function(h,...) {
      value <- svalue(h$obj)
      enabled(g) <- lst$depends.FUN(value) 
  }))
    enabled(g) <- lst$depends.FUN(svalue(widget))
  }
  
  ## fix label adjust
  if(is.null(lst$label.pos))
    lst$label.pos <- "left"
  if(lst$label.pos == "top") {
    label.anchor <- c(-1,0)
  } else {
    if(is.null(lst$label.just) || lst$label.just == "left")
      label.anchor <- c(-1,1)
    else if(lst$label.just == "center")
      label.anchor <- c(0,1)
    else
      label.anchor <- c(1,1)
  }
  
  if(is.null(lst$columns)) 
    no.columns <- 1
  else
    no.columns <- lst$columns
  
  ## add children
  for(i in 1:length(lst$children)) {
    l <- lst$children[[i]]
    ## each child is a list with name, label, type, then arguments
    ## make new list for do.call
    tmp <- l;
    tmp$name <- tmp$label <- tmp$type <- NULL
    tmp$depends.on <- tmp$depends.FUN <- tmp$depends.signal <- NULL
    tmp$container <- tbl

    newWidget <- do.call(l$type, tmp)

    ## store
    if(!is.null(l$name))
      e[[l$name]] <- newWidget
    ## do we enable new object
    if(!is.null(l$depends.on)) {
      widget <- e[[l$depends.on]]
      if(is.null(l$depends.signal))
        l$depends.signal <- "addHandlerChanged"
      do.call(l$depends.signal, list(obj = widget, handler =  function(h,...) {
        value <- svalue(h$obj)
        enabled(newWidget) <- l$depends.FUN(value)
      }))
      enabled(newWidget) <- l$depends.FUN(svalue(widget))
    }

    
    ## add to table
    col <- 1 + (i - 1) %% no.columns    #1, ..., no.columns
    row <- 1 + (i - 1) %/% no.columns   #1, ...
    newLabel <- glabel(l$label, container = tbl)
    if(!is.null(lst$label.font))
      font(newLabel) <- lst$label.font
    if(is.null(lst$label.pos) || lst$label.pos == "left") {
      tbl[row, 2 * (col - 1) + 1, anchor=label.anchor] <- newLabel
      if(l$type %in% c("gcombobox","gdroplist"))
        tbl[row, 2 * (col - 1) + 2, anchor=c(-1,1), expand=TRUE] <- newWidget
      else
        tbl[row, 2 * (col - 1) + 2, anchor=c(-1,1)] <- newWidget
    } else {
      tbl[2 * (row - 1) + 1, col, anchor=label.anchor] <- newLabel
      if(l$type %in% c("gcombobox","gdroplist"))
        tbl[2 * (row - 1) + 2, col, anchor=c(-1,1), expand=TRUE] <- newWidget
      else
        tbl[2 * (row - 1) + 2, col, anchor=c(-1,1)] <- newWidget
    }
  }
}


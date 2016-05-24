##' @include S4-methods.R
NULL

## Methods using a previously defined S3 generic
## Mostly data store methods

##' Length of object
##'
##' @inheritParams base::length
##' @return length of object
##' @export
##' @rdname gWidgets2-S3methods
##' @method length GComponent
##' @S3method length GComponent
length.GComponent <- function(x) {
  if(isExtant(x))
    x$get_length()
}

##' set length of object
##'
##' @param x component
##' @param value value to assign
##' @export
##' @usage \method{length}{GComponent} (x) <- value
##' @rdname gWidgets2-S3methods
##' @method length<- GComponent
##' @S3method length<- GComponent
"length<-.GComponent" <- function(x, value) {
  if(isExtant(x))
    x$set_length(value)
  x
}

##' dimension of object
##'
##' @export
##' @rdname gWidgets2-S3methods
##' @method dim GComponent
##' @S3method dim GComponent
dim.GComponent <- function(x) {
  if(isExtant(x))
    x$get_dim()
}


##' get names of object
##'
##' Names are used in many different contexts.
##' @export
##' @rdname gWidgets2-S3methods
##' @method names GComponent
##' @S3method names GComponent
names.GComponent <- function(x) {
  if(isExtant(x))
    x$get_names()
}

##' set names of object
##'
##' @export
##' @usage \method{names}{GComponent} (x) <- value
##' @rdname gWidgets2-S3methods
##' @method names<- GComponent
##' @S3method names<- GComponent
"names<-.GComponent" <- function(x, value) {
  if(isExtant(x))
    x$set_names(value)
  x
}

##' Get dimnames of object
##'
##' @export
##' @rdname gWidgets2-S3methods
##' @method dimnames GComponent
##' @S3method dimnames GComponent
dimnames.GComponent <- function(x) {
  if(isExtant(x))
    x$get_dimnames()
}

##' Set dimnames of object
##'
##' @export
##' @usage \method{dimnames}{GComponent} (x) <- value
##' @rdname gWidgets2-S3methods
##' @method dimnames<- GComponent
##' @S3method dimnames<- GComponent
"dimnames<-.GComponent" <- function(x, value) {
  if(isExtant(x))
    x$set_dimnames(value)
  x
}

##' Get items of object
##'
##' We use the extraction operator, \code{[}, typically to refer to
##' the underlying items from which a selection can be made. As well,
##' we overload this to containers to refer to the child components.
##' @param i index or row index if applicable
##' @param j column index if applicable
##' @param drop logical. Does return value get "dropped" down to something easier?
##' @param ... dots argument
##' @export
##' @rdname gWidgets2-S3methods
##' @method [ GComponent
##' @S3method [ GComponent
"[.GComponent" <- function(x, i, j, ..., drop=TRUE) {
  if(isExtant(x))
    x$get_items(i, j, ..., drop=drop)
}

##' Return children of a parent container
##'
##' @export
##' @rdname gWidgets2-S3methods
##' @method [ GContainer
##' @S3method [ GContainer
"[.GContainer" <- function(x, i, j, ..., drop=TRUE) {
  if(isExtant(x))
    x$get_items(i, j, ..., drop=TRUE)
}

##' Set object's items
##'
##' @export
##' @usage \method{[}{GComponent} (x, i, j, ...) <- value
##' @rdname gWidgets2-S3methods
##' @method [<- GComponent
##' @S3method [<- GComponent
"[<-.GComponent" <- function(x, i, j, ..., value) {
  if(isExtant(x)) 
    x$set_items(value, i, j, ...)
  x
}

### This is an issue with the redesign
## ##' $ -- get property from underlying widget
## "$.GComponent" <- function(x, key, ...) x$get_property(key, ...)

## ##' $<- set property of underlying widget
## "$<-.GComponent" <- function(x, key, ..., value) {
##   x$set_property(key, ..., value=value)
##   x
## }

## "[[.GComponent" <- function(x, i, ...) 
##   x$get_property(i, ...)

## "[[<-.GComponent" <- function(x, i, ..., value)  {
##   x$set_property(key, ..., value=value)
##   x
## }

##' Call widgets \code{update_widget} method
##'
##' The update method will ca
##' use a widget to recompute itself, if it is necessary.
##' @param object object to update
##' @export
##' @rdname gWidgets2-S3methods
##' @method update  GComponent
##' @S3method update  GComponent
update.GComponent <- function(object, ...) {
  if(isExtant(object))
    object$update_widget(...)
}


##' str method for widgets
##'
##' @export
##' @rdname gWidgets2-S3methods
##' @method str GComponent
##' @S3method str GComponent
str.GComponent <- function(object, ...) cat(sprintf("Object of class %s\n", class(object)))



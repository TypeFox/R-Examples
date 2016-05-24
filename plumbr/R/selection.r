## A selection is of a set of observations or a region in data space.

## How to handle weighted selections? This is only relevant to
## observations. It does not make sense to have a separate logical
## selection vector when there are weights, so we accept either a
## logical or numeric vector for a selection of observations. This
## yields a fairly natural, implicit equivalence, i.e., one could
## easily treat a logical selection as weights. Factors should also be
## allowed. Pretty much any atomic vector.

## Regions are matrices, one column per dimension, one row for each
## point. One row is a single point. Two rows implies a line or
## rectangle. More than two represents a polygon. If a row contains an
## NA, it begins a new polygon. If categorical variables are involved,
## the levels will be lost in coercion to numeric. Problem?

## Is it necessary to have a separate class for each? Or do we have a
## single class, with a separate accessor for a list of observations
## and a region, where only one return value is ever non-empty?
## Probably should go with two separate classes, so that it is
## formally indicated which type of selection will be occurring. Then
## we define coercion methods that yield vectors and matrices. The
## underlying selection (.selection) should be or yield something that
## is coercible to the appropriate type.

## Selection is a data model, and proxies may be defined that perform
## linking and other transformations of the selection. The root model
## typically stores the selection, or retrieves it from some other
## data structure. Note that the selection could also be stored in the
## data; here we are storing it on a per brush/widget basis. Is it
## necessary to subclass Selection to implement a proxy? Or can the
## user provide a delegate function? Probably easier for user to
## provide a delegate.

### Abstract Selection ###

## This active binding hides the dynamic nature of .selection
selection_binding <- function(val) {
  if (missing(val)) {
    if (is.function(.selection))
      .selection()
    else .selection
  } else {
    if (is.function(.selection))
      .selection(val)
    else {
      .self$.selection <- val
      changed$emit()
    }
  }
}

Selection.gen <- setRefClass("Selection",
                             fields = c(
                               .selection = "ANY",
                               declareSignal(changed()),
                               selection = selection_binding
                               ))

### Accessors ###

Selection.gen$methods(replace = function(x) {
  selection <<- x
})

### Linking and Scaling ###

Selection.gen$methods(link = function(linker) {
  if (is.integer(linker)) {
    matching <- linker
    linker <- function(x) as.vector(x)[matching]
  }
  if (is.function(linker)) {
    delegate <- function(val) {
      if (missing(val))
        linker(.self)
      else .self$replace(linker(.self, val))
    }
    proxy <- ItemSelection(delegate)
### FIXME: should probably see if the linked selection has really changed first
    changed$connect(proxy$changed$emit)
    proxy
  } else stop("'linker' must be either a function or an integer matching")
})

## the scaler takes a selection and modifies another object, usually a
## data pipeline or perhaps another selection model (plot-plot)
Selection.gen$methods(scale = function(scaler, data) {
  invisible(changed$connect(function() scaler(.self, data)))
})

### Item and Region Selections ###

ItemSelection.gen <- setRefClass("ItemSelection", contains = "Selection")
RegionSelection.gen <- setRefClass("RegionSelection", contains = "Selection")

##' @param delegate An object to which the model delegates for
##' obtaining the selection status. This can be either a logical
##' vector (\code{TRUE} where selected) or a function with one
##' argument. If the function is called with no arguments, it should
##' return the selection. Otherwise, the argument is the new
##' selection status, and the function should store it. This is the
##' same semantic as \link[=makeActiveBinding]{active bindings}.
##' @noRd
##' @exportClass Selection, ItemSelection, RegionSelection
##' @export ItemSelection RegionSelection DataSelection
ItemSelection <- function(delegate = NULL) {
  ItemSelection.gen$new(.selection = delegate, changed = Signal())
}
RegionSelection <- function(delegate = NULL) {
  RegionSelection.gen$new(.selection = delegate, changed = Signal())
}

##' @exportMethod as.integer, as.numeric, as.logical, as.factor, as.matrix, which
setAs("ItemSelection", "vector", function(from) as.vector(x))
setMethod("as.vector", "ItemSelection", function(x) as.vector(x$selection))

setAs("ItemSelection", "integer", function(from) as.integer(x))
setMethod("as.integer", "ItemSelection", function(x) as.integer(x$selection))

setAs("ItemSelection", "numeric", function(from) as.numeric(x))
setMethod("as.numeric", "ItemSelection", function(x) as.numeric(x$selection))

setAs("ItemSelection", "logical", function(from) as.logical(x))
setMethod("as.logical", "ItemSelection", function(x) as.logical(x$selection))

setAs("ItemSelection", "factor", function(from) as.factor(x))
setMethod("as.factor", "ItemSelection", function(x) as.factor(x$selection))

setAs("RegionSelection", "matrix", function(from) as.matrix(x))
setMethod("as.matrix", "RegionSelection", function(x) as.matrix(x$selection))

### Convenience ###

setMethod("which", "ItemSelection", function(x) which(as.logical(x$selection)))

### Selection Calculus ###

## Do we want a binary operator syntax like +, -, & and ^ (instead of
## xor)?  This would imply a copy. What we really want are assignment
## operators like %+=%. A copy *might* be useful when trying to track
## selection history, but that could be recorded through a pipeline
## history mechanism or by listening to selection changes and saving
## the vector representation.

## Or we could simply use methods like: x$add(i), x$subtract(i),
## x$toggle(i) and x$intersect(i). This is probably better than a
## weird %+=%. Also need x$replace(i), because we do not want to user
## accessing the 'selection' field directly (a selection of a
## selection?) or using []<- (subassignment not replacement, in general).

ItemSelection.gen$methods(add = function(x) {
  if (is.logical(selection) && is.logical(x)) # most common case
    selection[x] <<- TRUE
  else
    selection <<- as(as.numeric(selection) + as.numeric(x), class(selection))
})

ItemSelection.gen$methods(subtract = function(x) {
  if (is.logical(selection) && is.logical(x)) # most common case
    selection[x] <<- FALSE
  else
    selection <<- as(as.numeric(selection) - as.numeric(x), class(selection))
})

ItemSelection.gen$methods(toggle = function(x) {
  x <- as.logical(x)
  selection[x] <<- !as.logical(selection[x])
})

ItemSelection.gen$methods(intersect = function(x) {
  x <- as.logical(x)
  selection <<- x & as.logical(selection)
})

## anything beyond 'add' is too much work without e.g. Qt
## we just add another region, possibly overlapping
RegionSelection.gen$methods(add = function(x) {
  selection <<- rbind(as.matrix(selection), NA, as.matrix(x))  
})

### Data Selection ###

##' Implement a selection model against a dataset/pipeline
##'
##' @title Selection in Data
##' @param data \code{\link{mutaframe}} of the dataset/pipeline
##' @param column Column index of selection variable in data
##' @return An \code{\link{ItemSelection}} reflecting the selection in the data
##' @author Michael Lawrence
DataSelection <- function(data, column = 1L) {
  selection <- ItemSelection(function(val) {
    if (missing(val))
      data[[column]]
    else data[[column]] <- as.vector(selection)
  })
  if (is.mutaframe(data)) {
    add_handler(data, function(i, j) {
      if (column %in% j)
        selection$changed$emit()
    })
  }
  selection
}


### =========================================================================
### PropertySet objects
### -------------------------------------------------------------------------

##' The \code{PropertySet} class is a collection of properties and is
##' useful as a data model, e.g., for storing the parameters of some
##' operation.
##'
##' \code{PropertySet} object has following methods, where \code{x} is
##' a \code{PropertySet} object:
##' 
##' \describe{
##'   \item{}{\code{x$properties()}: Return the classes of the properties as a
##'     named character vector. Compare to the \code{fields} method on
##'     \link[methods:setRefClass]{a reference class generator}}.
##'   \item{}{\code{as.list(x)}: Returns a named list of the property values.}
##' }
##'
##' When any property in the set changes, the \code{changed(name)}
##' signal is emitted, where \code{name} is the name of the property
##' that changed.
##' @example objectProperties/inst/examples/PropertySet.R
##' @rdname PropertySet-class
##' @aliases setPropertySet
##' @name PropertySet-class
##' @title PropertySet-class
##' @author Michael Lawrence, Tengfei Yin
##' @rdname PropertySet-class
setRefClass("PropertySet", contains = "VIRTUAL",
            fields = c(declareSignal(changed(name))),
            methods = list(
              properties = function() {
                fieldClasses <- getRefClass()$fields()
                fieldNames <- names(fieldClasses)
                fieldNames <- grep("^\\.init\\.",
                                   grep("^\\.", fieldNames, value = TRUE),
                                   invert = TRUE, value = TRUE)
                fieldClasses <- fieldClasses[fieldNames]
                names(fieldClasses) <- sub("^\\.", "", fieldNames)
                fieldClasses[!fieldClasses %in% "Signal"]
              },
              initialize = function(...) {
                ### qualifier is necessary for subclasses in other
                ### packages to work, unless we export it
                objectProperties:::PropertySet_initSignals(.self)
                callSuper(...)
              }))

## Connect and synchronize global/individual property changed signals.
##
## Factored this out to avoid annoying static code analysis warning
## from sharing a variable between closures.
PropertySet_initSignals <- function(x) {
  ## forward individual signals to global signal
  propNames <- names(x$properties())
  if (!length(propNames))
    return()
  signalNames <- paste(propNames, "Changed", sep = "")
  lapply(signalNames, function(signalName) {
    signal <- x[[signalName]]
    signal$connect(function() {
      skipOneEmission <<- TRUE # avoids infinite recursion
      x$changed$emit(signalName)
    })
  })
  ## forward global signal to individual signals
### FIXME: is.environment(x) is TRUE, but coercion seems necessary here - why?
  signals <- mget(signalNames, as.environment(x))
  skipOneEmission <- FALSE
  x$changed$connect(function(name) {
    signal <- signals[[name]]
    if (is.null(signal))
      stop("Unknown signal: ", name)
    if (!skipOneEmission)
      signal$emit()
    else skipOneEmission <<- FALSE
  })
}

##' Convenience function for defining a set of reference class fields
##' that signals when set.
##' 
##' When constructing signaling fields in this way, each field has the
##' ability to register its own signal and at the same time, there is
##' one top level signal which could be emitted no matter which field
##' changes. Please see the example to learn to register global signal
##' and individual signal.
##' @title Properties signaling fileds 
##' @param fields list of names of the field and associated fields class 
##' @param prototype A list of values declaring a default value for a field.
##' @return A list that is easily concatenated into the field list
##' @author Michael Lawrence, Tengfei Yin
##' @example objectProperties/inst/examples/properties.R
##' @export
properties <- function(fields = list(), prototype = list())
{
  if (!length(fields))
    return(list())
  .fieldNames <- paste(".", names(fields), sep = "")
  .initNames <- paste(".init.", names(fields), sep = "")
  hasPrototype <- names(fields) %in% names(prototype)
  if (any(sapply(fields, is.function) & hasPrototype))
    stop("An active binding field cannot have a prototype")
  if (length(invalidPrototypes <- setdiff(names(prototype), names(fields))))
    stop("Prototypes without matching property: ",
         paste(invalidPrototypes, collapse = ", "))
  activeFields <- mapply(function(fieldClass, fieldName, .fieldName, initName,
                                  hasPrototype, prototype, thisSignal)
  {
    as.function(c(alist(val=), substitute({
      if (missing(val)) {
        if (hasPrototype && !length(initName)) {
          .fieldName <<- prototype
          initName <<- TRUE
        }
        .fieldName
      } else {
        if (!is.function(fieldClass)) {
          coercedVal <- try(as(val, fieldClass, strict = FALSE), silent = TRUE)
          if (inherits(coercedVal, "try-error"))
            stop("Cannot set an object of type '", class(val), "' on '",
                 fieldName, "', a property of type '", fieldClass, "'")
          else if (!isTRUE(msg <- validObject(coercedVal, TRUE)))
            stop("Attempt to set invalid value on '", fieldName, "': ", msg)
          else val <- coercedVal
        }
        ## careful here; if field is active binding, it might not change
        oldVal <- .fieldName
        .fieldName <<- val
        if (hasPrototype) initName <<- TRUE
        if (!identical(oldVal, .fieldName)) {
          thisSignal$emit()
        }
      }
    }, list(.fieldName = as.name(.fieldName),
            fieldClass = fieldClass, fieldName = fieldName,
            thisSignal = as.name(thisSignal),
            initName = as.name(initName),
            hasPrototype = hasPrototype,
            prototype = prototype))))
  }, fields, names(fields), .fieldNames, .initNames, hasPrototype,
     prototype[names(fields)], paste(names(fields), "Changed", sep = ""))
  indSigs <- lapply(names(fields), function(nm) {
    nm <- paste(nm, "Changed", sep = "")
    fieldWithPrototype(nm, "Signal", Signal())
  })
  c(activeFields, structure(fields, names = .fieldNames),
    structure(rep("logical", sum(hasPrototype)),
              names = .initNames[hasPrototype]),
    unlist(indSigs))
}

##' A simple wrapper around \code{\link[methods]{setRefClass}} for
##' creating subclasses of \code{\linkS4class{PropertySet}}. It
##' ensures that all fields of the subclass are defined via
##' \code{\link{properties}}.
##'
##' @title Subclassing PropertySet
##' @param Class class name
##' @param fields list of fields
##' @param prototype list of default values, as in
##' \code{\link[methods]{setClass}}.
##' @param contains superclasses, one of which must extend PropertySet
##' @param ... additional arguments to \code{setRefClass}
##' @param where the environment in which to define the class
##' @return the class generator object
##' @rdname PropertySet-class
##' @author Michael Lawrence
setPropertySet <- function(Class, fields = list(), prototype = list(),
                           contains = "PropertySet", ...,
                           where = topenv(parent.frame()))
{
  if (!any(sapply(contains, extends, "PropertySet")))
    stop("At least one class in 'contains' must extend PropertySet")
  setRefClass(Class, properties(fields, prototype), contains, ...,
              where = where)
}

##' Coercion from \code{PropertySet} to \code{list}.
##'
##' This coersion only return a list of properties instances. 
##' filtering out singal function and other fields which are
##' not properties.
##' @title Coercion to \code{list}
##' @param x A \code{PropertySet} object.
##' @return A list of properties instance.
##' @author Tengfei Yin
##' @docType methods
##' @rdname as.list-methods
##' @aliases as.list
##' @aliases as.list,PropertySet-method
##' @aliases show,PropertySet-method
##' @examples
##' filt.gen <- setRefClass("Filter", properties(list(cutoff = "NonnegativeInteger",
##'                                                 weight = "PositiveInteger")),
##'                            contains = "PropertySet")
##' obj <- filt.gen$new(cutoff = NonnegativeInteger(0),
##'                     weight = PositiveInteger(1))
##' obj$properties()
##' as.list(obj)
setMethod("as.list", "PropertySet", function(x) {
  x <- as(x, "list")
  x
})

setAs("PropertySet", "list", function(from) {
  nms <- from$properties()
  lst <- lapply(names(nms), function(x){
    from$field(x)
  })
  names(lst) <- names(nms)
  lst
})

setMethod("show", "PropertySet", function(object) {
  show(as.list(object))
})

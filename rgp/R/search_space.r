## search_space.R
##   - Functions for defining the search space for Genetic Programming
##
## RGP - a GP system for R
## 2010 Oliver Flasch (oliver.flasch@fh-koeln.de)
## with contributions of Thomas Bartz-Beielstein, Olaf Mersmann and Joerg Stork
## released under the GPL v2
##

##' Functions for defining the search space for Genetic Programming
##'
##' The GP search space is defined by a set of functions, a set of
##' input variables, a set of constant constructor functions, and some
##' rules how these functions, input variables, and constants may be
##' combined to form valid symbolic expressions.  The function set is
##' simply given as a set of strings naming functions in the global
##' environment. Input variables are also given as strings.
##' Combination rules are implemented by a type system and defined by
##' assigning sTypes to functions, input variables, and constant
##' constructors.
##'
##' Function sets and input variable sets are S3 classes containing
##' the following fields: \code{\$all} contains a list of all
##' functions, or input variables, or constant factories.
##' \code{\$byRange} contains a table of all input variables, or
##' functions, or constant factories, indexed by the string label of
##' their sTypes for input variables, or by the string label of their
##' range sTypes for functions, or by the string label of their range
##' sTypes for constant factories. This field exists mainly for
##' quickly finding a function, input variable, or constant factory
##' that matches a given type.
##' 
##' Multiple function sets, or multiple input variable sets, or
##' multiple constant factory sets can be combined using the
##' \code{\link{c}} function. \code{functionSet} creates a function
##' set. \code{inputVariableSet} creates an input variable set.
##' \code{constantFactorySet} creates a constant factory set.
##'
##' Probability weight for functions, input variables, and constants
##' can be given by tagging constant names, input variables, and
##' constant factory functions via the \code{pw} function (see
##' the examples). The predicate \code{hasPw} can be used to check
##' if an object \code{x} has an associated probability weight.
##' The function \code{getPw} returns the probability weight
##' associated with an object \code{x}, if available.
##'
##' @param ... Names of functions or input variables given as strings.
##' @param list Names of functions or input variables given as a list of strings.
##' @param parentEnvironmentLevel Level of the parent environment used to resolve
##'   function names.
##' @param recursive Ignored when concatenating function- or input variable sets.
##' @param x An object (function name, input variable name, or constant
##'   factory) to tag with a probability \code{pw}.
##' @param pw A probability weight.
##' @param default A default probability weight to return iff no probability
##'   weight is associated with an object.
##' @return A function set or input variable set.
##'
##' @examples
##' # creating an untyped search space description...
##' functionSet("+", "-", "*", "/", "exp", "log", "sin", "cos", "tan")
##' inputVariableSet("x", "y")
##' constantFactorySet(function() runif(1, -1, 1))
##' # creating an untyped function set with probability weights...
##' functionSet(pw("+", 1.2), pw("-", 0.8), pw("*", 1.0), pw("/", 1.0))
##' 
##' @rdname searchSpaceDefinition
##' @export
functionSet <- function(..., list = NULL, parentEnvironmentLevel = 1) {
  ll <- if (missing(list)) list(...) else c(list, ...)
  funcset <- list()
  class(funcset) <- c("functionSet", "list")
  funcset.all <- ll
  funcset.byType <- sortByType(funcset.all)
  funcset.byRange <- sortByRange(funcset.all)
  funcset$all <- extractAttributes(funcset.all, "probabilityWeight", default = 1.0)
  funcset$byType <- Map(function(set) extractAttributes(set, "probabilityWeight", default = 1.0),
                        funcset.byType)
  funcset$byRange <- Map(function(set) extractAttributes(set, "probabilityWeight", default = 1.0),
                         funcset.byRange)
  funcset$nameStrings <- as.character(funcset$all)
  funcset$arities <- numeric()
  funcset$envir <- parent.frame(n = parentEnvironmentLevel)
  for (nameString in funcset$nameStrings) {
    funcset$arities <- c(funcset$arities, arity(get(nameString, envir = funcset$envir)))
  }
  funcset
}

##' @rdname searchSpaceDefinition
##' @export
inputVariableSet <- function(..., list = NULL) {
  ll <- if (missing(list)) list(...) else c(list, ...)
  inset <- list()
  class(inset) <- c("inputVariableSet", "list")
  inset.all <- ll
  inset.byType <- sortByType(inset.all)
  inset.byRange <- sortByRange(inset.all) # TODO remove this field
  inset$all <- extractAttributes(inset.all, "probabilityWeight", default = 1.0)
  collectFormals <- function(inputVariable) {
    if (is.character(inputVariable)) {
      list(as.symbol(inputVariable)) 
    } else {
      extractLeafSymbols(inputVariable)
    }
  }
  inset$allFormals <- unique(flatten(Map(collectFormals, inset$all)))
  inset$byType <- Map(function(set) extractAttributes(set, "probabilityWeight", default = 1.0),
                      inset.byType)
  inset$byRange <- Map(function(set) extractAttributes(set, "probabilityWeight", default = 1.0),
                       inset.byRange)
  inset$nameStrings <- as.character(inset$all)
  inset
}

##' @rdname searchSpaceDefinition
##' @export
constantFactorySet <- function(..., list = NULL) {
  ll <- if (missing(list)) list(...) else c(list, ...)
  constfacset <- list()
  class(constfacset) <- c("constantFactorySet", "list")
  constfacset.all <- ll
  constfacset.byType <- sortByType(constfacset.all)
  constfacset.byRange <- sortByRange(constfacset.all) # TODO remove this field
  constfacset$all <- extractAttributes(constfacset.all, "probabilityWeight", default = 1.0)
  constfacset$byType <- Map(function(set) extractAttributes(set, "probabilityWeight", default = 1.0),
                            constfacset.byType)
  constfacset$byRange <- Map(function(set) extractAttributes(set, "probabilityWeight", default = 1.0),
                             constfacset.byRange)
  constfacset
}

##' @rdname searchSpaceDefinition
##' @export
pw <- function(x, pw) {
  if (hasPw(x))
    stop("pw: Object ", x, " already has an probability weight (of ", pw, ").")
  attr(x, "probabilityWeight") <- pw
  x
}

##' @rdname searchSpaceDefinition
##' @export
hasPw <- function(x) !is.null(attr(x, "probabilityWeight"))

##' @rdname searchSpaceDefinition
##' @export
getPw <- function(x, default = 1.0)
  if (hasPw(x)) attr(x, "probabilityWeight") else default

##' @rdname searchSpaceDefinition
##' @method c functionSet 
##' @S3method c functionSet 
##' @export
c.functionSet <- function(..., recursive = FALSE) {
  fSets <- list(...)
  combinedFsets <- list()
  for (fSet in fSets) combinedFsets <- append(fSet$all, combinedFsets)
  functionSet(list = combinedFsets, parentEnvironmentLevel = 2)
}

##' @rdname searchSpaceDefinition
##' @method c inputVariableSet
##' @S3method c inputVariableSet
##' @export
c.inputVariableSet <- function(..., recursive = FALSE) {
  iSets <- list(...)
  combinedIsets <- list()
  for (iSet in iSets) combinedIsets <- append(iSet$all, combinedIsets)
  inputVariableSet(list = combinedIsets)
}

##' @rdname searchSpaceDefinition
##' @method c constantFactorySet
##' @S3method c constantFactorySet
##' @export
c.constantFactorySet <- function(..., recursive = FALSE) {
  cSets <- list(...)
  combinedCsets <- list()
  for (cSet in cSets) combinedCsets <- append(cSet$all, combinedCsets)
  constantFactorySet(list = combinedCsets)
}

##' Tabulate a list of functions or input variables by their sTypes
##'
##' @param x A list of functions or input variables to sort by sType.
##' @return A table of the objects keyed by their sTypes.
sortByType <- function(x) {
  byTypeTable <- list()
  for (o in x) {
    o <- if (is.character(o)) as.name(o) else o
    if (hasStype(o)) {
      oStype <- sType(o)
      if (is.null(byTypeTable[[oStype$string]])) byTypeTable[[oStype$string]] <- list()
      byTypeTable[[oStype$string]] <- append(byTypeTable[[oStype$string]], list(o))
    }
  }
  byTypeTable
}

##' Tabulate a list of functions or input variables by the range part of their sTypes
##'
##' @param x A list of functions or input variables to sort by range sType.
##' @return A table of the objects keyed by their range sTypes.
sortByRange <- function(x) {
  byRangeTable <- list()
  for (o in x) {
    o <- if (is.character(o)) as.name(o) else o
    if (hasStype(o)) {
      oStype <- sType(o)
      oStypeRange <- if (inherits(oStype, "sFunctionType")) oStype$range else oStype
      if (is.null(byRangeTable[[oStypeRange$string]])) byRangeTable[[oStypeRange$string]] <- list()
      byRangeTable[[oStypeRange$string]] <- append(byRangeTable[[oStypeRange$string]], list(o))
    }
  }
  byRangeTable
}

##' Extract a given attribute of all objects in a list and tag that list with the
##' list of extracted attributes
##'
##' @param x A list with objects containing the attribute \code{attribute}.
##' @param extractAttribute The attribute to extract from all objects in the list \code{x}.
##' @param tagAttribute The name of the attribute for \code{x} holding the list of
##'   extracted attributes.
##' @param default A default value to return if an object in \code{x} has no attribute
##'   \code{attribute}.
##' @return The list \code{x}, tagged with a new attribute \code{tagAttribute}.
extractAttributes <- function(x, extractAttribute, tagAttribute = extractAttribute, default = NULL) {
  extractedAttributes <- Map(function(o) if (!is.null(attr(o, extractAttribute)))
                             attr(o, extractAttribute) else default, x)
  attr(x, tagAttribute) <- extractedAttributes
  x
}

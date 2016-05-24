##twList.r
##2013-09-60 dmontaner@cipf.es
##TiddlyWikiR library

##' @name twList-class
##' @docType class
##' 
##' @author David Montaner \email{dmontaner@@cipf.es}
##'
##' @aliases twList
##' @aliases initialize,twList-method
##' @aliases elements elements<-
##' @aliases level level<-
##' @aliases type type<-
##' 
##' @keywords wiki list
##' @seealso \code{\link{twLink}} and \code{\link{twImage}} \code{\link{twTable}}
##' 
##' @title A class to handle TiddlyWiki lists.
##'
##' @description Some utilities to insert ordered and unordered lists into TiddlyWiki.
##'
##' @section Usage:
##' new ("twList", elements, level, type)
##' 
##' twList (elements, ...)
##' 
##' @section Slots: 
##'  \describe{
##'    \item{\code{elements}:}{a character vector of the elements of the list.}
##'    \item{\code{level}:}{a numeric vector indicating the level of indentation of each of the elements of the list.}
##'    \item{\code{type}:}{a character vector indicating the desired bullet type. Allowed values are: "o" for "ordered" type elements and "u" for unordered. }
##'  }
##'
##' @examples
##' list1 <- new ("twList", elements = c("line one", "line two", "line three"),
##'               level = c(1,2,2), type = c("u", "o", "o"))
##' list2 <- twList (LETTERS[1:5])
##'
##' wikify (list1)
##' wikify (list2)
##' 
##' @export

setClass ("twList",
          representation = representation (
            elements = "character",
            level = "numeric", #"integer",
            type = "character"
            ))

##' @export
twList <- function (...) {
  new ("twList", ...)
}

################################################################################

##' @export
setMethod ("initialize", "twList",
           function (.Object,
                     elements,
                     level,
                     type) { ##the variable has to be called .Object
             
             if (missing (elements)) {
               stop ("elements is missing with no default.")
             } else {
               .Object@elements <- elements
             }
             
             ## ################################################################
             
             if (missing (level)) {
               .Object@level <- rep (1L, times = length (elements))
             } else {
               if (length (level) != length (elements)) {
                 stop ('"level" and "elements" have different dimensions.')
               } else {
                 .Object@level <- level
               }
             }

             ## ################################################################
             
             if (missing (type)) {
               .Object@type <- rep ("u", times = length (elements))
             } else {
               if (! all (type %in% c("u", "o"))) {
                 stop ('"type" values must be one of "u", "o".')
               }
               if (length (type) != length (elements)) {
                 if (length (elements) %% length (type) == 0) {
                   .Object@type <- rep (type, times = length (elements) / length (type))
                 } else {
                   warning ("elements length is not a multiple of type length")
                   .Object@type <- rep (type, times = ceiling (length (elements) / length (type)))[1:length (elements)]
                 }
               } else {
                 .Object@type <- type
               }
             }
             
             ## ################################################################             
             
             return (.Object)
           })

################################################################################

## does not need @export
setValidity ("twList", function (object) { ##the variable can be named other than "object"
  out <- NULL
  
  if (length (object@elements) != length (object@type)) {
    out <- c (out, '"elements" and "type" have different length.')
  }
  
  if (length (object@elements) != length (object@level)) {
    out <- c (out, '"elements" and "level" have different length.')
  }
  
  if (! all (object@type %in% c("u", "o"))) {
    out <- c (out, '"type" values must be one of "u", "o".')
  }
  
  ##RETURN
  if (is.null (out)) {
    out <- TRUE
  }
  return (out)
})


################################################################################
### ACCESSOR and REPLACEMENTS
################################################################################

##' @export
elements <- function (x) {
  x@elements
}

##' @export
'elements<-' <- function (x, value) {
  x@elements <- value
  return (x)
}

########################################

##' @export
level <- function (x) {
  x@level
}

##' @export
'level<-' <- function (x, value) {
  x@level <- value
  return (x)
}

########################################

##' @export
type <- function (x) {
  x@type
}

##' @export
'type<-' <- function (x, value) {
  x@type <- value
  return (x)
}

################################################################################

###Wikify

wikify.list <- function (object) {
  
  valid.o <- c ("o") # "ordered"     ##valid tags for ordered and unordered elements
  valid.u <- c ("u") # "unordered"
  
  ##FORMAT DATA
  elements <- object@elements
  level    <- object@level
  type     <- object@type
  
  bullet <- rep ("", times = length (elements))
  
  bullet[tolower (type) %in% valid.o] <- "#"
  bullet[tolower (type) %in% valid.u] <- "*"
  
  for (i in 1:length(elements)) {
    bullet[i] <- paste (rep (bullet[i], times = level[i]), collapse = "")
  }

  ##RETURN
  out <- paste (bullet, elements, sep = " ")
  return (out)
}

###Wikify method
setMethod ("wikify", "twList", wikify.list)

################################################################################

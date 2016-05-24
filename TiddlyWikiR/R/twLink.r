##twLink.r
##2013-06-13 dmontaner@cipf.es

##' @name twLink-class
##' @docType class
##' @author David Montaner \email{dmontaner@@cipf.es}
##' 
##' @aliases twLink
##' @aliases initialize,twLink-method
##' @aliases label,twLink-method  label<-,twLink-method
##' @aliases ref,twLink-method  ref<-,twLink-method
##' @aliases show,twLink-method
##' 
##' @keywords wiki link
##' @seealso \code{\link{twTable}} and \code{\link{twImage}}
##' 
##' @title A class to handle TiddlyWiki links.
##'
##' @description Some utilities to insert internal or external links into TiddlyWiki.
##'
##' @section Usage:
##' new ("twLink", label, ref)
##' 
##' twLink ("the text here", ref = "the target here")
##' 
##' @section Slots: 
##'  \describe{
##'    \item{\code{label}:}{the visible text or anchor.}
##'    \item{\code{ref}:}{the target or URL to be redirected to.}
##'  }
##'
##' @section Methods:
##'  \describe{
##'    \item{\code{label} and \code{label (object) <- value}:}{Access and set \code{label}.}
##'    \item{\code{ref} and \code{ref (object) <- value}:}{Access and set \code{ref}.}
##' }
##'
##' @examples
##' new ("twLink", label = "the text here", ref = "aTiddlerNameHere")
##' twLink (label = "the text here", ref = "http://www.r-project.org/")
##'
##' twLink ("here", ref = "http://www.dmontaner.com")
##' paste ("see my home page",
##'        twLink ("here", ref = "http://www.dmontaner.com"), ", I hope you like it.")
##' 
##' @export
setClass ("twLink",
          representation = representation (
            label = "character",   ##could be a vector of any type?
            ref   = "character"))

##' @export
twLink <- function (...) {
  new ("twLink", ...)
}

################################################################################

##' @export
setMethod ("initialize", "twLink",
           function (.Object, label, ref) { ##the variable has to be called .Object

             if (missing (label)) {
               stop ("label is missing with no default.")
             } else {
               .Object@label <- as.character (label)
             }
             
             if (missing (ref)) {
               .Object@ref <- NA_character_
             } else {
               .Object@ref <- as.character (ref)
             }
             
             return (.Object)
           })

################################################################################

setValidity ("twLink", function (object) {
  separador <- "\n   "
  out <- ""
  ##
  if (length (object@label) != length (object@ref)) {
    out <- paste (out, "label and ref should have the same length.", sep = separador)
  }
  ##
  if (out == "") out <- TRUE
  return (out)
})


################################################################################
### ACCESSOR and REPLACEMENTS
################################################################################

setMethod ("label", "twLink", function (object) {
  return (object@label)
})

setReplaceMethod ("label", "twLink", function (object, value) {
  object@label <- value
  return (object)
})

########################################

setMethod ("ref", "twLink", function (object) {
  return (object@ref)
})
##
setReplaceMethod ("ref", "twLink", function (object, value) {
  object@ref <- value
  return (object)
})


################################################################################

setMethod ("show", "twLink", function (object) {  ##the variable has to be called object as in the generic; see: getGeneric ("show")
  out <- paste (label (object), "===>", ref (object))
  cat (out, sep = "\n", fill = TRUE)
})

################################################################################

setMethod ("wikify", "twLink", function (object, leStr = "[[", miStr = "|", riStr = "]]") {
  if (is.na (object@ref)) {
    out <- paste (leStr, object@label, riStr, sep = "")
  } else {
    out <- paste (leStr, object@label, miStr, object@ref, riStr, sep = "")
  }
  return (out)
})

################################################################################

setMethod ("as.character", "twLink", function (x) {
  wikify (x)
})

################################################################################

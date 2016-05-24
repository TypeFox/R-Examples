##tiddler.r
##2013-06-13 dmontaner@cipf.es
##TiddlyWikiR library

## Tiddler header structure:
## <div title="New Tiddler" creator="YourName" modifier="YourName"
##  created="201306131256" modified="201306131631"
## tags="tag1 [[tag 2]] tag3" changecount="2">

##' @name tiddler-class
##' @docType class
##' @author David Montaner \email{dmontaner@@cipf.es}
##' 
##' @aliases tiddler
##' @aliases initialize,tiddler-method
##' 
##' @aliases tdTitle tdTitle<- 
##' @aliases tdContent tdContent<-
##' @aliases tdTags tdTags<-
##' @aliases tdCreator tdCreator<-
##' @aliases tdModifier tdModifier<-
##' @aliases tdCreated tdCreated<-
##' @aliases tdModified tdModified<-
##' @aliases tdChangecount tdChangecount<-
##'
##' @aliases newTiddler
##' 
##' @keywords wiki tiddler
##' 
##' @title A class to handle TiddlyWiki tiddlers.
##' 
##' @description Object to store a TiddlyWiki tiddler information.
##'
##' @section Usage:
##' new ("tiddler", title, content, ...)
##' 
##' newTiddler (title, content, ...)
##'
##' @section Slots:
##' \describe{
##'   \item{\code{title}:}{title of the tiddler.}
##'   \item{\code{content}:}{content of the tiddler.}
##'   \item{\code{tags}:}{tags of the tiddler.}
##'   \item{\code{creator}:}{indicates who created the tiddler.}
##'   \item{\code{modifier}:}{indicates who modified the tiddler.}
##'   \item{\code{created}:}{timestamp of the tiddler creation.}
##'   \item{\code{modified}:}{timestamp of the tiddler modification.}
##'   \item{\code{changecount}:}{number of times the tiddler was changed.}
##'  }
##'
##' @section Methods:
##' Methods are implemented to access and set all slots form a tiddler object.
##' The most useful ones may be:
##  HERE Y WOULD RATHER USING describe INSTEAD OF itemize, BUT SOME WARNINGS ARE RETURNED BY oxygenize
##' \itemize{                        
##' \item{\code{tdTitle}}
##' \item{\code{tdContent}}
##' \item{\code{tdTags}}
##' \item{\code{tdCreator}}
##' }
##' But \code{tdModifier}, \code{tdCreated}, \code{tdModified} and \code{tdChangecount}
##' are also implemented.
##' 
##' @examples
##' tiddler1 <- new ("tiddler", title = "MyNewTiddler1", content = "the first content")
##' tiddler2 <- newTiddler (title = "MyNewTiddler2", content = "the second content", creator = "me")

##' @export
setClass ("tiddler",
          representation = representation (
            title       = "character",
            content     = "list",       ## the only slot named by me
            tags        = "character",  ## may be longer than 1
            creator     = "character",
            modifier    = "character",
            created     = "character",
            modified    = "character",
            changecount = "numeric"))

##' @export
newTiddler <- function (...) {
  new ("tiddler", ...)
}

################################################################################

##' @export
setMethod ("initialize", "tiddler",
           function (.Object, title, content, tags, creator, modifier, created, modified, changecount) { ##the variable has to be called .Object
             
             if (missing (title)) {
               stop ("title is missing with no default.")
             } else {
               .Object@title <- title
             }
             
             if (!missing (content)) {
               if (is.list (content)){
                 .Object@content <- content
               } else {
                 .Object@content <- list (content)
               }
             }
             
             if (!missing (tags)){
               if (is.character (tags)) {
                 .Object@tags <- tags
               } else {
                 .Object@tags <- as.character (tags)
               }
             }
             
             if (missing (creator)) {
               .Object@creator <- "TiddlyWikiR"
             } else {
               .Object@creator <- creator
             }
             
             if (missing (modifier)) {
               .Object@modifier <- .Object@creator
             } else {
               .Object@modifier <- modifier
             }
             
             if (missing (created)) {
               .Object@created <- format (Sys.time(), "%Y%m%d%H%M")
             } else {
               .Object@created <- created
             }
             
             if (!missing (modified)) {
               .Object@modified <- modified
             }
             if (missing (changecount)) {
               .Object@changecount <- 1
             } else {
               .Object@changecount <- changecount
             }

             return (.Object)
           })

################################################################################

## does not need @export
setValidity ("tiddler", function (object) {
  separador <- "\n"
  out <- ""
  if (length (object@title) == 0) {
    out <- paste (out, "Any tiddler must have a title.", sep = separador)
  }
  if (length (object@creator) == 0) {
    out <- paste (out, "Any tiddler must have a creator.", sep = separador)
  }
  if (length (object@created) == 0) {
    out <- paste (out, "Any tiddler must have a creation date; fill the created slot of the tiddler.", sep = separador)
  }
  if (length (object@changecount) == 0) {
    out <- paste (out, "Any tiddler must have a number of times that has been accessed; fill the changecount slot of the tiddler with a value of 1.", sep = separador)
  }
  ##
  if (out == "") out <- TRUE
  return (out)
})


################################################################################
### ACCESSOR and REPLACEMENTS
################################################################################

##' @export
tdTitle <- function (x) {
  x@title
}

##' @export
'tdTitle<-' <- function (x, value) {
  x@title <- value
  return (x)
}

########################################

##' @export
tdContent <- function (x) {
  x@content
}

##' @export
'tdContent<-' <- function (x, value) {
  if (!is.list (value)){
    value <- list (value)
  }
  x@content <- value
  return (x)
}

########################################

##' @export
tdTags <- function (x) {
  x@tags
}

##' @export
'tdTags<-' <- function (x, value) {
  x@tags <- value
  return (x)
}

########################################
##' @export
tdCreator <- function (x) {
  x@creator
}

##' @export
'tdCreator<-' <- function (x, value) {
  x@creator <- value
  return (x)
}

########################################
##' @export
tdModifier <- function (x) {
  x@modifier
}

##' @export
'tdModifier<-' <- function (x, value) {
  x@modifier <- value
  return (x)
}

########################################

##' @export
tdCreated <- function (x) {
  x@created
}

##' @export
'tdCreated<-' <- function (x, value) {
  x@created <- value
  return (x)
}

########################################

##' @export
tdModified <- function (x) {
  x@modified
}

##' @export
'tdModified<-' <- function (x, value) {
  x@modified <- value
  return (x)
}

########################################

##' @export
tdChangecount <- function (x) {
  x@changecount
}

##' @export
'tdChangecount<-' <- function (x, value) {
  x@changecount <- value
  return (x)
}

########################################

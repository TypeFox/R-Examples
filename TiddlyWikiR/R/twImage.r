##twImage.r
##2013-06-20 dmontaner@cipf.es
##TiddlyWikiR library

##' @name twImage-class
##' @docType class
##' @author David Montaner \email{dmontaner@@cipf.es}
##' 
##' @aliases twImage
##' @aliases initialize,twImage-method
##' @aliases label,twImage-method  label<-,twImage-method
##' @aliases ref,twImage-method  ref<-,twImage-method
##' @aliases imgf imgf<-
##' @aliases align,twImage-method  align<-,twImage-method
##' @aliases width width<-
##' @aliases height height<-
##' 
##' @keywords wiki image
##' @seealso \code{\link{twTable}} and \code{\link{twLink}}
##' 
##' @title A class to handle TiddlyWiki images.
##'
##' @description Some utilities to insert images into TiddlyWiki.
##'
##' @details Resizing is done using http://www.tiddlytools.com/#ImageSizePlugin
##' where (w,h) indicates the desired width and height (in CSS units, e.g., px, em, cm, in, or \%).
##' [img(w+,h+)[...][...]]
##'
##' @section Usage:
##' new ("twImage", imgf, ...)
##' 
##' twImage (imgf, ...)
##' 
##' @section Slots: 
##'  \describe{
##'    \item{\code{imgf}:}{the image file.}
##'    \item{\code{label}:}{label of the image.}
##'    \item{\code{ref}:}{optional; a target tiddler or URL to be redirected to.}
##'    \item{\code{align}:}{alignment of the image. May be NA, "r" or "l".}
##'    \item{\code{width}:}{width scaling factor.}
##'    \item{\code{height}:}{height scaling factor.}
##'  }
##'
##' @section Methods:
##'  \describe{
##'    \item{\code{imgf}:}{access the imgf slot.}
##'    \item{\code{label}:}{access the label slot.}
##'    \item{\code{ref}:}{access the ref slot.}
##'    \item{\code{align}:}{access the align slot.}
##'    \item{\code{width}:}{access the width slot.}
##'    \item{\code{height}:}{access the height slot.}
##' }
##'
##'
##' @examples
##' new ("twImage", imgf = "myImageFile.png")
##' twImage ("myImageFile.png")
##' wikify (twImage ("myImageFile.png"))
##' 
##' @export

##newWikiImage <- setClass ("wikiImage",
setClass ("twImage",
          representation = representation (
            imgf   = "character",
            label  = "character",
            ref    = "character",
            align  = "character",
            width  = "character",
            height = "character"))

##' @export
twImage <- function (...) {
  new ("twImage", ...)
}

################################################################################

##' @export
setMethod ("initialize", "twImage",
           function (.Object, imgf, label, ref, align, width, height) { ##the variable has to be called .Object

             if (missing (imgf)) {
               stop ("imgf is missing with no default.")
             } else {
               .Object@imgf <- imgf
             }

             if (missing (label)) {
               .Object@label <- NA_character_
             } else {
               .Object@label <- label
             }

             if (missing (ref)) {
               .Object@ref <- NA_character_
             } else {
               .Object@ref <- ref
             }

             if (missing (align)) {
               .Object@align <- NA_character_
             } else {
               .Object@align <- align
             }
             
             if (missing (width)) {
               .Object@width <- NA_character_
             } else {
               .Object@width <- width
             }
             
             if (missing (height)) {
               .Object@height <- NA_character_
             } else {
               .Object@height <- height
             }
             
             return (.Object)
           })

################################################################################

## does not need @export
setValidity ("twImage", function (object) { ##the variable can be named other than "object"
  separador <- "\n   "
  out <- ""
  ##
  if (!object@align %in% c(NA, "l", "r")) {
    out <- paste (out, 'align should be one of NA, "l", "r"', sep = separador)
  }
  ##
  if (out == "") out <- TRUE
  return (out)
})


################################################################################
### ACCESSOR and REPLACEMENTS
################################################################################

setMethod ("label", "twImage", function (object) {
  return (object@label)
})

setReplaceMethod ("label", "twImage", function (object, value) {
  object@label <- value
  return (object)
})

########################################

setMethod ("ref", "twImage", function (object) {
  return (object@ref)
})
##
setReplaceMethod ("ref", "twImage", function (object, value) {
  object@ref <- value
  return (object)
})

########################################

##' @export
imgf <- function (x) {
  x@imgf
}

##' @export
'imgf<-' <- function (x, value) {
  x@imgf <- value
  return (x)
}

########################################

setMethod ("align", "twImage", function (object) {
  return (object@align)
})

setReplaceMethod ("align", "twImage", function (object, value) {
  object@align <- value
  return (object)
})

########################################

##' @export
width <- function (x) {
  x@width
}

##' @export
'width<-' <- function (x, value) {
  x@width <- value
  return (x)
}

########################################

##' @export
height <- function (x) {
  x@height
}

##' @export
'height<-' <- function (x, value) {
  x@height <- value
  return (x)
}

################################################################################

setMethod ("show", "twImage", function (object) {  ##the variable has to be called object as in the generic; see: getGeneric ("show")
  out <- imgf (object)
  out <- c(out, paste (label (object), "===>", ref (object)))

  if (!is.na (align (object))) {
    out <- c(out, paste ("align:", align (object)))
  }

  if (!is.na (width (object)) | !is.na (width (object))) {
    out <- c(out, paste ("resize:", "(", width (object), ",", height (object) ,")"))
  }

  cat (out, sep = "\n", fill = TRUE)
})

################################################################################

setMethod ("wikify", "twImage", function (object) {
  
  if (is.na (label (object))) {
    out <- paste ("[", imgf (object), "]", sep = "")
  } else {
    out <- paste ("[", label (object), "|", imgf (object), "]", sep = "")
  }
  
  if (is.na (ref (object))) {
    out <- paste (out, "]", sep = "")
  } else {
    out <- paste (out, "[", ref (object), "]]", sep = "")
  }
  
  al <- ""
  if (align (object) %in% c ("l", "left")) {
    ##al <- "<" ##not working???
    al <- ""
  }
  if (align (object) %in% c ("r", "right")) {
    al <- ">"
  }
  
  if (is.na (width (object)) & is.na (height (object))) {
    sc <- ""
  } else {
    sc <- "("
    if (!is.na (width (object))) sc <- paste (sc, width (object), sep = "")
    sc <- paste (sc, ",", sep = "")
    if (!is.na (height (object))) sc <- paste (sc, height (object), sep = "")
    sc <- paste (sc, ")", sep = "")
  }

  out <- paste ("[", al, "img", sc, out, sep = "")
  return (out)
})

################################################################################

setMethod ("as.character", "twImage", function (x) {
  wikify (x)
})

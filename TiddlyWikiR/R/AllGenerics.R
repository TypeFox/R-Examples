##AllGenerics.R
##2013-07-01 dmontaner@cipf.es
##TiddlyWikiR library: Generic functions.

## It seems that the file has to be called exactly AllGenerics.R
## (including the capitalized extension).
## It may be due to the "collate" description
## but I do not quiet understand why.


##The 'title' function exists in the 'graphics' package.
##Set it as generic to be used in twTable.
##Needs to be exported.
##' @export
setGeneric ("title") 

################################################################################

##' @name label
##' @docType methods
##' @author David Montaner \email{dmontaner@@cipf.es}
##' 
##' @aliases label label<-
##' 
##' @keywords label wiki image link
##' @seealso \code{\link{twLink}}, \code{\link{twImage}}
##' 
##' @title Extract and replace the "label" form \code{twLink} and \code{twImage} objects.
##' 
##' @description Generic function used to access the "label" slot in \code{twLink} and \code{twImage} objects.

##' @export
setGeneric ("label",   function (object)        standardGeneric ("label"))
##' @export
setGeneric ("label<-", function (object, value) standardGeneric ("label<-"))

################################################################################

##' @name ref
##' @docType methods
##' @author David Montaner \email{dmontaner@@cipf.es}
##' 
##' @aliases ref ref<-
##' 
##' @keywords ref wiki image link
##' @seealso \code{\link{twLink}}, \code{\link{twImage}}
##' 
##' @title Extract and replace the "ref" form \code{twLink} and \code{twImage} objects.
##'
##' @description Generic function used to access the "ref" slot in \code{twLink} and \code{twImage} objects.

##' @export
setGeneric ("ref",   function (object)        standardGeneric ("ref"))
##' @export
setGeneric ("ref<-", function (object, value) standardGeneric ("ref<-"))

################################################################################

##' @name align
##' @docType methods
##' @author David Montaner \email{dmontaner@@cipf.es}
##' 
##' @aliases align align<-
##' 
##' @keywords align wiki image table
##' @seealso \code{\link{twTable}}, \code{\link{twImage}}
##' 
##' @title Extract and replace the "alignment" form \code{twTable} and \code{twImage} objects.
##'
##' @description Generic function used to access the "align" slot in \code{twTable} and \code{twImage} objects.

##' @export
setGeneric ("align",   function (object)        standardGeneric ("align"))
##' @export
setGeneric ("align<-", function (object, value) standardGeneric ("align<-"))

################################################################################

##' @name wikify
##' @docType methods
##' @author David Montaner \email{dmontaner@@cipf.es}
##' 
##' @aliases wikify
##' @aliases wikify,character-method
##' @aliases wikify,tiddler-method
##' @aliases wikify,twImage-method
##' @aliases wikify,twLink-method
##' @aliases wikify,twList-method
##' @aliases wikify,twTable-method
##' 
##' @aliases as.character,twLink-method
##' @aliases as.character,twImage-method
##' 
##' @keywords wiki convert paste
##' @seealso \code{\link{writeTiddlers}}, \code{\link{writeTags}}, \code{capture.output}
##' 
##' @title Format R objects to be inserted into TiddlyWiki.
##'
##' @description Generic function used to do the final formatting of R objects
##' so that they can be introduced into a TiddlyWiki template file.
##'
##' @details At the moment, specific methods are described for classes:
##' \code{tiddler}, \code{twImage}, \code{twLInk}, \code{twList} and \code{twTable}.
##' A character vector will be inserted "as is" into the TiddlyWiki file;
##' each element of the vector will be a line within the file.
##' Any other R object will be introduced as "code" chunks into the file
##' using the standard R results display.
##'
##' The function is not intended to be called directly
##' but from within the functions \code{writeTiddlers} and \code{writeTags}.
##'
##' @return A character vector. Each element is intended to be a line within the TiddlyWiki template file.

##' @export
#setGeneric ("wikify", function (object, ...) standardGeneric ("wikify")) #no default method
setGeneric ("wikify", function (object, ...) {
  out <- capture.output (object) ##default display in R
  out <- c("{{{", out, "}}}")
})

setMethod ("wikify", "character", function (object) { ##I am not sure if this is the standard place to put this.
  return (object)
})

################################################################################

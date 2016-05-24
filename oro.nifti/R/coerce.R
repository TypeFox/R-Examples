############################################################################
## as("anlz", "nifti")
############################################################################
#' @title Force an Object to Belong to the ANALYZE or NIfTI Class
#' 
#' @description Methods for function \code{coerce} in Package \sQuote{methods}.
#' 
#' @name coerce-methods
#' @aliases coerce-methods coerce,array,anlz-method coerce<-,array,anlz-method
#' as,array,anlz-method as<-,array,anlz-method as<-,array,nifti-method
#' coerce,list,anlz-method coerce<-,list,anlz-method coerce,array,nifti-method
#' as,array,nifti-method coerce,list,nifti-method coerce,anlz,nifti-method
#' coerce<-,array,nifti-method as<-,array,nifti-method
#' coerce<-,list,nifti-method coerce<-,anlz,nifti-method
#' @docType methods
#' @param object is an object of class \code{array} or inherits from
#' \code{array}.
#' @param Class is the name of the class to which \sQuote{object} should be
#' coerced; i.e., \code{nifti}.
#' @param value is the values used to modify \sQuote{object} (see the
#' discussion below).  You should supply an object with class \code{nifti} in
#' order to pass NIfTI header information.
#' @aliases as.anlz as.nifti
#' @param from is the object to be converted.
#' @param value is the \code{nifti} class object to use as a template for
#' various ANALYZE/NIfTI header information.
#' @param verbose is a logical variable (default = \code{FALSE}) that allows
#' text-based feedback during execution of the function.
#' @return An object of class \code{anlz} or \code{nifti}.
#' @section Methods: 
#' \describe{ 
#' \item{from = "anlz", to = "nifti"}{An object of class \code{anlz} is coerced 
#' into a NIfTI object.} 
#' \item{from = "array", to = "anlz"}{An object of class \code{array} is 
#' coerced into an ANALYZE object.} 
#' \item{from = "array", to = "nifti"}{An object of class \code{array} is 
#' coerced into a NIfTI object.} 
#' \item{from = "list", to = "anlz"}{All objects of class \code{array} in the 
#' list are coerced into ANALYZE objects.  All other objects are left alone.  
#' The original list structure is retained.}
#' \item{from = "list", to = "nifti"}{All objects of class \code{array} in the
#' list are coerced into NIfTI objects.  All other objects are left alone.  The
#' original list structure is retained.} 
#' }
#' @author 
#' Andrew Thornton \email{zeripath@@users.sourceforge.net},\cr
#' Brandon Whitcher \email{bwhitcher@@gmail.com}
#' @seealso \code{\link{as}}
#' @keywords methods
#' @rdname coerce-methods
setAs("anlz", "nifti",
      function(from) { as.nifti(from) },
      function(from, value) { as.nifti(from, value) } )

############################################################################
## as("array", "nifti")
############################################################################
#' @rdname coerce-methods
#' @name coerce-methods
setAs("array", "nifti",
      function(from) { as.nifti(from) },
      function(from, value) { as.nifti(from, value) } )

############################################################################
## as("list", "nifti")
############################################################################
#' @rdname coerce-methods
#' @name coerce-methods
setAs("list", "nifti",
      function(from) { as.nifti(from) },
      function(from, value) { as.nifti(from, value) } )

############################################################################
## as("array", "anlz")
############################################################################
#' @rdname coerce-methods
#' @name coerce-methods
setAs("array", "anlz",
      function(from) { as.anlz(from) },
      function(from, value) { as.anlz(from, value) } )

############################################################################
## as("list", "anlz")
############################################################################
#' @rdname coerce-methods
#' @name coerce-methods
setAs("list", "anlz",
      function(from) { as.anlz(from) },
      function(from, value) { as.anlz(from, value) } )

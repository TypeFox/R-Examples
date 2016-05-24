##' Dimnames for hyperSpec objects
##'
##' hyperSpec objects can have row- and column names like data.frames. The "names" of the wavelengths
##' are treated separately: see \code{\link{wl}}
##' @param x the hyperSpec object
##' @aliases dimnames
##' @keywords methods
##' @rdname dimnames
##' @docType methods
##' @author C. Beleites
##' @seealso \code{\link{wl}} for the wavelength dimension
##'
##' \code{\link[base]{dimnames}}
##' @export
##' @examples
##' dimnames (flu)
setMethod ("dimnames", signature = signature (x = "hyperSpec"), function (x){
  validObject (x)

  list (row = rownames (x@data), data = colnames (x@data),
        wl = colnames (x@data$spc))
})

##'
##'
##' @rdname dimnames
##' @aliases rownames
##' @param do.NULL handed to \code{\link[base]{rownames}} or \code{\link[base]{colnames}}: logical.
##' Should this create names if they are \code{NULL}?
##' @param prefix handed to \code{\link[base]{rownames}} or \code{\link[base]{colnames}}
##' @seealso \code{\link[base]{rownames}}
##' @export
##' @examples
##' rownames (flu)
setMethod ("rownames", signature = signature (x = "hyperSpec"), function (x, do.NULL = TRUE, prefix = "row"){
  validObject (x)
  
  rownames (x@data, do.NULL = do.NULL, prefix = prefix)
})

##' @param value the new names
##' @usage
##' \S4method{rownames}{hyperSpec} (x) <- value
##' @aliases rownames<-,hyperSpec-method
##' @rdname dimnames
##' @name rownames<-
##' @export "rownames<-"
setReplaceMethod ("rownames", signature = signature (x = "hyperSpec"), function (x, value){
  validObject (x)
  
  rownames (x@data) <- value
  x
})

##' @rdname dimnames
##' @aliases colnames
##' @seealso \code{\link[base]{colnames}}
##' @export
##' @examples
##' colnames (chondro)

setMethod ("colnames", signature = signature (x = "hyperSpec"),
           function (x, do.NULL = TRUE, prefix = "col"){
  validObject (x)
  colnames (x@data, do.NULL = do.NULL, prefix = prefix)
})

##' @rdname dimnames
##' @usage
##' \S4method{colnames}{hyperSpec} (x) <- value
##' @aliases colnames<-,hyperSpec-method
##' @name colnames<-
##' @export "colnames<-"
setReplaceMethod ("colnames", signature = signature (x = "hyperSpec"),
                  function (x, value){
  validObject (x)

  names (x@label [colnames (x@data)]) <- value
  colnames (x@data) <- value
  
  validObject (x)                       # necessary: $spc could be renamed!
  x
})

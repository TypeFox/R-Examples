## Overloading binary operators for anlz object
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @title Operations for Objects in the ANALYZE and NIfTI classes
#' @name anlz-nifti-ops
#' @rdname anlz-nifti-ops
#' @aliases Ops,anlz,anlz-method
#' @param e1 object
#' @param e2 object
#' @author John Muschellli \email{muschellij2@@gmail.com}
#' @examples
#' 
#' img01 <- anlz(array(1:64, c(4,4,4,1)), datatype=4)
#' img02 <- anlz(array(64:1, c(4,4,4,1)), datatype=4)
#' is.anlz(img01 + img02)
#' is.anlz(sqrt(2) * img01)
#' is.anlz(img02 / pi)
#' 
setMethod("Ops", signature(e1="anlz", e2="anlz"),
          function(e1, e2) {
            # e1 = drop_img_dim(e1)            
            # e2 = drop_img_dim(e2)            
            e1@.Data <- callGeneric(e1@.Data, e2@.Data)
            e1 <- resetSlopeIntercept(e1)
            e1 <- calibrateImage(e1, infok=TRUE)
            # e1 = drop_img_dim(e1, warn=FALSE)
            ## creating the datatype that is the maximal precision
            new.dtype <- max(datatype(e1), datatype(e2))
            datatype(e1) <- new.dtype
            bitpix(e1) <- convert.bitpix.anlz()[[convert.datatype.anlz(new.dtype)]]
            # validObject(e1)
            return(e1)
          }
)
#' @rdname anlz-nifti-ops
#' @aliases Ops,anlz,numeric-method
setMethod("Ops", signature(e1="anlz", e2="numeric"),
          function(e1, e2) {
            # e1 = drop_img_dim(e1)            
            e1@.Data <- callGeneric(e1@.Data, e2)
            e1 <- resetSlopeIntercept(e1)            
            e1 <- calibrateImage(e1, infok=TRUE)
            # e1 = drop_img_dim(e1)
            # validObject(e1)
            return(e1)
          }
)
#' @rdname anlz-nifti-ops
#' @aliases Ops,numeric,anlz-method
setMethod("Ops", signature(e1="numeric", e2="anlz"),
          function(e1, e2) {
            # e2 = drop_img_dim(e2)
            e2@.Data <- callGeneric(e1, e2@.Data)
            # reset to keep the same code (fewer copy/paste errors)
            e1 <- e2
            e1 <- resetSlopeIntercept(e1)            
            e1 <- calibrateImage(e1, infok=TRUE)
            # e1 = drop_img_dim(e1)
            # validObject(e1)
            return(e1)
          }
)

#' @aliases Ops,nifti,anlz-method
#' @rdname anlz-nifti-ops
setMethod("Ops", signature(e1="nifti", e2="anlz"),
          function(e1, e2) {
            # e1 = drop_img_dim(e1) 
            # e2 = drop_img_dim(e2)            
            e1@.Data <- callGeneric(e1@.Data, e2@.Data)
            e1 <- resetSlopeIntercept(e1)
            e1 <- calibrateImage(e1, infok=TRUE)
            # e1 = drop_img_dim(e1, warn=FALSE)
            ## creating the datatype that is the maximal precision
            # new.dtype = max(datatype(e1), datatype(e2))
            # datatype(e1) = new.dtype
            # bitpix(e1) = convert.bitpix.anlz()[[new.dtype]]
            # validObject(e1)
            return(e1)
          }
)
#' @rdname anlz-nifti-ops
#' @aliases Ops,anlz,nifti-method
setMethod("Ops", signature(e1="anlz", e2="nifti"),
          function(e1, e2) {
            # e1 = drop_img_dim(e1)            
            # e2 = drop_img_dim(e2)            
            e2@.Data <- callGeneric(e1@.Data, e2@.Data)
            e1 <- e2
            e1 <- resetSlopeIntercept(e1)
            e1 <- calibrateImage(e1, infok=TRUE)
            # e1 = drop_img_dim(e1, warn=FALSE)
            ## creating the datatype that is the maximal precision
            # new.dtype = max(datatype(e1), datatype(e2))
            # datatype(e1) = new.dtype
            # bitpix(e1) = convert.bitpix.anlz()[[new.dtype]]
            # validObject(e1)
            return(e1)
          }
)

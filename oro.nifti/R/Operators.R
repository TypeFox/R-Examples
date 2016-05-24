## Overloading binary operators
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @title Operations for NIfTI Objects
#' @name nifti-operators
#' @rdname niftiops
#' @aliases Ops,nifti,nifti-method
#' @param e1 is an object of class \code{nifti}.
#' @param e2 is an object of class \code{nifti}.
#' @author John Muschellli \email{muschellij2@@gmail.com}
#' @examples
#' 
#' img01 <- nifti(array(1:64, c(4,4,4,1)), datatype=4)
#' img02 <- nifti(array(64:1, c(4,4,4,1)), datatype=4)
#' is.nifti(img01 + img02)
#' is.nifti(sqrt(2) * img01)
#' is.nifti(img02 / pi)
#' 
setMethod("Ops", signature(e1="nifti", e2="nifti"),
          function(e1, e2) {
            ## either use drop_img_dim and validObject or take out both
            # e1 = drop_img_dim(e1)            
            # e2 = drop_img_dim(e2)            
            e1@.Data <- callGeneric(e1@.Data, e2@.Data)
            e1 <- resetSlopeIntercept(e1)
            e1 <- calibrateImage(e1, infok=TRUE)
            # e1 = drop_img_dim(e1, warn=FALSE)
            ### creating the datatype that is the maximal precision
            new.dtype <- max(datatype(e1), datatype(e2))
            datatype(e1) <- new.dtype
            bitpix(e1) <- convert.bitpix()[[convert.datatype(new.dtype)]]
            # validObject(e1)
            return(e1)
          }
)
#' @rdname niftiops
#' @aliases Ops,nifti,numeric-method
setMethod("Ops", signature(e1="nifti", e2="numeric"),
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
#' @rdname niftiops
#' @aliases Ops,numeric,nifti-method
setMethod("Ops", signature(e1="numeric", e2="nifti"),
          function(e1, e2) {
            # e2 = drop_img_dim(e2)
            e2@.Data <- callGeneric(e1, e2@.Data)
            e1 <-  e2
            e1 <- resetSlopeIntercept(e1)            
            e1 <- calibrateImage(e1, infok=TRUE)
            # e1 = drop_img_dim(e1)
            # validObject(e1)
            return(e1)
          }
)

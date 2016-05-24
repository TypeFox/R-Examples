.labels <- function (object, which = bquote(), drop = TRUE, ..., use.colnames = TRUE){
  validObject (object)                  

  if (! missing (which) & ! is.character (which))
    warning ("labels are not guaranteed to have the same order as columns.",
             "Consider using indexing by name.")

  ## fill in colnames?
  if (use.colnames){                            
    label <- structure (as.list (c(colnames (object@data), ".wavelength")),
                        names = c(colnames (object@data), ".wavelength"))
    label <- modifyList (label, object@label)

    label <- label [which]
  } else {
    label <- object@label [which]
  }
  
  if (drop && length (label) == 1L)
    label <- label [[1]]

  label
}

##' @include hyperspec-package.R
.test (.labels) <- function (){
  .sort <- function (x)
    x [order (names (x))]
  
  checkEquals (.sort (labels (flu, use.colnames = FALSE)), .sort (flu@label))
  tmp <- c (flu@label, file = "file")
  checkEquals (.sort (labels (flu)), .sort (tmp))
  
  tmp <- flu
  tmp@label$file <- NULL
  checkEquals (labels (tmp, c ("c", "file", "fil"), use.colnames = FALSE),
               structure(list(c = "c / (mg / l)", `NA` = NULL, `NA` = NULL),
                         .Names = c("c", NA, NA)))

  checkEquals (labels (tmp, c ("c", "file", "fil")),
               structure(list(c = "c / (mg / l)", file = "file", `NA` = NULL),
                         .Names = c("c", "file", NA)))
}

##' @rdname labels
##' @usage
##' labels (object, which = NULL, ...) <- value
##' 
##' @aliases labels<-,hyperSpec-method
##' @export "labels<-"
##' @param value the new label(s)
##' @return  \code{labels<-} returns a \code{hyperSpec} object.
##' @examples
##' 
##' labels (flu, "c") <- expression ("/" ("c", "mg / l"))
##' 
`labels<-` <- function (object, which = NULL, ..., value){
  chk.hy (object)
  validObject (object)

  if (is.null (which))
    object@label <- value
  else {
    if ((is.character (which) && !which %in% colnames (object@data)) &&
         which != ".wavelength" ||      # neither a colname nor .wavelength
        (is.numeric (which) && (which < 1 || which > ncol (object@data) + 1)) ||
        (is.logical (which) && length (which) != ncol (object@data) + 1)
                                        # or outside valid indices
        )
      stop ("Column to label does not exist!")

    object@label [[which]] <- value

  }

  validObject (object)
  object
}


##' Get and Set Labels of a hyperSpec Object
##' \code{value} may be a list or vector of labels giving the new label for
##' each of the entries specified by \code{which}.
##' 
##' The names of the labels are the same as the colnames of the
##' \code{data.frame}.  The label for the wavelength axis has the name
##' \code{.wavelength}.
##' 
##' The labels should be given in a form ready for the text-drawing functions
##' (see \code{\link[grDevices]{plotmath}}), e.g. as \code{expression} or a
##' \code{character}.
##' 
##' @param object a hyperSpec object
##' @param which numeric or character to specify the label(s)
##' @param ... ignored
##' @param drop if the result would be a list with only one element, should the
##'   element be returned instead?
##' @param use.colnames should missing labels be replaced by column names of
##'   the extra data?
##' @return \code{labels} returns a list of labels.  If \code{drop} is
##'   \code{TRUE} and the list contains only one element, the element is
##'   returned instead.
##' @docType methods
##' @rdname labels
##' @author C. Beleites
##' @seealso \code{\link[base]{labels}}
##' @export
##' @examples
##' 
##' labels (chondro)
##' 
setMethod ("labels", signature = signature (object = "hyperSpec"), .labels)


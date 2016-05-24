.na.if.different <- function (x) {
  if (length (unique (x)) > 1) NA else x[1]
}


.apply <- function (data, MARGIN, FUN, ...){

  if (length (data$spc) == 0)
    stop ("empty spectra matrix.")

  spc <- apply (data [, "spc", drop = FALSE], MARGIN, FUN, ...)

  if (MARGIN == 1){
    if (is.null (spc))
      spc <- matrix (ncol = 0, nrow = nrow (data))
    else if (is.vector (spc))
      dim (spc) <- c(length (spc), 1)
    else if (is.matrix (spc))
      spc <- t (spc)

    data$spc <- I(spc)
  } else if (MARGIN == 2){
    if (is.null (spc))
      return (data [0, ])
    if (is.null (dim (spc)))
      dim (spc) <- c(1, ncol (data$spc))

    if (all(dim (spc) == dim (data$spc))){
      data$spc <- spc
    }  else {
      nrow <- nrow (spc)

      cols <- colnames (data)
      cols <- which (cols != "spc")
      if (length (cols) > 0) {
        data [1,cols] <- lapply (data [,cols,drop = FALSE], .na.if.different)
      }

      data <- data [rep (1, nrow), , drop = FALSE]

      data$spc <- I (spc)
      rownames (data) <- rownames (spc)
    }
  }

  data$spc <- unclass (data$spc)
  data
}



##' apply
##' Computes summary statistics for the spectra of a \code{hyperSpec} object.
##' 
##' \code{apply} gives the functionality of \code{\link[base]{apply}} for
##' \code{hyperSpec} objects.
##' 
##' The generic functions of group \code{\link[methods]{Math}} are not definded
##' for \code{hyperSpec} objects. Instead, \code{apply} can be used. For
##' functions like \code{log} that work on scalars, \code{MARGIN = 1 : 2} gives
##' the appropriate behaviour.
##' 
##' \code{spcapply} does the same as \code{apply} with \code{MARGIN = 1}, but
##' additionally allows to set a new wavelength axis and adjust the labels.
##' 
##' \code{wlapply} does the same as \code{apply} with \code{MARGIN = 2}, but
##' additionally allows to set a new wavelength axis and adjust the labels.
##'
##' @name apply
##' @rdname apply
##' @aliases apply apply,hyperSpec-method
##' @docType methods
##' @param X,spc a \code{hyperSpec} object
##' @param MARGIN The subscript which the function will be applied over.
##' 
##' \code{1} indicates rows (\code{FUN} is applied to each spectrum),
##' 
##' \code{2} indicates columns (\code{FUN} is applied to each wavelength),
##' 
##' \code{1 : 2} indicates that \code{FUN} should be applied to each single
##'   element of the spectra matrix. Note that many basic mathematical
##'   functions are already defined for hyperSpec objects (see
##'   \code{\link{Math}}).
##' 
##' If \code{MARGIN} is missing, the whole spectra matrix is handed to
##'   \code{FUN}, see also the examples.
##' @param FUN function to compute the summary statistics
##' @param ... further arguments passed to \code{FUN}
##' @param label.wl,label.spc new labels for wavelength and spectral intensity
##'   axes
##' @param new.wavelength for \code{MARGIN = 2}: numeric vector or name of the
##'   argument in \dots{} that is to be used (character) as wavelength axis of
##'   the resulting object.
##' @return A \code{hyperSpec} object
##' @author C. Beleites
##' @seealso \code{\link[base]{apply}}, for applying \code{FUN} to subgroups of
##'   the \code{hyperSpec} object: \code{\link[hyperSpec]{aggregate}}.
##' @keywords methods iteration
##' @export
##' @examples
##' 
##' 
##' plotspc (apply (chondro, 2, range))
##' 
##' avgflu <- apply (flu, 1, mean,
##'                  label.spc = expression (bar (I)),
##'                  new.wavelength = mean (wl (flu)))
##' avgflu
##' 
##' flu[[,,405:407]]
##' apply (flu, 1:2, "*", -1)[[,,405:407]]
##' 
##' ## without MARGIN the whole matrix is handed to FUN
##' apply (flu [,,405:407], , print) [[]]
##' 
##' ## whereas MARGIN = 1 : 2 leads to FUN being called for each element separately
##' apply (flu [,,405:407], 1 : 2, print) [[]]
##' 
setMethod ("apply", signature = signature (X = "hyperSpec"),
           function (X, MARGIN, FUN, ...,
                     label.wl = NULL, label.spc = NULL, new.wavelength = NULL){
  validObject (X)

  if (missing (MARGIN)){                # apply for functions that the complete spectra matrix
    ## is easier: tmp <- apply (x, , FUN, ...)
    ## does: 
    ## tmp <- x
    ## tmp [[]] <- FUN (x [[]], ...)

    X@data$spc <- do.call (FUN, list (X@data$spc, ...))
    
  } else if (all (MARGIN == 1 : 2)){    # apply for functions that take scalar arguments. 

    tmp <- apply (X@data$spc, MARGIN = MARGIN, FUN, ...)
    tmp <- as.numeric (tmp)             # otherwise surprises will be waiting

    dim (tmp) <- dim (X@data$spc)       

    X@data$spc <- tmp
    
  } else {
    ## the usual: for each row / for each column

    X@data <- .apply (X@data, MARGIN = MARGIN, FUN = FUN, ...)

    if (all (MARGIN == 1)) {
      
      ## if the number of data points per spectrum is changed, the wavelength vector needs to be
      ## adapted, too
      
      if (ncol (X@data$spc) != length (X@wavelength)) {
        
        ## only internal functions here: the validation will fail until the wavelength axis is
        ## adjusted
        
        if (!is.null (new.wavelength)){    # vector with new wavelength is given
          if (is.numeric (new.wavelength)) # either directly,
            .wl (X) <- new.wavelength
          else {
            dots <- list (...)
            .wl (X) <- dots [[new.wavelength]] # or as name of the argument that becomes the new
                                               # wavelength vector
          }
        } else if (ncol (X@data$spc) != length (X@wavelength)){
          wl <- as.numeric (colnames (X@data$spc)) # if not given, try to make from colnames of the 
                                                   # spectra matrix
        
          if (length (wl) != ncol (X@data$spc) || any (is.na (wl)))
            wl <- seq_len (ncol (X@data$spc)) # or just number sequentially

          .wl (X) <- wl
        }
      }
    }
  }

  if (!is.null (label.wl))
    X@label$.wavelength <- label.wl

  if (!is.null (label.spc))
    X@label$spc <- label.spc

  validObject (X)

  X
})



.test (apply) <- function (){
  ## check whether .na.if.different is working correctly
  rm (flu)
  flu$equal <- 1
  tmp <- apply (flu, 2, mean)$..
  checkEquals (is.na (tmp),
               structure(c (TRUE, TRUE, FALSE),
                         .Dim = c(1L, 3L),
                         .Dimnames = list (NULL, c("file", "c", "equal")))
               )
  checkEquals (tmp$equal, 1)

  ## POSIXct
  flu$ct <- as.POSIXct(Sys.time()) 
  checkEquals (apply (flu, 2, mean)$ct, flu$ct [1])

  ## POSIXlt
  flu$lt <- as.POSIXlt(Sys.time()) 
  checkEquals (apply (flu, 2, mean)$lt, flu$lt [1])

  rm (flu)
}
